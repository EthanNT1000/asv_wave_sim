#!/usr/bin/env python3
"""
Gazebo Wave Simulation — Performance Checker

Measures:
  • Environment resolution time  (GZ path setup)
  • World load time              (spawn → first WorldStatistics message)
  • Sustained real-time factor   (RTF) over --duration seconds
  • Per-core CPU utilisation     during the sampling window
  • Linux perf: hardware counters + function-level call-graph (--perf)

Standalone: no ROS / ament dependency. Requires only gz_waves_models in the
install tree (i.e. a built gz-waves workspace).

Prerequisites:
    source install/setup.bash
    pip install psutil                              # optional, per-core CPU
    sudo apt-get install linux-tools-generic        # for perf

Usage:
    python3 src/asv_wave_sim/gz-waves/test/performance/perf_check.py
    python3 src/asv_wave_sim/gz-waves/test/performance/perf_check.py --perf --duration 60
    python3 src/asv_wave_sim/gz-waves/test/performance/perf_check.py --speedup 3.0 --threads 8

If perf_event_paranoid > 2, run once to unlock:
    echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid
"""

import argparse
import ctypes
import os
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path

# ── Optional dependencies ──────────────────────────────────────────────────────

try:
    import psutil
    _HAVE_PSUTIL = True
except ImportError:
    _HAVE_PSUTIL = False

try:
    from gz.msgs10.world_stats_pb2 import WorldStatistics
    from gz.transport13 import Node as GzNode
    _HAVE_GZ = True
except ImportError:
    _HAVE_GZ = False

# ── Constants ──────────────────────────────────────────────────────────────────

WORLD_NAME      = 'perf_waves'
STATS_TOPIC     = f'/world/{WORLD_NAME}/stats'
GZ_LOAD_TIMEOUT = 90.0
GZ_STOP_TIMEOUT = 10.0

# perf settings
PERF_FREQ        = 99     # Hz — prime avoids aliasing with 100 Hz kernel tick
PERF_DWARF_SIZE  = 65528  # bytes — kernel hard limit (65536 - 8 header bytes)
PERF_EVENTS = (
    'cycles,instructions,'
    'cache-references,cache-misses,'
    'branch-instructions,branch-misses,'
    'context-switches,task-clock'
)
PERF_TOP_N  = 25    # number of top functions to show in the report

# ANSI colours
_G = '\033[92m'
_Y = '\033[93m'
_R = '\033[91m'
_B = '\033[94m'
_N = '\033[0m'


# ── Helpers ────────────────────────────────────────────────────────────────────

def _have_lib(name: str) -> bool:
    try:
        ctypes.CDLL(name)
        return True
    except OSError:
        return False


def _perf_paranoid() -> int:
    try:
        return int(Path('/proc/sys/kernel/perf_event_paranoid').read_text().strip())
    except Exception:
        return 4


def _ptrace_scope() -> int:
    try:
        return int(Path('/proc/sys/kernel/yama/ptrace_scope').read_text().strip())
    except Exception:
        return -1   # not present (non-Linux or Yama not loaded)


def _resolve_world_sdf() -> Path:
    """Return the bundled perf_waves.sdf next to this script."""
    return Path(__file__).resolve().parent / 'worlds' / 'perf_waves.sdf'


def _find_install_root() -> Path | None:
    """Walk up from this script to find the workspace install/ directory."""
    for parent in Path(__file__).resolve().parents:
        install = parent / 'install'
        if install.is_dir():
            return install
    return None


def _resolve_gz_waves_models() -> Path | None:
    """Return .../install/<pkg>/share/gz_waves_models/models, or None."""
    install = _find_install_root()
    if install is None:
        return None
    for candidate in install.glob('*/share/gz_waves_models/models'):
        if candidate.is_dir():
            return candidate
    return None


def _resolve_plugin_dirs() -> list[str]:
    install = _find_install_root()
    if install is None:
        return []
    return sorted(str(p) for p in install.glob('*/lib') if p.is_dir())


def _build_world(src: Path, tmp: Path, speedup: float) -> Path:
    text = src.read_text()
    text = re.sub(
        r'(<real_time_factor>)\s*[^<]*\s*(</real_time_factor>)',
        lambda m: f'{m.group(1)}{speedup:.6g}{m.group(2)}',
        text,
    )
    out = tmp / 'perf_waves.sdf'
    out.write_text(text)
    return out


def _build_gz_env(gz_models: Path | None, threads: int,
                  plugin_dirs: list[str]) -> dict:
    env = dict(os.environ)
    env['GZ_PARTITION']    = 'gz_perf_check'
    env['OMP_NUM_THREADS'] = str(threads)

    resource_paths = []
    if gz_models:
        resource_paths.append(str(gz_models))
    existing_res = env.get('GZ_SIM_RESOURCE_PATH', '')
    if existing_res:
        resource_paths.append(existing_res)
    if resource_paths:
        env['GZ_SIM_RESOURCE_PATH'] = ':'.join(resource_paths)

    if plugin_dirs:
        existing_plug = env.get('GZ_SIM_SYSTEM_PLUGIN_PATH', '')
        env['GZ_SIM_SYSTEM_PLUGIN_PATH'] = (
            ':'.join(plugin_dirs) + (':' + existing_plug if existing_plug else ''))

    return env


def _bar(value: float, lo: float, hi: float, width: int = 20) -> str:
    ratio  = max(0.0, min(1.0, (value - lo) / (hi - lo) if hi > lo else 0.0))
    filled = int(ratio * width)
    return '[' + '#' * filled + '-' * (width - filled) + ']'


# ── Checker ────────────────────────────────────────────────────────────────────

class PerfChecker:

    def __init__(self, world: Path | None, speedup: float,
                 duration: int, threads: int,
                 use_perf: bool, perf_out: Path) -> None:
        self._world_src = world
        self._speedup   = speedup
        self._duration  = duration
        self._threads   = threads
        self._use_perf  = use_perf
        self._perf_dir  = perf_out

        self._proc:         subprocess.Popen | None = None
        self._tmp:          tempfile.TemporaryDirectory | None = None
        self._gz_node:      GzNode | None = None
        self._ready         = threading.Event()
        self._perf_record:  subprocess.Popen | None = None
        self._perf_stat:    subprocess.Popen | None = None

        # Collected results
        self.env_time:      float        = 0.0
        self.load_time:     float        = 0.0
        self.load_ok:       bool         = False
        self.rtf_samples:   list[float]  = []
        self.cpu_samples:   list[list[float]] = []
        self.plugin_dirs:   list[str]    = []
        self.gz_models:     Path | None  = None
        self.world_path:    Path | None  = None
        self.have_fftw_omp: bool         = False
        self.have_libomp:   bool         = False
        self.have_perf:     bool         = False
        self.perf_paranoid: int          = 4
        self.ptrace_scope:  int          = -1
        self.perf_ok:       bool         = False
        self.perf_stat_text: str         = ''
        self.perf_top_lines: list[str]   = []

    # ── Phase 1: environment ──────────────────────────────────────────────────

    def check_environment(self) -> None:
        print(f'\n{_B}[1/4] Environment{_N}')
        t0 = time.monotonic()

        self.world_path    = self._world_src or _resolve_world_sdf()
        self.gz_models     = _resolve_gz_waves_models()
        self.plugin_dirs   = _resolve_plugin_dirs()
        self.have_fftw_omp = _have_lib('libfftw3_omp.so')
        self.have_libomp = any(_have_lib(n) for n in (
            'libgomp.so.1',
            'libomp.so.5',
            'libomp.so',
            'libomp.dylib',
        ))
        self.have_perf     = bool(shutil.which('perf'))
        self.perf_paranoid = _perf_paranoid()
        self.ptrace_scope  = _ptrace_scope()

        self.env_time = time.monotonic() - t0

        omp_env = os.environ.get('OMP_NUM_THREADS', f'{self._threads} (will set)')

        world_str = str(self.world_path) if self.world_path else _R + 'NOT FOUND' + _N
        models_str = str(self.gz_models) if self.gz_models else _Y + 'NOT FOUND (model:// may fail)' + _N

        print(f'  World SDF        : {world_str}')
        print(f'  gz_waves_models  : {models_str}')
        print(f'  Plugin dirs      : {len(self.plugin_dirs)} found')
        print(f'  CPU cores        : {os.cpu_count()}')
        print(f'  OMP_NUM_THREADS  : {omp_env}')
        print(f'  libfftw3_omp     : {_G+"yes"+_N if self.have_fftw_omp else _Y+"no"+_N}'
              f'  (multi-threaded FFT)')
        print(f'  libomp           : {_G+"yes"+_N if self.have_libomp else _Y+"no"+_N}'
              f'  (OpenMP runtime)')
        print(f'  gz.transport     : {_G+"yes"+_N if _HAVE_GZ else _Y+"no — RTF disabled"+_N}')
        print(f'  psutil           : {_G+"yes"+_N if _HAVE_PSUTIL else _Y+"no — CPU chart disabled"+_N}')

        if self._use_perf:
            p_colour = _G if self.have_perf else _R
            print(f'  linux perf       : {p_colour + ("yes" if self.have_perf else "NOT FOUND") + _N}')
            par_colour = _G if self.perf_paranoid <= 1 else (_Y if self.perf_paranoid <= 2 else _R)
            print(f'  perf_paranoid    : {par_colour}{self.perf_paranoid}{_N}'
                  + ('' if self.perf_paranoid <= 1
                     else '  ← echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid'))
            if self.ptrace_scope >= 0:
                pt_colour = _G if self.ptrace_scope == 0 else _R
                pt_note = ('' if self.ptrace_scope == 0
                           else '  ← echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope')
                print(f'  ptrace_scope     : {pt_colour}{self.ptrace_scope}{_N}{pt_note}')

        print(f'  Resolution time  : {self.env_time * 1000:.1f} ms')

    # ── Phase 2: launch ───────────────────────────────────────────────────────

    def start_sim(self) -> None:
        print(f'\n{_B}[2/4] World load{_N}')

        if self.world_path is None or not self.world_path.exists():
            print(f'  {_R}Cannot start — world SDF not found: {self.world_path}{_N}')
            return

        self._tmp = tempfile.TemporaryDirectory(prefix='gz_perf_')
        tmp = Path(self._tmp.name)

        world_sdf = _build_world(self.world_path, tmp, self._speedup)
        gz_env    = _build_gz_env(self.gz_models, self._threads, self.plugin_dirs)

        if _HAVE_GZ:
            os.environ['GZ_PARTITION'] = 'gz_perf_check'
            self._gz_node = GzNode()
            self._gz_node.subscribe(WorldStatistics, STATS_TOPIC, self._on_stats)

        print(f'  World SDF  : {world_sdf}')
        print(f'  Speedup    : {self._speedup}×')
        print(f'  Threads    : OMP_NUM_THREADS={self._threads}')
        print('  Spawning headless Gazebo  (gz sim -s -r) ...')

        t0 = time.monotonic()
        self._proc = subprocess.Popen(
            ['gz', 'sim', '-s', '-v4', '-r', str(world_sdf)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            preexec_fn=os.setsid,
            env=gz_env,
        )
        print(f'  PID        : {self._proc.pid}')

        if _HAVE_GZ:
            print(f'  Waiting for /world/{WORLD_NAME}/stats '
                  f'(timeout {GZ_LOAD_TIMEOUT:.0f} s) ...', flush=True)
            self.load_ok   = self._ready.wait(timeout=GZ_LOAD_TIMEOUT)
            self.load_time = time.monotonic() - t0
        else:
            time.sleep(15.0)
            self.load_time = 15.0
            self.load_ok   = self._proc.poll() is None

        status = f'{_G}OK{_N}' if self.load_ok else f'{_R}TIMEOUT{_N}'
        print(f'  Load time  : {self.load_time:.2f} s  [{status}]')

    # ── Phase 3: sample ───────────────────────────────────────────────────────

    def sample(self) -> None:
        print(f'\n{_B}[3/4] Sampling ({self._duration} s){_N}')

        if not self.load_ok:
            print(f'  {_Y}Skipped — world did not load.{_N}')
            return

        if self._use_perf and self.have_perf and self._proc:
            self._start_linux_perf(self._proc.pid)

        if _HAVE_PSUTIL:
            psutil.cpu_percent(percpu=True)

        deadline = time.monotonic() + self._duration
        tick     = 0

        while time.monotonic() < deadline:
            time.sleep(1.0)
            tick += 1

            if _HAVE_PSUTIL:
                self.cpu_samples.append(psutil.cpu_percent(percpu=True))

            rtf_str = (f'RTF={self.rtf_samples[-1]:.3f}x'
                       if self.rtf_samples else 'RTF=...')
            eta = max(0, int(deadline - time.monotonic()))
            print(f'  [{tick:3d}/{self._duration}s]  {rtf_str}  eta={eta}s   ',
                  end='\r', flush=True)

        print()

        if self._use_perf:
            self._stop_linux_perf()
            self._parse_linux_perf()

    # ── Phase 4: stop ─────────────────────────────────────────────────────────

    def stop_sim(self) -> None:
        print(f'\n{_B}[4/4] Shutdown{_N}')
        if self._proc is None:
            print('  Nothing to stop.')
            return

        try:
            os.killpg(os.getpgid(self._proc.pid), signal.SIGTERM)
        except Exception:
            try:
                self._proc.terminate()
            except Exception:
                pass

        try:
            self._proc.wait(timeout=GZ_STOP_TIMEOUT)
            print('  Gazebo stopped cleanly.')
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(self._proc.pid), signal.SIGKILL)
            except Exception:
                self._proc.kill()
            self._proc.wait()
            print('  Gazebo force-killed.')

        self._proc = None
        if self._tmp:
            self._tmp.cleanup()
            self._tmp = None

    # ── Report ────────────────────────────────────────────────────────────────

    def report(self) -> None:
        W   = 62
        sep = '─' * W

        print(f'\n{sep}')
        print('  Gazebo Wave Sim — Performance Report')
        print(sep)

        print(f'\n  Environment')
        print(f'    Resolution time  : {self.env_time * 1000:6.1f} ms')
        print(f'    CPU cores        : {os.cpu_count()}')
        print(f'    OMP_NUM_THREADS  : {self._threads}')
        print(f'    libfftw3_omp     : {"yes" if self.have_fftw_omp else "NO"}')
        print(f'    libomp           : {"yes" if self.have_libomp else "NO"}')
        print(f'    Plugin dirs      : {len(self.plugin_dirs)}')
        print(f'    World SDF        : {self.world_path or "not found"}')

        print(f'\n  World Load')
        c = _G if self.load_ok else _R
        print(f'    Status           : {c}{"OK" if self.load_ok else "TIMEOUT / FAILED"}{_N}')
        print(f'    Load time        : {self.load_time:.2f} s')

        print(f'\n  Real-time Factor  (target {self._speedup}×)')
        if self.rtf_samples:
            avg = sum(self.rtf_samples) / len(self.rtf_samples)
            pct = avg / self._speedup * 100.0 if self._speedup else 0.0
            c   = _G if pct >= 90 else (_Y if pct >= 70 else _R)
            print(f'    Samples          : {len(self.rtf_samples)}')
            print(f'    Average RTF      : {c}{avg:.4f}{_N}  '
                  f'{_bar(avg, 0.0, self._speedup)}')
            print(f'    Min / Max RTF    : {min(self.rtf_samples):.4f} / '
                  f'{max(self.rtf_samples):.4f}')
            print(f'    vs target        : {c}{pct:.1f}%{_N}')
        elif not _HAVE_GZ:
            print(f'    {_Y}N/A — gz.transport13 not found{_N}')
        else:
            print(f'    {_Y}No samples received.{_N}')

        print(f'\n  CPU Utilisation  ({len(self.cpu_samples)} samples @ 1 Hz)')
        if self.cpu_samples:
            n_cores = len(self.cpu_samples[0])
            avgs = [
                sum(s[i] for s in self.cpu_samples) / len(self.cpu_samples)
                for i in range(n_cores)
            ]
            for i, avg_pct in enumerate(avgs):
                c = _G if avg_pct >= 60 else (_Y if avg_pct >= 30 else _N)
                print(f'    Core {i:<2d}          : {c}{avg_pct:5.1f}%{_N}  '
                      f'{_bar(avg_pct, 0.0, 100.0)}')
            all_vals = [v for s in self.cpu_samples for v in s]
            print(f'    Overall avg      : {sum(all_vals)/len(all_vals):.1f}%')
        elif not _HAVE_PSUTIL:
            print(f'    {_Y}N/A — pip install psutil{_N}')
        else:
            print(f'    {_Y}No samples (world did not load).{_N}')

        if self._use_perf:
            self._report_linux_perf()

        print(f'\n  Recommendations')
        self._advice()

        print(f'\n{sep}\n')

    # ── Linux perf report ─────────────────────────────────────────────────────

    def _start_linux_perf(self, pid: int) -> None:
        self._perf_dir.mkdir(parents=True, exist_ok=True)
        stat_out = self._perf_dir / 'stat.txt'
        data_out = self._perf_dir / 'perf.data'

        print(f'  Attaching linux perf to PID {pid} ...')

        if self.ptrace_scope > 0:
            print(f'  {_R}ptrace_scope={self.ptrace_scope}: perf cannot attach to a '
                  f'different session.{_N}')
            print('  Fix: echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope')
            return

        self._perf_stat = subprocess.Popen(
            ['perf', 'stat',
             '-p', str(pid), '--inherit',
             '-e', PERF_EVENTS,
             '-I', '1000',
             '-o', str(stat_out)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        self._perf_record = subprocess.Popen(
            ['perf', 'record',
             f'-F{PERF_FREQ}',
             '-p', str(pid), '--inherit',
             '--call-graph', f'dwarf,{PERF_DWARF_SIZE}',
             '-o', str(data_out)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

        time.sleep(2.0)
        if self._perf_record.poll() is not None:
            err = self._perf_record.stderr.read().decode(errors='replace').strip()
            print(f'  {_R}perf record exited early.{_N}')
            if err:
                for line in err.splitlines():
                    print(f'    perf: {line}')
            if 'ptrace' in err.lower() or 'permission' in err.lower():
                print('  Root cause: ptrace_scope blocks cross-session attachment.')
                print('  Fix: echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope')
            elif 'paranoid' in err.lower():
                print('  Fix: echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid')
            else:
                print('  Try running the script with sudo for full profiling access.')
            self._perf_record = None
            return

        self.perf_ok = True
        print(f'  perf record running  (output → {self._perf_dir})')
        print(f'  perf stat  running   (events: cycles, instructions, cache-misses, …)')

    def _stop_linux_perf(self) -> None:
        for proc in (self._perf_record, self._perf_stat):
            if proc and proc.poll() is None:
                try:
                    proc.send_signal(signal.SIGINT)
                except ProcessLookupError:
                    pass

        for proc in (self._perf_record, self._perf_stat):
            if proc:
                try:
                    proc.wait(timeout=15.0)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()

        self._perf_record = self._perf_stat = None

    def _parse_linux_perf(self) -> None:
        stat_path = self._perf_dir / 'stat.txt'
        data_path = self._perf_dir / 'perf.data'

        if stat_path.exists():
            raw = stat_path.read_text()
            summary_lines = []
            in_summary = False
            for line in raw.splitlines():
                if re.match(r'^\s*#', line) and not in_summary:
                    continue
                if 'Performance counter stats' in line:
                    in_summary = True
                if in_summary:
                    summary_lines.append(line)
            self.perf_stat_text = '\n'.join(summary_lines)

        if data_path.exists() and data_path.stat().st_size > 0:
            try:
                result = subprocess.run(
                    ['perf', 'report',
                     '--stdio',
                     '--sort=symbol,dso',
                     '-n',
                     '--no-children',
                     '-i', str(data_path)],
                    capture_output=True, text=True, timeout=60,
                )
                for line in result.stdout.splitlines():
                    if re.match(r'\s+\d+\.\d+%', line):
                        self.perf_top_lines.append(line.strip())
                        if len(self.perf_top_lines) >= PERF_TOP_N:
                            break
            except Exception as exc:
                print(f'  {_Y}perf report failed: {exc}{_N}')

    def _report_linux_perf(self) -> None:
        print(f'\n  Linux perf  (data → {self._perf_dir})')

        if not self.have_perf:
            print(f'    {_R}perf not found.{_N}')
            print('    sudo apt-get install linux-tools-common linux-tools-generic')
            return

        par = self.perf_paranoid
        if par > 2:
            print(f'    {_Y}perf_event_paranoid={par} — profiling requires elevated privileges.{_N}')
            print('      echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid')

        pt = self.ptrace_scope
        if pt > 0:
            print(f'    {_R}ptrace_scope={pt} — perf cannot attach to gz sim (different session).{_N}')
            print('    gz sim is spawned with os.setsid(); ptrace_scope=0 is required')
            print('    for a separate process to attach across sessions:')
            print('      echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope')
            print('    Permanent (/etc/sysctl.d/99-perf.conf):')
            print('      kernel.yama.ptrace_scope = 0')

        if not self.perf_ok:
            print(f'    {_Y}No perf data collected (perf record did not attach successfully).{_N}')
            return

        if self.perf_stat_text:
            print(f'\n  perf stat  —  hardware counters')
            for line in self.perf_stat_text.splitlines():
                stripped = line.rstrip()
                if stripped:
                    print(f'    {stripped}')
        else:
            print(f'    {_Y}perf stat output not available.{_N}')

        if self.perf_top_lines:
            print(f'\n  Top {len(self.perf_top_lines)} functions by CPU time'
                  f'  (perf record -F{PERF_FREQ} --call-graph dwarf)')
            print(f'    {"Overhead":<10} {"Samples":<9} Symbol / DSO')
            print(f'    {"─"*8:<10} {"─"*7:<9} {"─"*42}')
            for line in self.perf_top_lines:
                print(f'    {line}')
            print(f'\n  Explore further:')
            print(f'    perf report -i {self._perf_dir}/perf.data')
            print(f'    # Flamegraph (requires flamegraph.pl):')
            print(f'    perf script -i {self._perf_dir}/perf.data \\')
            print(f'      | stackcollapse-perf.pl | flamegraph.pl > flame.svg')
        else:
            print(f'\n    {_Y}No function symbols found in perf.data.{_N}')
            print('    Build with -DCMAKE_BUILD_TYPE=RelWithDebInfo for symbol names.')

    # ── Recommendations ───────────────────────────────────────────────────────

    def _advice(self) -> None:
        issues = False

        if not self.have_fftw_omp:
            issues = True
            print(f'    {_Y}• libfftw3_omp not found — multi-threaded FFT unavailable.{_N}')
            print('         sudo apt-get install libfftw3-dev')
            print('      Then: fftw_init_threads() + fftw_plan_with_nthreads(N) in')
            print('      LinearRandomFFTWaveSimulation.cc, link with -lfftw3_omp.')

        if not self.have_libomp:
            issues = True
            print(f'    {_Y}• OpenMP runtime not found.{_N}')
            print('         sudo apt-get install libomp-dev')

        if self.rtf_samples:
            avg = sum(self.rtf_samples) / len(self.rtf_samples)
            pct = avg / self._speedup * 100.0 if self._speedup else 100.0
            if pct < 70:
                issues = True
                print(f'    {_R}• RTF {pct:.0f}% of target — sim is heavily CPU-bound.{_N}')
                print('      Add #pragma omp parallel for collapse(2) to OceanTile.cc ~line 814')
                print('      and call fftw_init_threads() in LinearRandomFFTWaveSimulation.cc.')
            elif pct < 90:
                issues = True
                print(f'    {_Y}• RTF {pct:.0f}% of target — moderate CPU pressure.{_N}')
                print('      OpenMP on OceanTile vertex loop + aero drag reduction should help.')
            else:
                print(f'    {_G}• RTF {avg:.3f}x — healthy (≥90% of {self._speedup}× target).{_N}')

        if self.cpu_samples:
            n  = len(self.cpu_samples[0])
            avgs = [
                sum(s[i] for s in self.cpu_samples) / len(self.cpu_samples)
                for i in range(n)
            ]
            idle   = sum(1 for a in avgs if a < 20.0)
            spread = max(avgs) - min(avgs)
            if idle > 0:
                issues = True
                print(f'    {_Y}• {idle} core(s) idle (<20%) — load not distributed.{_N}')
                print('      OceanTile.cc ~line 814:  #pragma omp parallel for collapse(2)')
            if spread > 40.0:
                issues = True
                print(f'    {_Y}• Core imbalance {spread:.0f}% — try schedule(dynamic) on hydrodynamics loops.{_N}')

        if self._threads < (os.cpu_count() or 1):
            issues = True
            print(f'    {_Y}• {self._threads} thread(s) vs {os.cpu_count()} cores available.{_N}')
            print(f'      Re-run: --threads {os.cpu_count()}')

        if not issues:
            print(f'    {_G}• No issues detected.{_N}')

    # ── gz-transport callback ─────────────────────────────────────────────────

    def _on_stats(self, msg: 'WorldStatistics') -> None:
        self.rtf_samples.append(msg.real_time_factor)
        self._ready.set()


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Measure Gazebo wave simulation performance.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--duration', type=int, default=30,
                        help='Sampling window in seconds (default: 30)')
    parser.add_argument('--speedup', type=float, default=1.0,
                        help='Gazebo real_time_factor target (default: 1.0)')
    parser.add_argument('--threads', type=int, default=os.cpu_count() or 1,
                        help='OMP_NUM_THREADS for Gazebo (default: all cores)')
    parser.add_argument('--world', type=str, default=None,
                        help='Override world SDF path')

    perf_grp = parser.add_mutually_exclusive_group()
    perf_grp.add_argument('--perf', dest='perf', action='store_true',
                          default=False,
                          help='Enable linux perf stat + call-graph recording')
    perf_grp.add_argument('--no-perf', dest='perf', action='store_false',
                          help='Disable linux perf (default)')
    parser.add_argument('--perf-out', type=str, default='./perf_out',
                        help='Directory for perf data files (default: ./perf_out)')

    args = parser.parse_args()

    checker = PerfChecker(
        world    = Path(args.world) if args.world else None,
        speedup  = args.speedup,
        duration = args.duration,
        threads  = args.threads,
        use_perf = args.perf,
        perf_out = Path(args.perf_out),
    )

    try:
        checker.check_environment()
        checker.start_sim()
        checker.sample()
    except KeyboardInterrupt:
        print('\n  Interrupted.')
    finally:
        checker.stop_sim()

    checker.report()


if __name__ == '__main__':
    main()
