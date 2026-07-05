# Wave Sim

[![Ubuntu Jammy CI](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ubuntu-jammy-ci.yml/badge.svg)](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ubuntu-jammy-ci.yml)
[![macOS Ventura CI](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/macos13-ventura-ci.yml/badge.svg)](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/macos13-ventura-ci.yml)
[![Cpplint](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ccplint.yml/badge.svg)](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ccplint.yml)
[![Cppcheck](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ccpcheck.yml/badge.svg)](https://github.com/srmainwaring/asv_wave_sim/actions/workflows/ccpcheck.yml)

This package contains plugins that support the simulation of waves and surface vessels in [Gazebo](https://gazebosim.org/home).

![rs750_ardupilot_v3_upwind](https://user-images.githubusercontent.com/24916364/228044489-b434b1ae-c30f-4676-9415-1719ee75479b.gif)


The main branch targets [Gazebo Garden](https://gazebosim.org/docs/garden) and no longer has a dependency on ROS. 

There are new features including FFT wave generation methods, ocean tiling, and support for the [Ogre2](https://github.com/OGRECave/ogre-next) render engine. There are some changes in the way that the wave parameters need to be set, but as far possible we have attempted to retain compatibility with the earlier versions. Further details are described below where you can also find a section describing [support for legacy versions of Gazebo](#legacy-versions).

## Dependencies

- A working installation of [Gazebo Garden](https://gazebosim.org/docs/garden) or later including development symbols.

- The simulation uses the [CGAL](https://www.cgal.org/) library for mesh manipulation and [FFTW](http://www.fftw.org/) to compute Fourier transforms. Both libraries are licensed GPL-3.0.

- [OpenMP](https://www.openmp.org/) is used to parallelise wave mesh updates, hydrodynamics force calculations, and FFT execution across multiple CPU cores.

- This fork packages `gz-waves` and `gz-waves-models` as ROS 2 `ament_cmake` packages (see [ROS 2 / ament build](#ros-2--ament-build) below), and extends the hydrodynamics system with a Fossen-style per-DOF damping model, planing-hull foil lift, aerodynamic drag on above-waterline faces, a bulk water-current field, and telemetry/sensor plugins (InfluxDB export, speed-through-water topic, anemometer). See the [Changes](#changes) section and the physics write-ups in [`doc/`](doc/) for details.

## Ubuntu

- Ubuntu 22.04 (Jammy)
- Gazebo Sim, version 7.1.0 (Garden)

Install CGAL and FFTW:

```zsh
sudo apt-get update
sudo apt-get install libcgal-dev libfftw3-dev
```

### OpenMP (Ubuntu)

OpenMP is included with GCC and requires no extra package. Verify it is available:

```bash
echo '#include <omp.h>
int main() { return omp_get_max_threads(); }' | g++ -fopenmp -x c++ - -o /tmp/omp_check && echo "OpenMP OK"
```

If you are building with **Clang** instead of GCC, install the LLVM OpenMP runtime:

```bash
sudo apt-get install libomp-dev
```

To enable multi-threaded FFTW (optional but recommended for large wave grids), install the OpenMP-enabled FFTW variant:

```bash
sudo apt-get install libfftw3-dev
# libfftw3-dev already includes the threaded library (libfftw3_omp)
# verify with:
ls /usr/lib/x86_64-linux-gnu/libfftw3_omp*
```

## macOS

- macOS 12.6 (Monterey)
- Gazebo Sim, version 7.1.0 (Garden)

Install CGAL and FFTW:

```zsh
brew update
brew install cgal fftw
```

### OpenMP (macOS)

Apple Clang does not ship with OpenMP. Install the LLVM OpenMP runtime via Homebrew:

```zsh
brew install libomp
```

Then tell CMake where to find it by adding these flags to your `colcon build` command:

```zsh
-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
-DOpenMP_CXX_LIB_NAMES="omp" \
-DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib"
```

Full macOS build command with OpenMP:

```zsh
colcon build --symlink-install --merge-install --cmake-args \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DBUILD_TESTING=ON \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_MACOSX_RPATH=FALSE \
  -DCMAKE_INSTALL_NAME_DIR=$(pwd)/install/lib \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib"
```

## Installation

### Create a workspace

```bash
mkdir -p gz_ws/src
```

### Clone and build the package

Clone the `asv_wave_sim` repository:

```bash
cd ~/gz_ws/src
git clone https://github.com/srmainwaring/asv_wave_sim.git
```

Compile the package:

#### Ubuntu

```bash
colcon build --symlink-install --merge-install --cmake-args \
-DCMAKE_BUILD_TYPE=RelWithDebInfo \
-DBUILD_TESTING=ON \
-DCMAKE_CXX_STANDARD=17
```

Source the workspace:

```bash
source ./install/setup.bash
```

#### macOS

```bash
colcon build --symlink-install --merge-install --cmake-args \
-DCMAKE_BUILD_TYPE=RelWithDebInfo \
-DBUILD_TESTING=ON \
-DCMAKE_CXX_STANDARD=17 \
-DCMAKE_MACOSX_RPATH=FALSE \
-DCMAKE_INSTALL_NAME_DIR=$(pwd)/install/lib
```

Source the workspace:

```bash
source ./install/setup.zsh
```

### Build the GUI plugin (optional) 

There is an optional GUI plugin that controls the wave parameters.

```bash
cd ~/gz_ws/src/asv_wave_sim/gz-waves/src/gui/plugins/waves_control 
mkdir build && cd build
cmake .. && make
```

### ROS 2 / ament build

`gz-waves` and `gz-waves-models` are `ament_cmake` packages (see their `package.xml` files) and build as part of a normal ROS 2 / colcon workspace, with no other ROS 2 dependencies required.

Each package registers `ament_environment_hooks` (`gz-waves/hooks/gz-waves1.{dsv,sh}.in` and `gz-waves-models/hooks/gz_waves_models.{dsv,sh}.in`) so that sourcing the workspace's `install/setup.bash` automatically prepends:

- `GZ_SIM_SYSTEM_PLUGIN_PATH` and `GZ_GUI_PLUGIN_PATH` with the installed `gz-waves1` plugin and gui directories.
- `GZ_SIM_RESOURCE_PATH` with the installed `gz_waves_models` `models/` and `worlds/` directories.

This means the manual `export` commands in [Set environment variables](#set-environment-variables) below are only needed if you are not sourcing the workspace via colcon (e.g. building `gz-waves` standalone with plain CMake).

## Usage

### Set environment variables

If your workspace was built and sourced with `colcon` (see [ROS 2 / ament build](#ros-2--ament-build)), `GZ_SIM_RESOURCE_PATH`, `GZ_SIM_SYSTEM_PLUGIN_PATH` and `GZ_GUI_PLUGIN_PATH` are already set by the package's environment hooks and the exports below are not required.

```bash
# for future use - to support multiple Gazebo versions
export GZ_VERSION=garden

# not usually required as should default to localhost address
export GZ_IP=127.0.0.1

# ensure the model and world files are found
export GZ_SIM_RESOURCE_PATH=\
$GZ_SIM_RESOURCE_PATH:\
$HOME/gz_ws/src/asv_wave_sim/gz-waves-models/models:\
$HOME/gz_ws/src/asv_wave_sim/gz-waves-models/world_models:\
$HOME/gz_ws/src/asv_wave_sim/gz-waves-models/worlds

# ensure the system plugins are found
export GZ_SIM_SYSTEM_PLUGIN_PATH=\
$GZ_SIM_SYSTEM_PLUGIN_PATH:\
$HOME/gz_ws/install/lib

# ensure the gui plugin is found
export GZ_GUI_PLUGIN_PATH=\
$GZ_GUI_PLUGIN_PATH:\
$HOME/gz_ws/src/asv_wave_sim/gz-waves/src/gui/plugins/waves_control/build
```

### Ubuntu VM

If running on an Ubuntu virtual machine you may need to use software rendering if the hypervisor does not support hardware acceleration for OpenGL 4.2+. Install `mesa-utils` to enable llvmpipe:

```bash
sudo apt-get install mesa-utils
```

To use the llvmpipe software renderer, prefix Gazebo commands with the `LIBGL_ALWAYS_SOFTWARE` environment variable:

```bash
LIBGL_ALWAYS_SOFTWARE=1 gz sim waves.sdf
```

## Examples

On macOS the client and server must be launched separately. The commands may be combined on Ubuntu.

Launch a Gazebo session.

Server:

```bash
gz sim -v4 -s -r waves.sdf
```

Client:

```bash
gz sim -v4 -g
```

The session should include a wave field and some floating objects.

## Changes

There are some changes to the plugin SDF schema for hydrodynamics and waves.   

### Waves model and visual plugins

- The `filename` and `name` attributes for the wave model and visal plugins have changed.
- The `<size>` element has been renamed to `<tile_size>` and moved into `<waves>`
- The `<cell_count>` element has been moved into `<waves>`
- Add new element `<algorithm>` to specify the wave generation algorithm. Valid options are: `sinusoid`, `trochoid` and `fft`.
- Add new element `<wind_velocity>` for use with the `fft` algorithm.
- Add new element `<wind_speed>` for use with the `fft` algorithm.
- Add new element `<wind_angle_deg>` for use with the `fft` algorithm.

```xml
<plugin
    filename="gz-waves1-waves-model-system"
    name="gz::sim::systems::WavesModel">
    <static>0</static>
    <update_rate>30</update_rate>
    <wave>
      <!-- Grid dimensions
        - The tile_size and cell_count may be a single value
          for square grids, or a 2d vector if different resolution
          is desired along the x and y axis.
        - The cell_count must be a power of 2 for fft waves
      -->
      <!-- Either: single value for square grids -->
      <tile_size>256.0</tile_size>
      <cell_count>128</cell_count>

      <!-- Or: 2d vectors for different resolution in each axis -->
      <tile_size>256.0 64.0</tile_size>
      <cell_count>128 32</cell_count>

      <!-- Wave algorithms
        - These elements specify the wave generation method
          and wave spectrum parameters.
      -->

      <!-- Either: `fft` waves parameters -->
      <algorithm>fft</algorithm>
      <wind_speed>5.0</wind_speed>
      <wind_angle_deg>135</wind_angle_deg>
      <steepness>2</steepness>

      <!-- Or: `trochoid` waves parameters -->
      <algorithm>trochoid</algorithm>
      <number>3</number>
      <scale>1.5</scale>
      <angle>0.4</angle>
      <amplitude>0.4</amplitude>
      <period>8.0</period>
      <phase>0.0</phase>
      <steepness>1.0</steepness>
      <direction>1 0</direction>
    </wave>
</plugin>
```

The waves visual plugin has the same algorithm elements as the model plugin and extra elements to control the shading algorithm. Two approaches are available:

  - `DYNAMIC_GEOMETRY` uses PBS shaders and is suitable for small areas.
  - `DYNAMIC_TEXTURE` uses a custom shader and is suitable for tiled areas.

```xml
<plugin
    filename="gz-waves1-waves-visual-system"
    name="gz::sim::systems::WavesVisual">
  <static>0</static>

  <!-- set the mesh deformation method  -->
  <mesh_deformation_method>DYNAMIC_GEOMETRY</mesh_deformation_method>

  <!-- number of additional tiles along each axis -->
  <tiles_x>-1 1</tiles_x>
  <tiles_y>-1 1</tiles_y>
  <wave>
    <!-- `fft` wave parameters -->
    <algorithm>fft</algorithm>
    <tile_size>256.0</tile_size>
    <cell_count>128</cell_count>
    <wind_speed>5.0</wind_speed>
    <wind_angle_deg>135</wind_angle_deg>
    <steepness>2</steepness>
  </wave>

  <!--
    Shader parameters only apply when using DYNAMIC_TEXTURE
  -->

  <!-- shader program -->
  <shader language="glsl">
    <vertex>materials/waves_vs.glsl</vertex>
    <fragment>materials/waves_fs.glsl</fragment>
  </shader>
  <shader language="metal">
    <vertex>materials/waves_vs.metal</vertex>
    <fragment>materials/waves_fs.metal</fragment>
  </shader>

  <!-- vertex shader params -->
  <param>
    <shader>vertex</shader>
    <name>world_matrix</name>
  </param>
  <param>
    <shader>vertex</shader>
    <name>worldviewproj_matrix</name>
  </param>
  <param>
    <shader>vertex</shader>
    <name>camera_position</name>
  </param>
  <param>
    <shader>vertex</shader>
    <name>rescale</name>
    <value>0.5</value>
    <type>float</type>
  </param>
  <param>
    <shader>vertex</shader>
    <name>bumpScale</name>
    <value>64 64</value>
    <type>float_array</type>
  </param>
  <param>
    <shader>vertex</shader>
    <name>bumpSpeed</name>
    <value>0.01 0.01</value>
    <type>float_array</type>
  </param>
  <param>
    <shader>vertex</shader>
    <name>t</name>
    <value>TIME</value>
  </param>

  <!-- pixel shader params -->
  <param>
    <shader>fragment</shader>
    <name>deepColor</name>
    <value>0.0 0.05 0.2 1.0</value>
    <type>float_array</type>
  </param>
  <param>
    <shader>fragment</shader>
    <name>shallowColor</name>
    <value>0.0 0.1 0.3 1.0</value>
    <type>float_array</type>
  </param>
  <param>
    <shader>fragment</shader>
    <name>fresnelPower</name>
    <value>5.0</value>
    <type>float</type>
  </param>
  <param>
    <shader>fragment</shader>
    <name>hdrMultiplier</name>
    <value>0.4</value>
    <type>float</type>
  </param>
  <param>
    <shader>fragment</shader>
    <name>bumpMap</name>
    <value>materials/wave_normals.dds</value>
    <type>texture</type>
    <arg>0</arg>
  </param>
  <param>
    <shader>fragment</shader>
    <name>cubeMap</name>
    <value>materials/skybox_lowres.dds</value>
    <type>texture_cube</type>
    <arg>1</arg>
  </param>

</plugin>
```

### Hydrodynamics plugin

- The `filename` and `name` attributes for the hydrodynamics plugin have changed.
- The hydrodynamics parameters are now scoped in an additional `<hydrodynamics>` element.
- The buoyancy and hydrodynamics forces can be applied to specific entities
in a model using the `<enable>` element. The parameter should be a fully
scoped model entity (model, link or collision name).
- The `<wave_model>` element is not used.
- The single scalar linear/angular damping coefficients (`cDampL1/L2/R1/R2`) have been
replaced by a full per-DOF [Fossen](doc/Hydrodynamics%20Physics%20Algorithm.md)-style
damping matrix (`cDampU*`/`V*`/`W*`/`P*`/`Q*`/`N*` for surge/sway/heave/roll/pitch/yaw).
A `<randomize>` block can be added instead to draw each coefficient from a uniform
distribution, useful for domain randomization across training episodes.
- Foil lift (planing-hull dynamic lift) and above-waterline aerodynamic drag can be
enabled per-model. The plugin can also load a bulk water-current field, publish the
model's speed through water, and stream per-triangle telemetry to InfluxDB — see
[Telemetry, sensors and water current](#telemetry-sensors-and-water-current) below.

```xml
<plugin
  filename="gz-waves1-hydrodynamics-system"
  name="gz::sim::systems::Hydrodynamics">

  <!-- Apply hydrodynamics to the entire model (default) -->
  <enable>model_name</enable>

  <!-- Or apply hydrodynamics to named links -->
  <enable>model_name::link1</enable>
  <enable>model_name::link2</enable>

  <!-- Or apply hydrodynamics to named collisions -->
  <enable>model_name::link1::collision1</enable>
  <enable>model_name::link1::collision2</enable>

  <!-- Above-waterline aerodynamic drag (flat-plate approximation) -->
  <aerodynamic_drag_on>1</aerodynamic_drag_on>
  <cAeroDrag>1.0</cAeroDrag>

  <!-- Publish speed through water (defaults to /model/<model>/speed_through_water) -->
  <SpeedThroughWater>
    <topic>/model/model_name/speed_through_water</topic>
    <link_name>base_link</link_name>
  </SpeedThroughWater>

  <!-- Stream per-triangle hydrodynamics telemetry to InfluxDB over UDP (line protocol) -->
  <influxDBUdp>
    <ip>127.0.0.1</ip>
    <port>8094</port>
    <measurement>asv_wave_sim_triangle</measurement>
    <update_rate>1000.0</update_rate>
  </influxDBUdp>

  <!-- Hydrodynamics -->
  <hydrodynamics>
    <damping_on>1</damping_on>
    <viscous_drag_on>1</viscous_drag_on>
    <pressure_drag_on>1</pressure_drag_on>

    <!-- Per-DOF Fossen damping: C<DOF><order>, DOF = U/V/W/P/Q/N
         (surge/sway/heave/roll/pitch/yaw), order = 1 linear, 2 quadratic -->
    <cDampU1>1.0E-6</cDampU1>
    <cDampU2>1.0E-6</cDampU2>
    <cDampV1>1.0E-4</cDampV1>
    <cDampV2>1.0E-4</cDampV2>
    <cDampW1>1.0E-4</cDampW1>
    <cDampW2>1.0E-4</cDampW2>
    <cDampP1>1.0E-3</cDampP1>
    <cDampP2>1.0E-3</cDampP2>
    <cDampQ1>1.0E-3</cDampQ1>
    <cDampQ2>1.0E-3</cDampQ2>
    <cDampN1>1.0E-4</cDampN1>
    <cDampN2>1.0E-4</cDampN2>

    <!-- Or, instead of fixed cDamp* values, randomize each coefficient
         (uniformly) once per load, e.g. for domain randomization -->
    <!-- <randomize>1</randomize>
    <dampUMin>1.0e-7</dampUMin><dampUMax>1.0e-5</dampUMax> -->

    <!-- 'Pressure' Drag -->
    <cPDrag1>1.0E+2</cPDrag1>
    <cPDrag2>1.0E+2</cPDrag2>
    <fPDrag>0.4</fPDrag>
    <cSDrag1>1.0E+2</cSDrag1>
    <cSDrag2>1.0E+2</cSDrag2>
    <fSDrag>0.4</fSDrag>
    <vRDrag>1.0</vRDrag>

    <!-- Planing-hull foil lift -->
    <foil_lift_on>1</foil_lift_on>
    <cLift1>1.0</cLift1>
    <alphaStall>0.2618</alphaStall>

    <!-- Optional bulk water-current field, preprocessed from HEC-RAS output
         with preprocess_hecras.py into the binary WCRG v1 grid format -->
    <water_current_grid>/path/to/current_grid.bin</water_current_grid>
  </hydrodynamics>

  <!-- Control visibility of markers -->
  <markers>
    <update_rate>10</update_rate>
    <water_patch>1</water_patch>
    <waterline>1</waterline>
    <underwater_surface>1</underwater_surface>
  </markers>
</plugin>
```

### Telemetry, sensors and water current

- **Water current** — `WaterCurrentGrid` (`gz-waves/include/gz/waves/WaterCurrentGrid.hh`) loads a pre-processed binary grid (format `WCRG v1`, produced from HEC-RAS output by `preprocess_hecras.py`) and bilinearly samples a bulk current velocity at each submerged triangle's centroid, added to the wave orbital velocity to give the total fluid velocity used by drag, damping and foil lift.
- **Speed through water** — the hydrodynamics plugin publishes the hull velocity relative to the fluid (wave orbital + current) as a `gz.msgs.Vector3d` on the topic set by `<SpeedThroughWater><topic>`, defaulting to `/model/<model_name>/speed_through_water`.
- **InfluxDB telemetry** — when `<influxDBUdp>` is configured, per-triangle and per-submerged-triangle properties (position, velocities, forces) are sent as InfluxDB line protocol over UDP at `<update_rate>` Hz, batched up to the maximum UDP packet size for efficiency.
- **Anemometer** — `gz::sim::systems::Anemometer` (`gz-waves/src/systems/anemometer/`), ported from [`srmainwaring/asv_sim`](https://github.com/srmainwaring/asv_sim), is a custom Gazebo sensor that publishes apparent wind speed and direction, configured under a sensor's `<gz:anemometer>` element (supports `<noise>`).
- **Aerodynamic drag** — a flat-plate drag model (`<aerodynamic_drag_on>`, `<cAeroDrag>`) applies wind force/torque to above-waterline triangles, complementing the anemometer and wave/current fields.

For the full derivation of the damping, drag and lift models see [`doc/Hydrodynamics Physics Algorithm.md`](doc/Hydrodynamics%20Physics%20Algorithm.md) and [`doc/Hydrodynamics physics functions.md`](doc/Hydrodynamics%20physics%20functions.md).

## Tests

```bash
# build with tests
$ colcon build --merge-install --cmake-args -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_MACOSX_RPATH=FALSE -DCMAKE_INSTALL_NAME_DIR=$(pwd)/install/lib -DBUILD_TESTING=ON --packages-select gz-waves1

# run tests
colcon test --merge-install 

# check results
colcon test-result --all --verbose 
```

Testing within a project build directory

```bash
$ cd ~/gz_ws/src/asv_wave_sim/gz-waves
$ mkdir build && cd build
$ cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTING=ON
$ make && make test
```

## Plots

Plots may be generated for some of the wave spectra and wave simulation methods:

```bash
./install/bin/PLOT_WaveSpectrum
```

## Legacy versions

There is no plan to back-port new features to Gazebo9 or Gazebo11. The following branches are maintained for legacy support:

- [`gazebo9`](https://github.com/srmainwaring/asv_wave_sim/tree/gazebo9) - for Gazebo9 / ROS Melodic / Ubuntu 18.04 (Bionic).

- [`gazebo11`](https://github.com/srmainwaring/asv_wave_sim/tree/gazebo11) - for Gazebo11 / ROS Noetic / Ubuntu 20.04 (Focal).

In addition there are three branches that contain development iterations of the FFT wave simulation - for Gazebo11 / ROS Noetic / Ubuntu 20.04 (Focal):

- [`feature/fft-waves-v1`](https://github.com/srmainwaring/asv_wave_sim/tree/feature/fft-waves-v1)
- [`feature/fft-waves-v2`](https://github.com/srmainwaring/asv_wave_sim/tree/feature/fft-waves-v2)
- [`feature/fft-waves-v3`](https://github.com/srmainwaring/asv_wave_sim/tree/feature/fft-waves-v3)


## License

This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the [GNU General Public License](LICENSE) for more details.

This project makes use of other open source software, for full details see the file [LICENSE_THIRDPARTY](LICENSE_THIRDPARTY).

## Acknowledgments

- Jacques Kerner's two part blog describing boat physics for games: [Water interaction model for boats in video games](https://www.gamasutra.com/view/news/237528/Water_interaction_model_for_boats_in_video_games.php) and [Water interaction model for boats in video games: Part 2](https://www.gamasutra.com/view/news/263237/Water_interaction_model_for_boats_in_video_games_Part_2.php).
- The [CGAL](https://doc.cgal.org) libraries are used for the wave field and model meshes.
- The [UUV Simulator](https://github.com/uuvsimulator/uuv_simulator) package for the orginal vertex shaders used in the wave field visuals.
- The [VMRC](https://bitbucket.org/osrf/vmrc) package for textures and meshes used in the wave field visuals.
- Jerry Tessendorf's paper on	[Simulating Ocean Water](https://people.cs.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf)
- Curtis Mobley's web book [Ocean Optics](https://www.oceanopticsbook.info/) in particular the section on [Modeling Sea Surfaces](https://www.oceanopticsbook.info/view/surfaces/level-2/modeling-sea-surfaces) and [example IDL code](https://www.oceanopticsbook.info/packages/iws_l2h/conversion/files/IDL-SurfaceGenerationCode.zip)  
