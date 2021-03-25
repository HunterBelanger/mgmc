# Multi-Group Monte Carlo Transport Code - MGMC
![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)

MGMC is a 3D, multi-group, Monte Carlo transport code which solves the
Boltzmann neutron transport equation for k-eigenvalue problems. 

All problem parameters such as geometry, materials, and simulation settings are
controlled through a single YAML input file. Several examples of such files are
provided in the `input_file` directory. A problem geometry may be constructed
using combinations of surface half-spaces to create cells. Like in other Monte
Carlo codes, Universes and Lattices may be used to construct more complex
geometries such as for fuel assemblies found in nuclear reactors. Material
cross sections may be provided with any number of energy groups, so long as all
materials in the problem have the same number of groups. Currently, only
isotropic scattering is supported.

Two different particle tracking methods are available to users: delta-tracking
[1], and carter-tracking [2]. Either of these may be selected in the settings
portion of the input file.  When using carter-tracking, if the majorant cross
section is under-estimated, the particle population will not be stable unless
particle weight cancellation is used. A mesh may be defined over which 3D
regional cancellation is used to annihilate positive and negative weight.

Tallies for the scalar flux and power (fission rate) on a rectilinear mesh are
available, both using a collision estimator. The scalar flux and power are
saved in binary Numpy files (`.npy`) to facilitate quick and easy analysis and
plotting in Python. Plots of the geometry may also be specified in the input
file, and are generated when the program is run with the `--plot` flag. Shared
memory parallelism is implemented with OpenMP, and is turned on by default.

[1] J. Leppänen, “On the use of delta-tracking and the collision flux estimator
in the Serpent 2 Monte Carlo particle transport code,” Ann Nucl Energy, vol.
105, pp. 161–167, 2017, doi:
[10.1016/j.anucene.2017.03.006](https://dx.doi.org/10.1016/j.anucene.2017.03.006).

[2] L. L. Carter, E. D. Cashwell, and W. M. Taylor, “Monte Carlo Sampling with
Continuously Varying Cross Sections Along Flight Paths,” Nucl Sci Eng, vol. 48,
no. 4, pp. 403–411, 1972, doi:
[10.13182/nse72-1](https://dx.doi.org/10.13182/nse72-1). 


## License
MGMC is released under the CeCILL-v2.1 license. This is a French open source
copyleft license.  It is very similar to the GPL, and is explicitly compatible
with with the GPL, AGPL, and EUPL licenses. For more information about this
license, please look at the `LICENSE` file for the French version, and the
`LICENSE-ENGLISH` file for the English version. You may also reference the
[CeCILL website](https://cecill.info/).

## Install
To build MGMC, a linux system with a C++17 compliant compiler is required (gcc
7 works), along with cmake >= 3.11. A few third-party compile-time dependencies
([docopt](http://docopt.org/), [pcg](https://www.pcg-random.org),
[ndarray](https://github.com/HunterBelanger/ndarray)) are shipped in the
`vendor/` directory. In addition to these two basic requirements, the
[yaml-cpp](https://github.com/jbeder/yaml-cpp) library with header files is
also needed to build the library. If not already installed, the library will be
downloaded and built as part of the build process. If you prefer, you may
optionally install the library using your distribution's package manager, and it
should be automatically picked up by the build system instead.

On Ubuntu/Debian based distros, yaml-cpp may be installed with the following
command:
```
$ sudo apt install libyaml-cpp-dev
```
For Fedora/CentOS distros:
```
$ sudo dnf install yaml-cpp-devel
```
On Arch based distros:
```
$ sudo pacman -S yaml-cpp
```

To build the MGMC executable, once inside the source directory, run
```
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make -j
```
This will produce an executable called `mgmc`. You can run the provided example
with the following commands:
```
$ cp ../input_files/ref_sqr_c5g7.yaml .
$ ./mgmc -i ref_sqr_c5g7.yaml
```

## Contributors
MGMC was developed by Hunter Belanger in the framework of his Ph.D. thesis at
the French Alternative Energies and Atomic Energy Commission (CEA).
