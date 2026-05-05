# ECOSTRESS Level 2 Surface Temperature

[![Ubuntu CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-ubuntu.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-ubuntu.yml)
[![macOS CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-macos.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-macos.yml)
[![Windows CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-windows.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-windows.yml)

This is the main repository for the ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) collection 3 level 2 surface temperature data product algorithm.

Glynn C. Hulley (he/him)<br>
[glynn.hulley@jpl.nasa.gov](mailto:glynn.hulley@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 321H

Robert Freepartner (he/him)<br>
Raytheon

Tinh La (he/him)<br>
[tinh.t.la@jpl.nasa.gov](mailto:tinh.t.la@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 321H

[Gregory H. Halverson](https://github.com/gregory-halverson-jpl) (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 321H

## Prerequisites

### mamba

This C package was designed to be deployed on Linux, but has been retrofitted to compile on macOS and Windows as well, using mamba to consistently install cross-platform dependencies. Continuous integration checks for all three platforms have been included with status badges at the top of the README.

Install [miniforge](https://github.com/conda-forge/miniforge) to obtain `mamba` or `micromamba`. Either is supported — the `MAMBA` variable in the root `Makefile` defaults to `mamba` but can be overridden:

```bash
make environment
```

Running `make environment` creates a conda environment named `ECOv003-L2-LSTE` and installs the following packages from `conda-forge`:

- `hdf4`, `hdf5` — HDF I/O libraries
- `libxml2` — XML configuration parsing
- `eccodes` — GRIB/BUFR meteorological data
- `pkg-config` — build-time dependency resolution

> **Note:** There is no `environment.yml` — packages are installed directly by the `Makefile`. `make install` calls `make environment` automatically, so running them separately is optional.

### RTTOV

This software requires the Radiative Transfer for TOVS (RTTOV) radiative transfer model for atmospheric correction. 

> **Caveat:** RTTOV is not open-source, but is free for registered users.

To obtain it:

1. [Register with the NWP SAF](https://nwp-saf.eumetsat.int/site/register/) (or [log in](https://nwp-saf.eumetsat.int/site/login/) if already registered).
2. Add RTTOV to your software preferences, then download **RTTOV v12** from the [RTTOV v12 page](https://nwp-saf.eumetsat.int/site/software/rttov/rttov-v12/). This package uses **RTTOV 12.2.0**, which is no longer supported by the NWP SAF but remains available for download.

#### Compiling the RTTOV forward model

This repository includes a Fortran 90 forward model driver (`src/rttov_ECOSTRESS_fwd.F90`) that must be compiled against the RTTOV v12 Fortran libraries. Compile it according to the RTTOV v12 build instructions to produce the executable `rttov_ECOSTRESS_fwd.exe`.

#### Coefficient file

The ECOSTRESS instrument coefficient file (`OSP/rtcoef_iss_1_ecostres_v7pred.dat`) is already included in this repository. You do not need to download it separately.

#### Configuring runtime paths

RTTOV is invoked as a subprocess at runtime — it is not linked into the `L2_PGE` binary. Before running the PGE, edit `OSP/PgeRunParameters.xml` to set the correct paths for your installation:

```xml
<scalar name="RttovExe">/path/to/rttov_ECOSTRESS_fwd.exe</scalar>
<scalar name="RttovCoef">/path/to/OSP/rtcoef_iss_1_ecostres_v7pred.dat</scalar>
```

`make install` will succeed without RTTOV present, but the PGE will exit with an error at runtime if `RttovExe` is not a valid path.

#### Licensing

`src/rttov_ECOSTRESS_fwd.F90` carries a EUMETSAT/Met Office copyright. Use of this file is subject to the [RTTOV license agreement](https://nwp-saf.eumetsat.int/site/software/rttov/) accepted upon NWP SAF registration.

## Cross-Platform Installation

A make target for generating a mamba environment has been supplied that will install HDF all other dependencies:

```bash
make environment
```

Activate the `ECOv003-L2-LSTE` mamba environment before compiling:

```bash
mamba activate ECOv003-L2-LSTE
```

Once the mamba environment has been activated on Linux, macOS, or Windows, you should be able to install:

```bash
make install
```
