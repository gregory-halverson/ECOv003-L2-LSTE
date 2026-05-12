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

Dr. Tanvir Islam<br>
NASA Jet Propulsion Laboratory

Dr. Nabin Malakar<br>
NASA Jet Propulsion Laboratory

Simon Latyshev<br>
Raytheon

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

## Algorithm

```mermaid
flowchart TD
	subgraph Inputs
		L1CG[ECOSTRESS L1CG TOA radiance]
		NWP[NWP profiles\n(temp, WV, pressure, skin temp)]
		ASTER[ASTER GED emissivity / NDVI]
		LUTs[Algorithm LUTs\n(Planck, cloud BT, SST coeffs, uncertainty coeffs)]
		GEO[Geolocation / viewing geometry]
		MASK[Land-water / snow masks]
	end

	L1CG --> AC
	NWP --> AC
	GEO --> AC
	AC[Atmospheric Correction\n(RTTOV)] --> WVS
	ASTER --> WVS
	WVS[TG-WVS humidity correction\n(applied when PWV is high)] --> TES

	TES[TES retrieval\n(NEM + MMD + LST)] --> LSTE
	LUTs --> TES
	MASK --> TES

	TES --> CLOUD
	LUTs --> CLOUD
	CLOUD[Cloud detection\n(BT-LUT + emissivity test)] --> CLOUD_OUT

	TES --> SST
	GEO --> SST
	LUTs --> SST
	MASK --> SST
	SST[Ocean split-window SST\n(replaces TES LST over water)] --> LSTE

	TES --> UQ
	GEO --> UQ
	LUTs --> UQ
	UQ[Uncertainty quantification\n(LST/emissivity error + quality tiers)] --> QC

	LSTE[ECOv003_L2G_LSTE\n(LST + emissivity fields)]
	CLOUD_OUT[ECOv003_L2G_CLOUD\n(cloud mask)]
	QC[QC plane\n(error estimates + quality flags)]
```



### Atmospheric Correction (RTTOV)

Top-of-atmosphere radiances from the ECOSTRESS L1CG product are atmospherically corrected using the Radiative Transfer for TOVS (RTTOV) model. NWP atmospheric profiles (temperature, water vapor mixing ratio, surface pressure, and skin temperature) from MERRA-2, GEOS-5, NCEP, or ECMWF are interpolated to each pixel and passed to RTTOV, which computes per-band transmittance (`t`), upwelling path radiance (`pathr`), and downwelling sky radiance (`skyr`).

The standard atmospheric correction is:

$$L_{surf} = \frac{L_{TOA} - L_{path}}{\tau}$$

where $L_{surf}$ is the atmospherically corrected bottom-of-atmosphere surface-leaving radiance in W/m²/sr/μm.

#### Water Vapor Scaling (TG-WVS)

When precipitable water vapor (PWV) exceeds a threshold (~2.0–2.5 cm), an additional Water Vapor Scaling (WVS) correction is applied. ASTER GED emissivity is used to estimate a ground brightness temperature `Tg`, which is then used to derive per-band gamma correction factors (`g`) that scale the effective water vapor path. RTTOV is run a second time with water vapor scaled by 0.7 to produce a corrected atmospheric state, and the transmittance and path radiance are blended between the two runs using the gamma factors. This improves surface radiance retrieval under high-humidity conditions.

### Temperature & Emissivity Separation (TES)

The TES algorithm simultaneously retrieves land surface temperature (LST) and per-band surface emissivity from the atmospherically corrected surface-leaving radiances. It operates entirely in radiance space (W/m²/sr/μm) using pre-computed Planck function lookup tables.

**NEM (Normalized Emissivity Method):** An assumed maximum emissivity (`emax`) is used as a starting point. For each pixel and band, the reflected sky component is removed:

$$R_i = L_{surf,i} - (1 - \epsilon_{max}) \cdot L_{sky,i}$$

The temperature in each band is retrieved by inverting the Planck function via LUT. The maximum band temperature is taken as the NEM temperature estimate `Tnem`.

**MMD (Max-Min Difference):** Per-band emissivities are computed from the ratio of surface radiance to the Planck function evaluated at `Tnem`. The spectral contrast (MMD) is used to refine the emissivity estimate through an empirical relationship between MMD and minimum emissivity. This step is iterated until convergence.

**LST:** The final LST is derived from the clearest (highest-emissivity) band using the refined emissivity and the inverted Planck function. Over water pixels, a fixed emissivity is applied.

### Cloud Detection

Cloud detection is performed by `process_cloud()` after TES retrieval. The algorithm applies two complementary tests:

1. **BT-LUT test:** Band 4 brightness temperature (TB4) is compared against a climatological clear-sky BT threshold from a time-of-day-dependent LUT (four 6-hourly files covering 00, 06, 12, and 18 UTC). Pixels colder than the threshold are flagged as cloud.

2. **BT-difference / emissivity test (Collection 3):** The mean of Bands 4 and 5 emissivity (`Emismm`, spatially smoothed with a 5×5 window) provides an additional discriminator between cloud and low-emissivity surfaces.

Cloud flags are spatially extended by `cloud_extend` pixels and written to the `ECOv003_L2G_CLOUD` product.

### Sea Surface Temperature

Over ocean pixels, a dedicated Sea Surface Temperature (SST) algorithm replaces the TES-derived LST. SST is computed from Band 4 and Band 5 brightness temperatures and the satellite zenith angle using a split-window regression:

$$SST = c_1 \cdot TB_4 + c_2 \cdot (TB_4 - TB_5) + c_3 \cdot (TB_4 - TB_5) \cdot \sec(\theta) + c_4$$

Coefficients are loaded from monthly/6-hourly NetCDF LUT files (`ECOSTRESS_SSTv3_Coeffs_MM_HH.nc`) and bilinearly interpolated to the granule grid using a geolocation reference file.

### Uncertainty Quantification

Per-pixel LST and emissivity error estimates are computed as a function of PWV and satellite zenith angle using empirical polynomial coefficients (`xe`). Error thresholds are tiered into three quality levels (low/medium/high uncertainty), with limits of `[1.0, 1.5, 2.5]` K for LST and `[0.013, 0.015, 0.017]` for emissivity. These errors and the corresponding quality flags are written to the `QC` data plane in the L2G LSTE product.

## References

### Core algorithm and product references

- Hulley, G. C., Göttsche, F. M., Rivera, G., Hook, S. J., Freepartner, R. J., Martin, M. A., Cawse-Nicholson, K., & Johnson, W. R. (2022). Validation and quality assessment of the ECOSTRESS Level-2 land surface temperature and emissivity product. *IEEE Transactions on Geoscience and Remote Sensing, 60*, 1–23. https://doi.org/10.1109/TGRS.2021.3079879

- Gillespie, A., Rokugawa, S., Matsunaga, T., Cothern, J. S., Hook, S., & Kahle, A. B. (1998). A temperature and emissivity separation algorithm for Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) images. *IEEE Transactions on Geoscience and Remote Sensing, 36*(4), 1113–1126. https://doi.org/10.1109/36.700995

- Sabol, D. E., Jr., Gillespie, A. R., Abbott, E., & Yamada, G. (2009). Field validation of the ASTER Temperature–Emissivity Separation algorithm. *Remote Sensing of Environment, 113*(11), 2328–2344. https://doi.org/10.1016/j.rse.2009.06.008

- Hulley, G. C., Hook, S. J., Abbott, E., Malakar, N., Islam, T., & Abrams, M. (2015). The ASTER Global Emissivity Dataset (ASTER GED): Mapping Earth's emissivity at 100 meter spatial scale. *Geophysical Research Letters, 42*(19), 7966–7976. https://doi.org/10.1002/2015GL065564

- Saunders, R., Matricardi, M., & Brunel, P. (1999). An improved fast radiative transfer model for assimilation of satellite radiance observations. *Quarterly Journal of the Royal Meteorological Society, 125*(556), 1407–1425. https://doi.org/10.1002/qj.1999.49712555615

- Saunders, R., Hocking, J., Turner, E., Rayer, P., Rundle, D., Brunel, P., Vidot, J., Roquet, P., Matricardi, M., Geer, A., Bormann, N., & Lupu, C. (2018). An update on the RTTOV fast radiative transfer model (currently at version 12). *Geoscientific Model Development, 11*(7), 2717–2737. https://doi.org/10.5194/gmd-11-2717-2018

- Meng, X., Cheng, J., Yao, B., & Guo, Y. (2022). Validation of the ECOSTRESS land surface temperature product using ground measurements. *IEEE Geoscience and Remote Sensing Letters, 19*, 1–5. https://doi.org/10.1109/LGRS.2021.3123816

