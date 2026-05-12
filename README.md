# ECOSTRESS Level 2 Surface Temperature

[![Ubuntu CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-ubuntu.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-ubuntu.yml)
[![macOS CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-macos.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-macos.yml)
[![Windows CI](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-windows.yml/badge.svg?branch=main)](https://github.com/ECOSTRESS-Collection-3/ECOv003-L2-LSTE/actions/workflows/ci-windows.yml)

This is the main repository for the ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) collection 3 level 2 surface temperature data product algorithm.

Glynn C. Hulley (he/him)<br>
[glynn.hulley@jpl.nasa.gov](mailto:glynn.hulley@jpl.nasa.gov)<br>
NASA Jet Propulsion Laboratory 321H

Robert Freepartner (he/him)<br>
[robert.freepartner@jpl.nasa.gov](robert.freepartner@jpl.nasa.gov)<br>
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

To utilize the cross-platform installation, install conda/mamba using [miniforge](https://github.com/conda-forge/miniforge).

### RTTOV

This software requires installation of Radiative Transfer for TOVS (RTTOV), which can be obtained from [Eumetsat](https://nwp-saf.eumetsat.int/site/software/rttov/).

This radiative transfer model is used for atmospheric correction.

## Cross-Platform Installation

This C package was designed to be deployed on Linux, but has been retrofitted to compile on macOS and Windows as well, using mamba to consistently install cross-platform dependencies. Continuous integration checks for all three platforms have been included with status badges at the top of the README.

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
