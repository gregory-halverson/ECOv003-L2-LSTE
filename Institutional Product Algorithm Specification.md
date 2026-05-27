# Institutional Product Algorithm Specification Document: Land Surface Temperature & Emissivity (LSTE) Processing Pipeline

## Document Metadata
* **Core Software Pipeline Version:** 3.0.3 (ECOSTRESS / NPP SIPS Framework)
* **Classification:** Open Science Processing Software Specification
* **Document Version:** 2.0 (Consolidated Operational Architecture)
* **Document Purpose:** Language-independent, deterministic blueprint combining abstract physical formulations with low-level implementation details to ensure exact binary replication across heterogeneous computing environments.

---

## 1. Pipeline Architecture & High-Level Control Flow

The LSTE pipeline operates as a modular data processing network, executing sequential multi-spectral matrix transformations over observed top-of-atmosphere (TOA) thermal infrared radiances to isolate thermodynamic skin components from atmospheric attenuation.

### 1.1 Structural Module Decomposition
A functional implementation must isolate functional components into distinct programmatic modules to ensure testing boundary containment:
* **Config Loader Module:** Decodes runtime process parameters and runtime properties from XML frameworks.
* **L1 Reader Ingestion Interface:** Parses raw unscaled radiance tracks and sensor look-geometries across varying data layouts (HDF4, HDF5, NetCDF).
* **NWP Adapter Normalization Block:** Adapts non-uniform weather models, processes constraints, and implements water path fallback integrals.
* **RTTOV Interface Engine:** Formats arrays into linear binary streams and executes clear-sky simulations via specialized shell commands.
* **LUT Service Table Resolver:** Handles coordinate linear interpolations over non-uniform sensor calibration data arrays.
* **TG/WVS Modulating Scaler:** Resolves horizontal atmospheric water vapor paths and isolates surface-leaving intensities.
* **TES Kernel Separation Engine:** Executes iterative Planck corrections and minimum-emissivity curve fitting.
* **Cloud Integration & Flagging Matrix:** Dynamically updates 16-bit uncertainty tracking variables and processes scanner data anomalies.
* **SST Split-Window Block:** Processes localized regressions over open-water locations.
* **Product Writer Packing Node:** Quantizes floating-point outputs into standardized structural formats.

### 1.2 Core Data and State Tensors

| Tensor Variable | Scientific Meaning | Memory Layout / Dimensions |
| :--- | :--- | :--- |
| `Y[b, r, c]` | Observed TOA Radiance from L1 product | `[n_channels, n_lines, n_pixels]` |
| `t1r[b, r, c]` | Atmospheric Transmittance from unscaled Pass 1 | `[n_channels, n_lines, n_pixels]` |
| `t2r[b, r, c]` | Atmospheric Transmittance from scaled Pass 2 | `[n_channels, n_lines, n_pixels]` (Optional) |
| `pathr[b, r, c]` | Upwelling Atmospheric Path Radiance | `[n_channels, n_lines, n_pixels]` |
| `skyr[b, r, c]` | Downwelling Reflected Sky Radiance | `[n_channels, n_lines, n_pixels]` |
| `pwv[r, c]` | Mapped Total Column Water Vapor | `[n_lines, n_pixels]` |
| `surfradi[b, r, c]` | Surface-leaving Radiance after atmospheric correction | `[n_channels, n_lines, n_pixels]` |
| `Tg[b, r, c]` | TG/WVS Brightness Temperature Surrogate | `[n_channels, n_lines, n_pixels]` |
| `g[b, r, c]` | Unfiltered Moisture Scaling Gamma Factor | `[n_channels, n_lines, n_pixels]` |
| `gi[b, r, c]` | Smoothed Gamma Filter Field | `[n_channels, n_lines, n_pixels]` |
| `Ts[r, c]` | Evaluated Land Surface Temperature (LST) | `[n_lines, n_pixels]` |
| `emisf[b, r, c]` | Decoupled Narrowband Surface Emissivity Spectra | `[n_channels, n_lines, n_pixels]` |
| `QC[r, c]` | Unified 16-Bit Quality Control Flags Field | `[n_lines, n_pixels]` |

---

## 2. Mathematical Infrastructure and Primitive Utility Functions

### 2.1 Piecewise Linear 1D Interpolation over Non-Uniform Spacing (`interp1d_npts`)
Used extensively for physical mapping inversions (e.g., inverting Planck functions or mapping sensor radiances), this function evaluates an independent query target against non-uniformly spaced lookup vector coordinates using a deterministic binary search tree.

Given an array of independent sample points $X = [x_0, x_1, \dots, x_{N-1}]$ sorted in strict monotonic order, a matching array of dependent coordinates $Y = [y_0, y_1, \dots, y_{N-1}]$, and a target query coordinate $x_q$:

1. **Orientation Check:** If $x_{N-1} < x_0$, the lookup framework is *descending*. If $x_{N-1} > x_0$, it is *ascending*.
2. **Boundary Clamping:**
   * For **ascending** vectors: If $x_q \le x_1$, set the localized baseline index $k = 0$. If $x_q \ge x_{N-2}$, set $k = N - 2$.
   * For **descending** vectors: If $x_q \ge x_1$, set the localized baseline index $k = 0$. If $x_q \le x_{N-2}$, set $k = N - 2$.
3. **Binary Search Traversal:** For any query falling within the remaining interior domain, execute a strict binary search to locate the boundary interval index $k$ such that:
   * $\text{Ascending: } x_k \le x_q \le x_{k+1}$
   * $\text{Descending: } x_k \ge x_q \ge x_{k+1}$
4. **Linear Evaluation:** Compute the segment slope ($m$) and intercept ($b$) to calculate the final interpolated scalar result $y_q$:
$$m = \frac{y_{k+1} - y_k}{x_{k+1} - x_k}$$
$$b = y_{k+1} - m \cdot x_{k+1}$$
$$y_q = m \cdot x_q + b$$

### 2.2 Streaming 2D Mean Filter Optimization (`smooth2d`)
Spatial filtering of multi-channel moisture scale parameters ($\gamma$) prevents high-frequency pixel noise from creating artifact boundaries in final data layers. Standard sliding-window implementations scale linearly with window geometry dimensions, causing computation bottlenecks. This module optimizes performance by tracking running row and column updates to constrain runtime complexity to $O(\text{nrows} \cdot \text{ncols} \cdot 2)$.

For a target matrix $M$ of dimensions $R \times C$, with horizontal and vertical filtering radii designated as $N_r$ and $N_c$, the boundary evaluation window expands to $(2N_r + 1) \times (2N_c + 1)$. The logic must strictly follow these constraints:
1. **NaN Isolation:** Any entry initialized to a Not-a-Number flag ($NaN$) must be omitted from all localized moving summary pools.
2. **State Protection:** If an active index $M(r, c)$ evaluates to $NaN$, or if the localized coordinate window contains exclusively $NaN$ values, the pixel assignment must directly preserve the $NaN$ state into the output matrix.
3. **Streaming Vector Updates:** Track column structures using an array of `SmoothingCol` structures maintaining a running `total` and active `count`. Shift execution windows sequentially by adding the leading boundary row/column elements and subtracting trailing edge parameters.

### 2.3 Spherical Distance Tracking (`distance`)
Distance queries between coordinates utilize a double-precision Haversine equation to account for planetary curvature:
$$\Delta \phi = \text{lat}_2 - \text{lat}_1, \quad \Delta \lambda = \text{lon}_2 - \text{lon}_1$$
$$a = \sin^2\left(\frac{\Delta \phi}{2}\right) + \cos(\text{lat}_1) \cdot \cos(\text{lat}_2) \cdot \sin^2\left(\frac{\Delta \lambda}{2}\right)$$
$$d = 2 \cdot R_{\text{earth}} \cdot \text{atan2}\left(\sqrt{a}, \sqrt{1 - a}\right)$$
Where the planetary radius constant is defined as $R_{\text{earth}} = 6371229.0\text{ m}$.

---

## 3. Global Science Primitives & Constant Coefficients

### 3.1 Numerical and Physics Constraints
* **Equivalence Tolerance ($\epsilon$):** $1.0 \times 10^{-9}$ (Threshold for absolute floating-point evaluations)
* **Temperature Hard Bounds ($T_{min}, T_{max}$):** $150.0 \text{ K}, 380.0 \text{ K}$ (Boundary cutoffs for validity)
* **Vegetation Max Emissivity Base ($e_{max, veg}$):** $0.985$
* **Bare Rock/Soil Max Emissivity Base ($e_{max, bare}$):** $0.970$
* **WVS Exponent Modifiers ($g_1, g_2$):** $1.0, 0.7$
* **Planck Constant Vectors:** $c_1 = 3.7418 \times 10^{-22}\text{ W}\cdot\text{m}^2$, $c_2 = 0.014388\text{ m}\cdot\text{K}$
* **Band Model Calibration Vector ($W_{bmp}$):**
  $$W_{bmp} = [1.2818, 1.5693, 1.6595, 1.8217, 1.8031]$$

### 3.2 Regression Parameter Arrays
* **Wide-Band Emissivity Linear Regression Weights ($\alpha_{wb}$):**
  $$\alpha_{wb} = [0.0715, 0.0657, 0.1970, 0.3384, 0.3703, -0.0443]$$
* **Emissivity Uncertainty Error Grid Coefficients ($X_e$ matrix of dimensions $5 \times 3$):**
  $$X_e = \begin{bmatrix} 
  0.0153 &  0.0155 & -0.0018 \\
  0.0121 &  0.0055 &  0.0004 \\
  0.0128 &  0.0036 &  0.0009 \\
  0.0110 &  0.0017 &  0.0005 \\
  0.0114 & -0.0038 &  0.0025 
  \end{bmatrix}$$
* **Land Surface Temperature Error Estimation Weights ($X_t$ vector):**
  $$X_t = [0.3842, 0.5307, 0.0055]$$

### 3.3 Dynamic Spectral Band Conventions
The pipeline dynamically configures internal matrix pointers based on the evaluated instrument configuration metadata:

| Attribute Dimension | 5-Channel Mode Alignment | 3-Channel Mode Alignment |
| :--- | :--- | :--- |
| **Active Bands** | Bands 1, 2, 3, 4, 5 | Bands 2, 4, 5 |
| **Processing Vector `band[]`** | `[0, 1, 2, 3, 4]` | `[1, 3, 4]` |
| **Reverse Mapping `band_index[]`**| `[0, 1, 2, 3, 4]` | `[-1, 0, -1, 1, 2]` |
| **Reference Thermal Index** | Index 3 (Band 4) | Index 1 (Band 4) |

---

## 4. Comprehensive Stage-by-Stage Operational Specification

### Stage 1: Runtime Configuration And Runtime Parameter Loading
1. Parse the primary XML run configuration file provided via the execution command argument.
2. Ingest `PgeRunParameters.xml` from the configured Operational Support Product (OSP) location.
3. Validate string compliance: Terminate data processing immediately if `PGEVersion` inside the parameter file does not exactly evaluate to compiled constant `PGE_VERSION` (`"3.0.3"`).
4. Parse system context properties: Isolate path handles for `NWP_DIR`, `L2_OSP_DIR`, `ProductPath`, and execution version values (`ProductCounter`). Initialize automated metadata descriptors.
5. Capture system environment signatures using Unix pipe utilities (`date -u` for processing timestamps, and `uname -a` to flag the physical processing cluster architecture).
6. Resolve deterministic product destination file structures derived from orbital telemetry indices (`OrbitNumber`, `SceneID`), version counters, and timestamps embedded in the file headers.

### Stage 2: Spectral Band Set Selection
1. Inspect the input radiance layer header parameters via custom file interfaces.
2. Evaluate active channel bounds: If the file specifies any architecture layout other than 3 or 5 operational channels, drop down to the 3-band configuration mode.
3. Configure the active band lookup parameters (`band` and `band_index` maps) and assign target forward simulation script locations to align with the active sensor channel widths.
4. Define the temperature-driving index configuration string: `lst_band_index = band_index[BAND_4]`. Band 4 serves as the constant reference track for temperature derivations.

### Stage 3: L1 Radiance Ingest and Spatial Grid Boundary Mapping
1. Initialize the internal core sensor metadata structure (`RAD`).
2. Read raw multi-spectral channel radiances into the matrix layers. For system legacy configurations flagged as Collection 2, parse geometry coordinates, tracking views, and terrain parameters out of a separate geolocation dataset (`L1B_GEO`).
3. Construct the central double-precision radiance data tensor $Y(b, r, c)$ spanning active channel rows $r$ and columns $c$.
4. Clean anomalous elements: Scan the array structure, and convert any pixel tracking a radiance value below $0.0 \text{ W}\cdot\text{m}^{-2}\cdot\text{sr}^{-1}\cdot\mu\text{m}^{-1}$ to a standard double-precision Not-a-Number flag ($NaN$).
5. Geolocation Frame Scanning: Interrogate the spatial lat/lon arrays to track extreme limits ($\text{minLat}, \text{maxLat}, \text{minLon}, \text{maxLon}$). This bounding coordinate envelope establishes the geographic clipping template for the subsequent NWP ingestion step.

### Stage 4: Optional ASTER GED Land Surface Baseline Referencing
1. Check process conditions: If Water Vapor Scaling logic is globally active (`run_tgwvs = true`), pass the scene bounding boundaries and coordinate matrices into the auxiliary ASTER data loader module (`read_aster_ged`).
2. Map surface climatology properties: Project the $1^\circ \times 1^\circ$ absolute horizontal land grids to match the target scene layout. Assign ocean pixels a fixed water-surface emissivity constant ($0.990$) and replace land-surface $NaN$ grid gaps with a standard nominal background value ($0.962$).
3. Output the combined continuous 2D surface matrix layout to layer variable `emis_aster`. If `run_tgwvs` is flagged as false, omit the subgrid queries entirely to conserve system execution paths.

### Stage 5: Ingestion & Normalization of NWP Meteorological Profiles
1. Identify the meteorological file structure type by inspecting the path string contained in `NWP_DIR` for key identifiers (`MERRA`, `GEOS`, `NCEP`, or `ECMWF`).
2. **MERRA Processing:** Extract multi-layer absolute variables directly. Apply physical data boundaries to vertical tracking elements. Invert the raw pressure level indexing order to conform to top-to-bottom requirements of the forward radiative model.
3. **GEOS Processing:** Read the pre-cropped spatial weather variables. Record absolute input file pathways inside the global metadata tracking parameters and invert the pressure levels to match forward transfer formats.
4. **NCEP Processing:** Parse the standard spatial meteorology matrix. Double the horizontal grid resolution profile by running multi-set bilinear spatial interpolation functions over the data layers, preserving the vertical pressure indexing order.
5. **Humidity Matrix Inversion:** If the weather archive provides moisture profiles $Q_{profile}$ formatted as Relative Humidity ($RH$), translate the values into an absolute mass mixing ratio ($w_{mr}$) over liquid supercooled water via the Murphy and Koop (2005) vapor pressure formulation before modeling:
$$\ln(e_s) = 54.842763 - \frac{6763.22}{T} - 4.210 \ln(T) + 0.000367 T + \tanh(0.0415(T - 218.8)) \cdot \left(53.878 - \frac{1331.22}{T} - 9.44523 \ln(T) + 0.014025 T\right)$$
$$e = \frac{RH}{100.0} \cdot e_s$$
$$w_{mr} = \left(\frac{e}{P_{levels} - e}\right) \cdot \frac{18.0152}{28.9644}$$
6. **Total Column Water Fallback Integration:** If the weather data format completely lacks explicit precipitable water tracks, execute numerical trapezoidal integration across the vertical layers to resolve total columnar content:
$$\text{TCW} = \frac{1.0}{100.0 \cdot 9.8} \sum_{z=0}^{Z_{\text{surf}}-2} \left[ \frac{q_z + q_{z+1}}{2.0} \cdot \left(\frac{18.0152}{28.9644 \cdot 1000.0}\right) \cdot (P_{z+1} - P_z) \right]$$
Constrain integration loops to execute exclusively on indices where localized pressures map below the true surface pressure limit.

### Stage 6: Structural Forward Model Input Synthesis
1. Construct 2D coordinate meshgrids from the discrete atmospheric lat/lon vectors.
2. Clip the weather fields to match the active scene domain. For non-GEOS datasets, constrain bounds to the geographic scene extrema plus an extra spatial padding margin of $\pm 2.0^\circ$ in both latitude and longitude.
3. Slice the weather profile states into cropped data blocks (`cropT`, `cropQ`, `cropSP`, `cropTCW`).
4. Project satellite look geometries: Run nearest-neighbor coordinate mapping from the fine satellite tracks onto the coarse weather subgrid. Extract localized surface terrain elevations (`cropHsurf`) and view zeniths (`cropSatZen`).
5. Select the surface skin temperature tracker ($T_{\text{skin}}$): Ingest the explicit weather surface skin temperature layer if present; otherwise, substitute the data array from the lowest vertical atmospheric temperature profile level.
6. Initialize the background surface emissivity placeholder matrix (`Bemis`) as a static array filled entirely with the constant value $1.0 \times 10^{-6}$ across all dimensions.

### Stage 7: Profile Stream Flattening and Validation
1. Reshape the cropped weather data structures into continuous, unrolled arrays configured for the binary interface parameters of the forward radiative model script via the `set_rttov_atmos` routine.
2. **Execute Validation and Data Repair Rules:**
   * If a temperature value drops below zero or equals the data error code ($0.0001 \text{ K}$), override the entry using the global spatial mean calculated across valid profile segments.
   * If a skin temperature element drops below the absolute low boundary threshold ($90.0 \text{ K}$), overwrite the entry using the valid grid temperature mean.
   * If near-surface moisture elements report missing values, pull data from the lowest valid positive level in the vertical humidity profile.
3. Apply a modulo mapping step ($\text{lon} \leftarrow \text{fmod}(\text{lon} + 360.0, 360.0)$) to transform longitudes into a $[0^\circ, 360^\circ)$ coordinate system. This transformation must execute *exclusively* during the initial binary profile stream construction pass.
4. Export the validated data stream blocks to disk as a raw binary file (`prof_in.bin`).

### Stage 8: Forward Simulation Execution Loop
The system executes the forward radiative transfer engine script using a dual-pass simulation loop:

1. **Pass 1 Simulation (`wvs_case = 0`):** Run the simulation against the unscaled, baseline binary profile `prof_in.bin`. Save the resulting forward model calculations to disk as a text dataset (`rad_out.dat`).
2. **Pass 2 Simulation (`wvs_case = 1`):** If Water Vapor Scaling laboratory adjustments are active, re-scale the humidity vectors across all dimensions by a factor of 0.7:
$$Q_{\text{scaled}} = Q_{\text{profile}} \cdot 0.7, \quad Qt_{\text{scaled}} = Qt_{\text{profile}} \cdot 0.7, \quad Q2_{\text{scaled}} = Q2_{\text{profile}} \cdot 0.7$$
Export this humidity-suppressed dataset to `prof_in.bin` and execute the simulation script again using the same configuration parameters to generate the updated `rad_out.dat` text dataset.

### Stage 9: Swath Geolocation Remapping & Atmospheric Lookup Synthesis
1. Parse the text data streams (`rad_out.dat`) generated during the simulation passes.
2. Convert internal wavenumber results into pure micrometer radiance tracking fields ($\text{mW}\cdot\text{m}^{-2}\cdot\text{sr}^{-1}\cdot\mu\text{m}^{-1}$) by running continuous piecewise linear interpolations over the radcon conversion table.
3. Map the coarse simulation grid outputs onto the fine satellite swath tracks by running multi-set bilinear interpolation functions (`multi_interp2`).
4. Isolate the pixel-level spatial matrices: $t_{1r}$ (Pass 1 transmittance), $t_{2r}$ (Pass 2 transmittance), $path_r$ (upwelling atmospheric path radiance), $sky_r$ (downwelling sky radiance), and the scene column water vapor grid (`pwv`). If the calculated transmittance array evaluates to non-positive elements across all cells, flag the track as unusable and terminate processing immediately.

### Stage 10: Radiative Calibration Lookup Ingestion
1. Read the static multi-spectral radiance-to-temperature calibration file (`rad_lut_file`) into memory.
2. Enforce constraint checks: Stop execution if the lookup matrix contains fewer than two distinct row records, as this violates boundary conditions for numeric interpolation stability.
3. Map columns to structural parameters: Index 0 tracks temperature scale boundaries ($K$), and Indices 1 through 5 track multi-spectral radiance profiles mapped for sensor bands 1 to 5. The matrix is maintained in memory as an indexed array to optimize continuous piecewise search requests within the primary solver loops.

### Stage 11: Atmospheric and Moisture Attenuation Correction

#### 11A. Standard Atmospheric Inversion Path (No WVS)
If WVS processing flags are disabled, calculate the surface-leaving radiances ($surfradi$) directly from the Pass 1 forward model parameters across all valid coordinates:
$$surfradi(b, r, c) = \frac{Y(b, r, c) - path_r(b, r, c)}{t_{1r}(b, r, c)}$$

#### 11B. Water Vapor Scaling (WVS) Correction Path
If WVS operations are active, resolve moisture variations across the spatial grids using an empirical dual-pass optimization loop:

1. Invoke the `tg_wvs` sub-module to calculate surrogate surface brightness temperatures ($T_g$) using observed radiances $Y$, mapped tracking angles, and the baseline ASTER-GED array layer.
2. Map the extracted $T_g$ fields into equivalent blackbody-leaving radiances ($B$) via linear lookup interpolation:
$$B(b, r, c) = \text{interp1d\_npts}(\text{LUT}_{\text{temp}}, \text{LUT}_{\text{rad\_band}[b]}, T_g(b, r, c))$$
3. Compute the structural moisture parameters for each individual band channel:
$$g_f = g_2^{\text{bmp}[\text{band}[b]]}$$
$$\text{term}_1 = \frac{t_{2r}}{t_{1r}^{g_f}}, \quad \text{term}_{2t} = \frac{B - \left( \frac{path_r}{1.0 - t_{1r}} \right)}{Y - \left( \frac{path_r}{1.0 - t_{1r}} \right)}, \quad \text{term}_3 = \frac{t_{2r}}{t_{1r}}$$
If any pixel maps to a non-positive or undefined $NaN$ state for variable $\text{term}_{2t}$ within any band, invalidate the entire coordinate column index across all processing dimensions.
4. Evaluate the raw uncalibrated moisture parameter scaling matrix ($g$):
$$\text{term}_2 = \text{term}_{2t}^{(g_1 - g_f)}$$
$$g = \frac{\ln(\text{term}_1 \cdot \text{term}_2)}{\ln(\text{term}_3)}$$
If the logarithmic operations yield an undefined, imaginary, or complex scalar result, substitute a static value of $1.0$.
5. Because the internal land surface type mask defaults to land cover (`gp_water = 1`), apply a greybody structural tracking step that overrides the scaling values across all bands using the data from Channel 5 (the primary moisture validation track):
$$g(b, r, c) = g(\text{Index 4}, r, c)$$
6. Convert cloudy pixels to $NaN$, clamp the remaining valid land grid parameters to the range $[-2.0, 3.0]$, and smooth the matrix by passing it through the optimized `smooth2d` utility to compile the smoothed scaling array $g_i$.
7. Interrogate the land-cover average: Calculate the average of Channel 4 $T_g$ values across all land coordinates. If this spatial mean drops below threshold limit `TGthresh` ($290.0\text{ K}$), select the conservative correction constraints (`smooth_scale1 = 750`, `PWVthresh1 = 2.5`). Otherwise, apply the standard processing configurations (`smooth_scale2 = 150`, `PWVthresh2 = 2.0`).
8. If a pixel reports a low total column water vapor value ($PWV < \text{PWV}_{\text{thresh}}$), force the scaling parameters to a value of $1.0$ to prevent correction oversampling artifacts.
9. Reconstruct the blended atmospheric transmission ($t_i$) and path radiance ($path_i$) parameters to output the true surface-leaving radiance matrices ($surfradi$):
$$t_i(b, r, c) = [t_{1r}(b, r, c)]^{\frac{g_i(b, r, c) - gf}{1.0 - gf}} \cdot [t_{2r}(b, r, c)]^{\frac{1.0 - g_i(b, r, c)}{1.0 - gf}}$$
$$path_i(b, r, c) = path_r(b, r, c) \cdot \left( \frac{1.0 - t_i(b, r, c)}{1.0 - t_{1r}(b, r, c)} \right)$$
$$surfradi(b, r, c) = \frac{Y(b, r, c) - path_i(b, r, c)}{t_i(b, r, c)}$$
Any element evaluating to a non-positive surface radiance value is converted to a standard $NaN$ flag.

### Stage 12: Core TES Temperature/Emissivity Separation Retrieval
The core separation solver (`apply_tes_algorithm`) evaluates the spatial surface-leaving radiance matrices to resolve individual emissivity spectra and isolate true surface skin temperatures. The algorithm iterates through each pixel footprint independently using the following modules:

#### 12.1 Module 1: Normalized Emissivity Method (NEM) Core
1. Inspect the surface type mask to set the initial maximum emissivity constraint ($e_{max}$). Because the surface water mask defaults to land cover across all indices, initialize calculations using the exposed rock/soil constraint: $e_{max} = 0.970$.
2. Calculate the initial surface emission intensity by subtracting the downwelling sky reflection contribution:
$$R_b(b) = surfradi(b) - (1.0 - e_{max}) \cdot sky_r(b)$$
3. Convert the radiance values into an initial brightness temperature spectrum $T(b)$ by performing piecewise linear interpolations over the sensor's calibration table:
$$T(b) = \text{interp1d\_npts}(\text{LUT}_{\text{rad\_band}[b]}, \text{LUT}_{\text{temp}}, \frac{R_b(b)}{e_{max}})$$
4. Extract the baseline target temperature from the warmest evaluated channel: $T_{nem} = \max(T(b))$.
5. **Iterative Relaxation Loop:** Execute refinement loops to optimize the temperature calculation against an iteration limit of 13.
   * Interpolate expected blackbody emission intensities: $B(b, T_{nem}) = \text{interp1d\_npts}(\text{LUT}_{\text{temp}}, \text{LUT}_{\text{rad\_band}[b]}, T_{nem})$.
   * Update the spectral emissivity estimate: $e(b) = \frac{R_b(b)}{B(b, T_{nem})}$.
   * Update the corrected surface emission intensity: $R_{new}(b) = surfradi(b) - (1.0 - e(b)) \cdot sky_r(b)$.
   * If the convergence delta satisfies $\Delta R = |R_{new}(b) - R_b(b)| < 0.05$ across all bands and at least 3 loops have occurred, the loop terminates. Otherwise, the target temperature is updated using the updated emission spectrum ($T_{nem} = \max(\text{interp1d\_npts}(\text{LUT}_{\text{rad\_band}[b]}, \text{LUT}_{\text{temp}}, R_{new}/e))$) and the cycle repeats. If the loop reaches 13 iterations without stabilizing, mark the pixel as a processing failure.

#### 12.2 Module 2: Min-Max Difference (MMD) Operational Calibration
1. Normalize the raw emissivity components into scaling ratio parameters $\beta(b)$ to remove calibration offsets:
$$\beta(b) = \frac{e(b)}{\frac{1}{3}\sum_{i=\text{Index 1}}^{\text{Index 3}} e(i)}$$
2. Isolate the maximum spectral contrast range ($MMD_2$) across the operational bands:
$$MMD_2 = \max(\beta(b)) - \min(\beta(b))$$
3. Evaluate the empirical power-law regression equations using the bare soil calibration constants ($co = [0.9950, 0.7264, 0.8002]$) to calculate the minimum absolute surface emissivity floor value ($e_{min2}$):
$$e_{min2} = co[0] - co[1] \cdot MMD_2^{co[2]}$$
4. Re-scale the relative $\beta$ spectra to compute the final calibrated surface emissivities ($emisf$):
$$emisf(b) = \beta(b) \cdot \left( \frac{e_{min2}}{\min(\beta(b))} \right)$$
5. The pipeline selects Channel 4 as the clearest spectral channel. It calculates the final Land Surface Temperature variable ($T_s$) by running an inverse lookup interpolation over the decoupled thermal radiance field:
$$R_{final} = \frac{surfradi(b_{clear}) - (1.0 - emisf(b_{clear})) \cdot sky_r(b_{clear})}{e_{final}(b_{clear})}$$
$$T_s = \text{interp1d\_npts}(\text{LUT}_{\text{rad\_band}[b_{clear}]}, \text{LUT}_{\text{temp}}, R_{final})$$
If the estimated emissivity for the temperature-driving channel drops below 0.0, both $T_s$ and the emissivity array elements for that pixel are set to 0.0.

### Stage 13: Cloud Mask Integration & Summary Statistics Creation
1. The pipeline delegates cloud processing to the external module `process_cloud`. It passes observed radiances, geometry metadata, collection numbers, and the preliminary TES emissivity arrays to execute multi-threshold cloud tests.
2. The cloud screening module converts radiances to brightness temperatures, evaluates regional meteorological limits, and exports the final binary cloud mask matrix (`Cloud_final`) to file.
3. The main processing loop re-opens the cloud dataset and extracts the spatial tracking layer to update quality control arrays and metadata attributes.
4. **Compile Cloud Summary Statistics:** For pixels flagged as cloud-covered, calculate an approximate cloud-top temperature proxy by applying an environmental lapse-rate adjustment over surface terrain heights ($\text{BT}_{11} = \text{Base\_Threshold} - \text{Height} \cdot \text{slapse}$). Record the final statistical parameters (cloud cover percentage; mean, minimum, maximum, and standard deviation of the cloud temperature proxy) into the global product metadata structures.

### Stage 14: Wideband Emissivity Synthesis
The pipeline calculates the integrated wideband emissivity product ($EmisWB$) by executing a linear combination of the final narrowband emissivity layers:
$$EmisWB = c_0 + \sum_{b=0}^{N-1} c_b \cdot emisf(b)$$
The regression weights ($c_b$) and the constant offset parameter ($c_0 = \alpha_{wb}[N]$) are loaded dynamically from the configuration files to match either the 3-band or 5-band sensor mode.

### Stage 15: Error Estimation And QA Flag Refinement
The system calculates absolute uncertainty parameters across each pixel using empirical polynomial regressions:

#### 15.1 Parametric Uncertainty Calibration
* **Narrowband Emissivity Uncertainty ($d\epsilon_b$):**
$$d\epsilon_b = X_e[b][0] + X_e[b][1] \cdot \text{TCW} + X_e[b][2] \cdot \text{TCW}^2$$
* **Surface Temperature Uncertainty ($dT$):**
$$dT = X_t[0] + X_t[1] \cdot \text{TCW} + X_t[2] \cdot \text{SVA}$$
Where $\text{TCW}$ represents the local total column water vapor value from the `pwv` matrix, and $\text{SVA}$ represents the sensor zenith look angle from `vRAD.Satzen`. Calculate the total emissivity root-mean-square error ($RMSE_{\epsilon}$):
$$RMSE_{\epsilon} = \sqrt{\frac{1}{N} \sum_{b=0}^{N-1} (d\epsilon_b)^2}$$

#### 15.2 Incomplete Input Scan Profiling
The system scans the input radiance line quality codes (`DataQ`). If a band records a missing scan line that was filled using neural-network spatial inpainting loops (`DataQ = 1`), the pipeline inflates the local uncertainty metrics to reflect the reduced data confidence:
$$dT \leftarrow dT + 0.34\text{ K}$$
$$d\epsilon_b \leftarrow d\epsilon_b + \text{Emis\_ScanErr}[b]$
Where $\text{Emis\_ScanErr} = [0.0099, 0.0053, 0.0050, 0.0047, 0.0023]$. If any band records an un-fillable, corrupt, or missing data line (`DataQ` value of 2, 3, or 4), the pixel is marked as unproduced, and the data layers are forced to the standard fill state.

#### 15.3 Core Validation Gatekeeper Constraints
* If $T_s < 100.0\text{ K}$, flag the calculation as anomalous, force the data layers to $NaN$, and update the QC mask to an unproduced state (`0x0F`).
* If $T_s > 380.0\text{ K}$, retain the nominal evaluation but modify the QC mask to flag the pixel as suspicious.
* If $T_s$ evaluates to $NaN$ while the QC mask reports a valid retrieval, override the flag and force it to an unproduced state.
* For any pixel marked as unproduced, ensure that both $T_s$ and all narrowband emissivity parameters are forced to $NaN$.

### Stage 16: Sea Surface Temperature (SST) Processing
The pipeline calculates Sea Surface Temperature ($SST$) values independently from the land core routines across all geolocated pixels:

1. The pipeline isolates the current scene temporal properties to open the matching monthly/hourly NetCDF coefficient file (`ECOSTRESS_SSTv3_Coeffs_MM_HH.nc`) corresponding to the nearest 6-hour interval.
2. The coarse coefficient structures are cropped and bilinearly interpolated onto the satellite track coordinates to produce the continuous local coefficient arrays $x_{eco1}$, $x_{eco2}$, $x_{eco3}$, and $x_{eco4}$.
3. Brightness temperatures $TB_4$ and $TB_5$ are resolved by passing the raw channel radiances through the lookup functions.
4. The pipeline evaluates the non-linear split-window regression equation to output the sea surface temperature field ($T_{sea}$):
$$\sec(\theta) = \frac{1.0}{\cos\left(\theta \cdot \frac{\pi}{180.0}\right)}$$
$$T_{sea} = x_{eco1} + x_{eco2} \cdot TB_4 + x_{eco3} \cdot (TB_4 - TB_5) + x_{eco4} \cdot (1.0 - \sec(\theta)) \cdot (TB_4 - TB_5)$

### Stage 17: Quantization and Output Data Packing
The completed double-precision floating-point arrays are compressed into fixed-point integer products before file serialization. The quantization step executes a standard rounding truncation mapping:
$$\text{Packed\_Int} = \text{floor}\left( \frac{\text{Double\_Value} - \text{Offset}}{\text{Scale}} + 0.5 \right)$$

The data packing must follow the specifications listed in the table below:

| HDF5 Product Variable Path | Output DataType | Scale Factor | Base Offset | Operational Packing Domain |
| :--- | :--- | :--- | :--- | :--- |
| `/Data Fields/LST` | `UInt16` | $0.02$ | $0.0$ | $7500 \dots 65535\text{ K}$ |
| `/Data Fields/SST` | `UInt16` | $0.02$ | $0.0$ | $7500 \dots 65535\text{ K}$ |
| `/Data Fields/LST_err` | `UInt8` | $0.04$ | $0.0$ | $1 \dots 255\text{ K}$ |
| `/Data Fields/Emisb` | `UInt8` | $0.002$ | $0.49$ | $1 \dots 255$ (Dimensionless) |
| `/Data Fields/Emisb_err` | `UInt16` | $1.0 \times 10^{-4}$ | $0.0$ | $1 \dots 65535$ (Dimensionless) |
| `/Data Fields/PWV` | `UInt16` | $0.001$ | $0.0$ | $1 \dots 65535\text{ cm}$ |

Pixels flagged as unproduced, missing, or obscured by clouds are assigned a value of `0` within the packed integer structures.

### Stage 18: HDF-EOS Structure Serialization
1. Construct the core HDF-EOS group layout hierarchy under path `/HDFEOS/GRIDS/ECO_L2G_LSTE_70m`.
2. Write all quantized output integer datasets and floating-point geometry layers.
3. If running in 3-band mode, populate dummy data arrays for bands 1 and 3 to ensure schema stability across products.
4. Export the global file metadata attributes, appending processing system profiles (`AncillaryNWP`, `CollectionLabel`, `BuildID`, and bounding box coordinates).
5. Generate the descriptive external CAS XML metadata files (`.met`) for both the LSTE and Cloud data products to facilitate automated catalog ingest operations.

---

## 5. Unified Quality Control Bit-Mapping Architecture

The Quality Control layer maps processing metadata into a spatial matrix of 16-bit unsigned integer bitmasks (`uint16`). The bit configurations are defined below: