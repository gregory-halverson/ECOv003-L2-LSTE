# ECOL2 User Guide V3

JPL D-103137

ECOsystem Spaceborne Thermal Radiometer
Experiment on Space Station (ECOSTRESS) Mission
Level 2 Product User Guide for Collection 3
12 December, 2026

Glynn Hulley, Tinh La, Robert Freepartner
ECOSTRESS Algorithm Development Team
Jet Propulsion Laboratory
California Institute of Technology

© 2024 California Institute of Technology. Government sponsorship acknowledged.

Paper copies of this document may not be current and should not be relied on for official purposes. The current version is in the ECOSTRESS DocuShare Library (*) at [https://bravo-lib.jpl.nasa.gov/docushare/dsweb/View/Library-509](https://www.google.com/search?q=https://bravo-lib.jpl.nasa.gov/docushare/dsweb/View/Library-509)
(*) Access limited to user group

National Aeronautics and Space Administration
Jet Propulsion Laboratory
4800 Oak Grove Drive
Pasadena, California 91109-8099
California Institute of Technology

This research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration. Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise, does not constitute or imply its endorsement by the United States Government or the Jet Propulsion Laboratory, California Institute of Technology.
© 2024. California Institute of Technology. Government sponsorship acknowledged.

Note: The users' guide is designed to be a living document that describes the ECOSTRESS Land Surface Temperature and Emissivity (LST&E) product. The document describes the current state of the art, and is revised as progress is made in the development and assessment of the LST product. The primary purpose of the document is to present an overview of the ECOSTRESS L2 data product to the potential user. For more detailed information on the physical basis and algorithm details please see the Algorithm Theoretical Basis Document (ATBD).

## Change History Log

| Revision | Effective Date | Prepared by | Description of Changes |
| --- | --- | --- | --- |
| Draft | 6/4/2018 | Glynn Hulley | User Guide first draft based on MxD21/VNP21 products |
| Version 1 | 10/23/2018 | Glynn Hulley | Remove bit 4 from cloud mask product (table 7). Other small edits and clarifications in document. |
| Version 2 | 06/18/2019 | Glynn Hulley | Updates to account for MSU failure anomaly. |
| Version 3 | 04/11/2022 | Robert Freepartner<br>

<br>

<br>Glynn Hulley | Revision for ECOSTRESS Collection 2 (build 7*) |
| Version 4 | 10/5/2022 | Glynn Hulley | Final revision for ECOSTRESS Collection 2 to include a confidence level mask (build 7*) |
| Version 4.2 | 08/12/2024 | Glynn Hulley | Corrections to QC description of bit masks to account for cloud information not getting propagated to the QC bit mask in Collection 2. |
| Version 5 | 12/10/2025 | Glynn Hulley | Updates for Collection 3 |

## Contacts

Readers seeking additional information about this product may contact the following:
Glynn C. Hulley (PI)
MS 183-509
Jet Propulsion Laboratory
4800 Oak Grove Dr.
Pasadena, CA 91109
Email: glynn.hulley@jpl.nasa.gov
Office: (818) 354-2979

## Contents

* [Tables](https://www.google.com/search?q=%23tables)
* [1 Introduction](https://www.google.com/search?q=%231-introduction)
* [2 L2 products in Collection 3](https://www.google.com/search?q=%232-l2-products-in-collection-3)
* [2.1 HDF5-EOS5 Scene Gridded Products](https://www.google.com/search?q=%2321-hdf5-eos5-scene-gridded-products)
* [2.2 Cloud-Optimized GeoTIFF (COG) Tile products](https://www.google.com/search?q=%2322-cloud-optimized-geotiff-cog-tile-products)
* [2.3 Product Availability](https://www.google.com/search?q=%2323-product-availability)


* [3 Algorithm Description](https://www.google.com/search?q=%233-algorithm-description)
* [3.1 Background](https://www.google.com/search?q=%2331-background)
* [3.2 5-band versus 3-band TES algorithm](https://www.google.com/search?q=%2332-5-band-versus-3-band-tes-algorithm)


* [4 Data Product Structure and Content](https://www.google.com/search?q=%234-data-product-structure-and-content)
* [4.1 Scientific Data Sets (SDS)](https://www.google.com/search?q=%2341-scientific-data-sets-sds)
* [4.2 Attributes](https://www.google.com/search?q=%2342-attributes)
* [4.3 Quality Assurance (QA)](https://www.google.com/search?q=%2343-quality-assurance-qa)


* [5 L2G Cloud Product](https://www.google.com/search?q=%235-l2g-cloud-product)
* [5.1 Algorithm and product description](https://www.google.com/search?q=%2351-algorithm-and-product-description)
* [5.2 Cloud Scientific Data Sets (SDS)](https://www.google.com/search?q=%2352-cloud-scientific-data-sets-sds)
* [5.3 Attributes](https://www.google.com/search?q=%2353-attributes)


* [6 References](https://www.google.com/search?q=%236-references)

---

## Tables

* [Table 1: Summary of the Level-2 ECOSTRESS LST&E and cloud products in Collection 3](https://www.google.com/search?q=%23table-1)
* [Table 2: ECOSTRESS input products and ancillary data required to produce the L2 LST&E products](https://www.google.com/search?q=%23table-2)
* [Table 3: The Scientific Data Sets (SDS) in the ECOSTRESS L2G product](https://www.google.com/search?q=%23table-3)
* [Table 4: Standard product metadata included in the ECOSTRESS L2G product](https://www.google.com/search?q=%23table-4)
* [Table 5: Product specific metadata for the ECOSTRESS L2 product](https://www.google.com/search?q=%23table-5)
* [Table 6: Bit flags defined in the QC SDS in the L2G and L2T products](https://www.google.com/search?q=%23table-6)
* [Table 7: The SDSs in the ECOSTRESS L2G Cloud product](https://www.google.com/search?q=%23table-7)
* [Table 8: The metadata definition in the ECOSTRESS L2G Cloud product](https://www.google.com/search?q=%23table-8)

---

## 1 Introduction


This is the user guide for the ECOSTRESS Level-2 (L2) Land Surface Temperature and Emissivity (LST&E) and Cloud mask products. The L2 product uses a physics-based algorithm to dynamically retrieve both the LST&E simultaneously for the 3 and/or five ECOSTRESS thermal infrared bands at a spatial resolution of ~70×70 m. The algorithm is based on the ASTER Temperature Emissivity Separation (TES) algorithm, which uses full radiative transfer simulations for the atmospheric correction, and an emissivity model based on the variability in the surface radiance data to dynamically retrieve both LST and spectral emissivity. Simulations and validation results available in the ATBD have shown consistent accuracies over all surfaces with a total RMSE of 1.07 K. In addition, a Sea Surface Temperature (SST) algorithm is run on all pixels regardless of land, ocean, or cloud. The SST algorithm uses a view-angle dependent split-window algorithm with coefficients optimized for time of day, location and time of year and has an accuracy of <0.5 K.

Major changes to L2 products in C3:

* Swath-based products (no gridding) in HDF5 format are no longer produced.
* Only gridded and tiled products are produced.
* The L2 LSTE product includes a Sea Surface Temperature (SST) Product valid over all water surfaces (ocean and inland water).
* The L2 Cloud algorithm has been significantly improved.

In Collection 3 the ECOSTRESS L2 products will no longer include the swath-based product in standard geographic (lat, lon tagged) HDF5 format. This limits the number of artifacts that get translated from the L1 swath products to the L2 swath products that arise from striping and geolocation inconsistencies from areas with overlapping scans. The L2 product is now only produced in gridded HDF5 (L2G) and Cloud-Optimized GeoTIFF (COG) tiled (L2T) formats. The algorithms and data content of both L2G and L2T for the LSTE and Cloud products are briefly described in this guide, with the purpose of providing a user with sufficient information about the content and structure of the data files to enable the user to access and use the data, in addition to understanding the uncertainties involved with the product and how to interpret the cloud mask information.

Descriptive overviews of the file formats are provided first followed by descriptions of the algorithm and product contents including all metadata. Publications and documents related to the ECOSTRESS LST&E and cloud products are listed in the final section. A description of the major components of the ECOSTRESS algorithm implemented in version 1 of the LST&E Product Generation Executive (PGE) code are shown in Table 1 and described in depth in the ATBD available at [https://ecostress.jpl.nasa.gov/products](https://ecostress.jpl.nasa.gov/products). The primary purpose of this document is to supply a user with sufficient information about the content and structure of the data files so that the users will be able to access and use the data with confidence.

---

## 2 L2 products in Collection 3

### 2.1 HDF5-EOS5 Scene Gridded Products

The Hierarchical Data Format Earth Observing System 5 (HDF-EOS5) format is used to distribute ECOSTRESS granules at the scene level. These product files have a `.h5` file extension and are internally organized using the HDF-EOS5 data standard. The L2 gridded products using this format include the letter G in their level identifiers: L2G.

The HDF5 format is utilized here for long-term archiving, and is not recommended for end-user analysis. These HDF-EOS5 files are compatible with HDF-View, Panoply, and the `h5py` package in Python. Information on Hierarchical Data Format 5 (HDF5) may be found at [https://www.hdfgroup.org/HDF5/](https://www.hdfgroup.org/HDF5/). The HDF format was developed by NCSA, and has been widely used in the scientific domain.

Each of the raster layers in the L2G gridded products are projected to a globally snapped 0.0006° grid in WGS84 latitude and longitude to approximate 70 m resolution. The spatial metadata defining the raster grid is found using a plain-text format defined by HDF-EOS5 under:
`/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/StandardMetadata/`

Each product contains a group containing the product layers shown below following the HDF-EOS5 format: `HDFEOS/GRIDS/ECO_L2G_LSTE_70m/Data Fields/`

```text
HDF5 File (.h5)
├── /Data_Fields/
│   ├── LST (Land Surface Temperature)
│   ├── LST_err (LST uncertainty)
│   ├── EmisWB (Wideband emissivity)
│   ├── Emis1, Emis2, Emis3, Emis4, Emis5 (Band emissivities)
│   ├── QC (Quality control flags)
│   ├── SST (Sea Surface Temperature)
│   ├── PWV (Precipitable Water Vapor)
│   ├── cloud_mask
│   ├── water_mask
│   ├── height
│   └── view_zenith
├── /Geolocation/
│   ├── latitude
│   └── longitude
└── /Metadata/
    ├── StandardMetadata
    ├── ProductMetadata
    └── Attributes

```

Example reading code in Python:

```python
import h5py
# Open file
f = h5py.File('ECOv003_L2G_LSTE_15801_001_20210419T213114_01.h5', 'r')
# Access LST data
lst = f['Data_Fields']['LST'][:]
qc = f['Data_Fields']['QC'][:]
lat = f['Geolocation']['latitude'][:]
lon = f['Geolocation']['longitude'][:]
# Read metadata
metadata = dict(f['Metadata']['StandardMetadata'].attrs)
f.close()

```

### 2.2 Cloud-Optimized GeoTIFF (COG) Tile products

To provide an analysis-ready data format that has applicability across a range of software including QGIS, ArcGIS, Google Earth, Python, Julia etc. all ECOSTRESS Collection 3 products including the L2 products (L1CT, L2T, L3T, L4T) are also distributed in a tiled form using the COG format. The tiling system used for ECOSTRESS Collection 2 is derived from the modified Military Grid Reference System (MGRS) tiling scheme used by Sentinel 2. These tiles divide the Universal Transverse Mercator (UTM) zones into square tiles 109,800 m across. ECOSTRESS uses a 70 m cell size with 1568 rows by 1568 columns in each tile. This allows the end user to assume that each 70 m ECOSTRESS pixel will remain in the same location at each timestep observed in analysis.

The COG format also facilitates end-user analysis as a universally recognized and supported format, compatible with open-source software, including QGIS, ArcGIS, GDAL, the Raster package in R, rioxarray in Python, and Rasters.jl in Julia. Each `.tif` COG data layer in each L2T product additionally contains a rendered browse image in GeoJPEG format with a `.jpeg` extension. This image format is universally recognized and supported, and these files are compatible with Google Earth. Each tile includes `.json` file containing the Product Metadata and Standard Metadata in JSON format.

Each float32 data layer occupies 4 bytes of storage per pixel, which amounts to an uncompressed size of 13.4 MB for each tiled data layer. The uint8 quality flag layers occupy a single byte per pixel, which amounts to an uncompressed size of 3.35 MB per tiled data quality layer.

Each data layer for a specific tile are distributed as separate GeoTIFF files as shown in the example below. In addition there is one browse image per tile (GeoJPEG), and one metadata file per tile (JSON).

Example File Set for Single Tile:

```text
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_LST.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_LST_err.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_EmisWB.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis1.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis2.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis3.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis4.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis5.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis1_err.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis2_err.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis3_err.tif
ECOv003_T2T_LSTE_11SQA_20250909T164037_0713_01_Emis4_err.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_Emis5_err.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_cloud.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_water.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_height.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_view_zenith.tif
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_browse.jpeg
ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_metadata.json

```

Example reading code in Python:

```python
import rasterio
import rioxarray

# Option 1: Using rasterio
with rasterio.open('ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_LST.tif') as src:
    lst = src.read(1) # Read band 1
    transform = src.transform
    crs = src.crs
    metadata = src.tags()

# Option 2: Using rioxarray (with xarray)
lst_da = rioxarray.open_rasterio('ECOv003_L2T_LSTE_11SQA_20250909T164037_0713_01_LST.tif')
lst_data = lst_da.values[0] # First band

```

### 2.3 Product Availability

The ECOSTRESS L2 products are available by searching for the short names above in Table 1 at the NASA Land Processes Distribution Active Archive Center (LPDAAC), accessed via the Earthdata search engine, or in The Application for Extracting and Exploring Analysis Ready Samples (AppEEARS).


**Table 1: Summary of the Level-2 ECOSTRESS LST&E and cloud products in Collection 3.**

| Short name | Long name | Data Dimension | Map Projection |
| --- | --- | --- | --- |
| ECO_L2G_LSTE | ECOSTRESS Gridded Land Surface Temperature and Emissivity Instantaneous L2 Global 70 m | 5632 lines by 5400 pixels per line | Geographic, WGS84 |
| ECO_L2G_CLOUD | ECOSTRESS Gridded Cloud Mask Instantaneous L2 Global 70 m | 5632 lines by 5400 pixels per line | Geographic, WGS84 |
| ECO_L2T_LSTE | ECOSTRESS Tiled Land Surface Temperature and Emissivity Instantaneous L2 Global 70 m (100x100 km Sentinel tiles in UTM projection, cloud-optimized GeoTIFF) | 1568 lines by 1568 pixels per line | UTM, cloud-optimized GeoTIFF |

---

## 3 Algorithm Description

### 3.1 Background


For a full detailed description of each module within the L2 PGE please see the ECOSTRESS L2 ATBD. The algorithm uses a physical-based Temperature and Emissivity Separation (TES) algorithm to retrieve the Land Surface Temperature and Emissivity (LST&E) products (Gillespie et al. 1998; Hulley and Hook 2011). The atmospheric correction of the ECOSTRESS thermal infrared (TIR) bands 1-5 are performed using the RTTOV radiative transfer model (Matricardi 2008; Saunders et al. 1999) with input atmospheric profiles from the GEOS5 reanalysis product produced by the NASA Global Modeling and Assimilation Office (GMAO) (Rienecker et al. 2011). The GEOS5 data are provided on a ~1/3 degree longitude, 1/4 degree latitude spatial grid every 3 hours, with data provided in near real-time via ftp.

In collection 3 we include an additional Sea Surface Temperature (SST) layer that uses a view-angle dependent split-window algorithm to compute surface temperatures over oceans, coastal waters, and inland water bodies. The algorithm is based on a multi-channel SST (MCSST) algorithm that uses a set of optimized split-window coefficients, individually tuned for the atmospheric conditions of any location on a 0.5 degree grid (using GOES5 data). The emissivity is set to that of a library spectra of water from the ECOSTRESS spectral library.

$$T_s = a_o + a_1 T_i + a_2 (T_i - T_j) + a_3 (T_i - T_j)(1 - \sec(\theta))$$

The ECOSTRESS LST&E and SST products will be produced for all acquired ECOSTRESS scenes and for every pixel of data regardless of cloud. Note: we leave it up to the user to apply the land-water mask included with the product to separate LST and SST data over land/water as needed for their studies. The L2 product also includes a full set of error estimates for both the LST and emissivity bands generated from an uncertainty model (Hulley et al. 2012).

### 3.2 5-band versus 3-band TES algorithm

Due to the MSU failure anomalies, L2 products generated between May 15th 2019 and April 28th, 2023 use a 3-band version of the TES algorithm with bands 2, 4 and 5. This results in emissivity only being produced in those bands and the remaining bands will have fill values. The dropped bands will have no effect on the cloud mask algorithm that only uses bands 4 and 5. The retrieved LST with a 3-band approach will also result in slightly degraded accuracy (~0.1-0.2 K) when compared to the 5-band approach. More details on these changes and validation and uncertainty estimates are available by the science team or in Hulley et al. (2021). In build 7.0 processing starting in 2023, ECOSTRESS reverted back to downlinking all 5 bands of data due to improvements in optimization of how data is stored on the SDRAM in the latest DPU-IO 3.0 firmware updates.


**Table 2: ECOSTRESS input products and ancillary data required to produce the L2 LST&E product.**

| Ancillary Data Set | Long Name | Data Layers Used |
| --- | --- | --- |
| **ECOv003_L1CG_RAD** | ECOSTRESS Level-1 gridded, calibrated and geolocated radiances | • radiance_1...5<br>• Land/ocean mask<br>• Elevation<br>• Sensor view angles |
| **ASTER GEDv3** | ASTER Global Emissivity Dataset v3 | • Emis bands 10..14<br>• NDVI |
| **GEOS5-FP** | Atmospheric reanalysis data from the Global Modeling and Assimilation Office (GMAO) | • Pressure and geopotential height<br>• Temperature<br>• Specific Humidity<br>• Surface Pressure |

---

## 4 Data Product Structure and Content

### 4.1 Scientific Data Sets (SDS)

The ECOSTRESS L2G product contains 20 scientific data sets (SDSs) highlighted in Table 3. All SDS data are output at native ECOSTRESS 70m resolution pixels. The `*_err` SDSs are calculated using a LST&E uncertainty simulator and includes the maximum total uncertainty for a specific pixel based on view angle, total water vapor, and land cover type (Hulley et al. 2012b). Furthermore, a spatially and temporally interpolated Precipitable Water Vapor (PWV) estimate from GEOS5 data is included in the SDS as an indicator for the amount of water vapor present in the atmosphere—the primary driving factor for atmospheric correction uncertainty in retrieving LST&E.

Due to requests from users we further include a final estimate of cloudy pixels (`cloud_mask`), water mask (`water_mask`), elevation (`height`), and sensor view angle (`view_zenith`). The final estimate of cloud is detailed in section 5.1 below and can be used instead of downloading and interpreting the cloud confidence flags in the cloud product (see details in section 5.1). The water mask, height, and view_zenith is extracted from the L1 product and can be useful for investigators using the L2 data for research related to water and water quality that don’t need to extract this information from the L1B products. Details of each SDS including fill and scale factors are shown in Table 3. See section 2.2 above for the file set of the L2T product.


**Table 3. The Scientific Data Sets (SDS) in the ECOSTRESS L2G product.**

| SDS | Long Name | Data type | Units | Valid Range | Fill Value | Scale Factor | Offset |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LST | Land Surface Temperature | uint16 | K | 7500-65535 | 0 | 0.02 | 0.0 |
| QC | Quality control for LST and emissivity | uint16 | n/a | 0-65535 | 0 | n/a | n/a |
| Emis1 | Band 1 emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| Emis2 | Band 2 emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| Emis3 | Band 3 emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| Emis4 | Band 4 emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| Emis5 | Band 5 emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| LST_Err | Land Surface Temperature error | uint8 | K | 1-255 | 0 | 0.04 | 0.0 |
| Emis1_Err | Band 1 emissivity error | uint16 | n/a | 0-65535 | 0 | 0.0001 | 0.0 |
| Emis2_Err | Band 2 emissivity error | uint16 | n/a | 0-65535 | 0 | 0.0001 | 0.0 |
| Emis3_Err | Band 3 emissivity error | uint16 | n/a | 0-65535 | 0 | 0.0001 | 0.0 |
| Emis4_Err | Band 4 emissivity error | uint16 | n/a | 0-65535 | 0 | 0.0001 | 0.0 |
| Emis5_Err | Band 5 emissivity error | uint16 | n/a | 0-65535 | 0 | 0.0001 | 0.0 |
| EmisWB | Wideband emissivity | uint8 | n/a | 1-255 | 0 | 0.002 | 0.49 |
| SST | Sea Surface Temperature | uint16 | K | 7500-65535 | 0 | 0.02 | 0.0 |
| PWV | Precipitable Water Vapor | uint16 | cm | 0-65535 | 0 | 0.001 | 0.0 |
| cloud_mask | Cloud Mask | uint8 | 1=cloud | 0-1 | 255 | 1 | 0 |
| water_mask | Water Mask | uint8 | 1=water | 0-1 | 255 | 1 | 0 |
| view_zenith | Sensor view zenith angle | Float32 | degrees | 0-35 | 0 | n/a | n/a |
| height | Ground elevation | Float32 | m | n/a | 0 | n/a | n/a |

### 4.2 Attributes

The L2G product `.h5` file contains two sets of product metadata:

* `HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/ProductMetadata` (see Table 5)
* `HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/StandardMetadata` (see Table 4)

While the L2G contains a custom set of ProductMetadata attributes, the StandardMetadata attributes are consistent across products at each orbit/scene, as listed in Table 4. The L2T tiles include a `.json` file containing the Product Metadata and Standard Metadata in JSON format.


**Table 4. Standard product metadata included in the ECOSTRESS L2G product.**

| Name | Type | Size | Examples |
| --- | --- | --- | --- |
| **Group** |  |  | **StandardMetadata** |
| AncillaryInputPointer | String | variable | Group name of ancillary file list |
| AutomaticQualityFlag | String | variable | PASS/FAIL (of product data) |
| BuildId | String | variable |  |
| CollectionLabel | String | variable |  |
| DataFormatType | String | variable | NCSAHDF5 |
| DayNightFlag | String | variable |  |
| EastBoundingCoordinate | LongFloat | 8 |  |
| HDFVersionId | String | variable | 1.8.16 |
| ImageLines | Int32 | 4 | 5632 |
| ImageLineSpacing | Float32 | 4 | 68.754 |
| ImagePixels | Int32 | 4 | 5400 |
| ImagePixelSpacing | Float32 | 4 | 65.536 |
| InputPointer | String | variable |  |
| InstrumentShortName | String | variable | ECOSTRESS |
| LocalGranuleID | String | variable |  |
| LongName | String | variable | ECOSTRESS |
| NorthBoundingCoordinate | LongFloat | 8 |  |
| PGEName | String | variable | L2_LSTE (L2_CLOUD) |
| PGEVersion | String | variable |  |
| PlatformLongName | String | variable | ISS |
| PlatformShortName | String | variable | ISS |
| PlatformType | String | variable | Spacecraft |
| ProcessingLevelID | String | variable | 1 |
| ProcessingLevelDescription | String | variable | Level 2 Land Surface Temperatures and Emissivity (Level 2 Cloud mask) |
| ProducerAgency | String | variable | JPL |
| ProducerInstitution | String | variable | Caltech |
| ProductionDateTime | String | variable |  |
| ProductionLocation | String | variable |  |
| CampaignShortName | String | variable | Primary |
| RangeBeginningDate | String | variable |  |
| RangeBeginningTime | String | variable |  |
| RangeEndingDate | String | variable |  |
| RangeEndingTime | String | variable |  |
| SceneID | String | variable |  |
| ShortName | String | variable | L2_LSTE (L2_CLOUD) |
| SISName | String | variable |  |
| SISVersion | String | variable |  |
| SouthBoundingCoordinate | LongFloat | 8 |  |
| StartOrbitNumber | String | variable |  |
| StopOrbitNumber | String | variable |  |
| WestBoundingCoordinate | LongFloat | 8 |  |


**Table 5. Product specific metadata for the ECOSTRESS L2 product.**

| Name | Type | Size | Example |
| --- | --- | --- | --- |
| **Group** |  |  | **L2 LSTE Metadata** |
| QAPercentCloudCover | Int | 4 | 80 |
| CloudMeanTemperature | LongFloat | 8 | 231 |
| CloudMaxTemperature | LongFloat | 8 | 275 |
| CloudMinTemperature | LongFloat | 8 | 221 |
| CloudSDevTemperature | LongFloat | 8 | 0.45 |
| QAFractionGoodQuality | Int | 4 | 0.7 |
| LSTGoodAvg | LongFloat | 8 | 285.4 |
| Emis1GoodAvg | LongFloat | 8 | 0.95 |
| Emis2GoodAvg | LongFloat | 8 | 0.95 |
| Emis3GoodAvg | LongFloat | 8 | 0.95 |
| Emis4GoodAvg | LongFloat | 8 | 0.95 |
| Emis5GoodAvg | LongFloat | 8 | 0.95 |
| AncillaryGEOS5 | Str | 255 | GEOS.fp.asm.inst3_3d_asm_Np.20140702_0000.V01 |
| BandSpecification | Float32 | µm | Wavelengths used in the L2 retrieval for bands 1-6: 1.6, 8.2, 8.7, 9.0, 10.5, 12.0; 0=fill data |

### 4.3 Quality Assurance (QA)

Indicators of quality are described exclusively in the quality control (QC) SDS generated during production. In addition to data quality, the QC SDS provides information on algorithm metrics for each pixel (e.g. convergence statistics). The QC SDS unsigned 16-bit data are stored as bit flags in the SDS. This QC information can be extracted by reading the bits in the 16-bit unsigned integer. The purpose of the QC SDS is to give the user information on algorithm results for each pixel that can be viewed in a spatial context. The bit flags in the QC SDS are listed in Table 6 and consist of flags related to missing data, algorithm diagnostics, and error estimates.

A value for bits 1&0 = 00 in the QC bit flags indicates that the LST and emissivity was retrieved, but the pixel may or may not be cloudy. The user should either apply the final `cloud_mask` layer in the standard product or explore the confidence flags in the L2 Cloud mask to determine if a pixel is cloudy or clear.

Users should be aware that in some cases the TES algorithm may produce data with degraded quality (bits 1&0 = 01) due to one or more of the following conditions:

1. The retrieved emissivity in both longwave bands 4 (10.6 micron) and 5 (12 micron) is < 0.95 indicating possible cloud contamination.
2. The pixel falls on a missing scan line in bands 1 and 5, in which the radiance was filled using a spatial neural net technique (see Appendix A for more details). The user should check error estimates for this pixel to see if they fall within tolerable bounds.
3. The pixel had transmissivity less than 0.4 indicating either possible cloud contamination or high humidity, which would result in higher uncertainty in the TES retrieval. The user is encourage to check error estimates before using this pixel for science analysis.

A value for bits 1&0=11 indicates that the pixel was not produced due to poorly calibrated or missing radiance data, or the TES algorithm failed to converge (rare).


**Table 6. Bit flags defined in the QC SDS in the L2G and L2T products. (Note: Bit 0 is the least significant bit).**

| Bits | Long Name | Description |
| --- | --- | --- |
| 1&0 | Mandatory QA flags | 00 = Pixel produced by TES<br><br>01 = Pixel produced by TES but either one or more of the following conditions are met:<br><br>&nbsp;&nbsp;&nbsp;&nbsp;1. Emissivity in both bands 4 and 5 < 0.95, i.e. possible cloud contamination<br><br>&nbsp;&nbsp;&nbsp;&nbsp;2. Low transmissivity due to high water vapor loading (<0.4), check PWV values and error estimates<br><br>&nbsp;&nbsp;&nbsp;&nbsp;3. Pixel falls on missing scan line in bands 1&5, and filled using spatial neural net. Check error estimates.<br><br>10 = not set<br><br>11 = Pixel not produced due to missing/bad data, or TES divergence, user should check data quality flag bits. |
| 3 & 2 | Data quality flag | 00 = Good quality L1B data<br><br>01 = Missing stripe pixel in bands 1 and 5<br><br>10 = not set<br><br>11 = Missing/bad L1B data |
| 5 & 4 | Cloud/Ocean Flag | Not set. Please check ECOSTRESS GEO and CLOUD products for this information. |
| 7 & 6 | Iterations | 00 = Slow convergence<br><br>01 = Nominal<br><br>10 = Nominal<br><br>11 = Fast |
| 9 & 8 | Atmospheric Opacity | 00 = >=3 (Warm, humid air; or cold land)<br><br>01 = 0.2 - 0.3 (Nominal value)<br><br>10 = 0.1 - 0.2 (Nominal value)<br><br>11 = <0.1 (Dry, or high altitude pixel) |
| 11 & 10 | MMD | 00 = > 0.15 (Most silicate rocks)<br><br>01 = 0.1 - 0.15 (Rocks, sand, some soils)<br><br>10 = 0.03 - 0.1 (Mostly soils, mixed pixel)<br><br>11 = <0.03 (Vegetation, snow, water, ice) |
| 13 & 12 | Emissivity accuracy (Average of all bands) | 00 = >0.02 (Poor performance)<br><br>01 = 0.015 - 0.02 (Marginal performance)<br><br>10 = 0.01 - 0.015 (Good performance)<br><br>11 = <0.01 (Excellent performance) |
| 15 & 14 | LST accuracy | 00 = >2 K (Poor performance)<br><br>01 = 1.5 - 2 K (Marginal performance)<br><br>10 = 1 - 1.5 K (Good performance)<br><br>11 = <1 K (Excellent performance) |

---

## 5 L2G Cloud Product

The ECOSTRESS Level-2 Cloud product is output in a separate file to the LST&E product and consists of two binary layers (one cloud confidence flag and one final cloud mask) shown in Table 7. Full details of the algorithm and thresholds employed are available in the Cloud mask ATBD at [https://ecostress.jpl.nasa.gov/data](https://ecostress.jpl.nasa.gov/data).

### 5.1 Algorithm and product description


LST data can only be reliably retrieved under clear-sky conditions, making accurate cloud identification an important component of accurate LST and evapotranspiration estimation. Detecting clouds is a challenging endeavor and depends on not only the type of cloud being detected, but also the type of surface over which the cloud is detected. Clouds are brighter and colder than the land surface they obscure and these properties can be exploited using a combination of visible, shortwave infrared and thermal data. We provide a separate L2 cloud mask product that includes confidence level estimate of cloudy and clear-sky conditions, and a final estimate of cloudy pixels on an ECOSTRESS scene shown in Table 3. Note that the L2 algorithm will run on all pixels regardless of cloud, primarily due to large uncertainties in cloud detection when using only thermal bands.

In Collection 1 of the L2 data (build 6.*), the L2 cloud mask was estimated by thresholding band 4 brightness temperatures and band 4-5 brightness temperature differences, however, these thresholds were not optimized for time of day and variations in land surface emissivity resulting in possible cloud omission and commission errors for certain difficult case scenarios (e.g. low warm clouds at night, cold clouds over cold surfaces such as ice/snow). For cases such as this the user would need to further explore the outputs from different cloud tests from the cloud mask and/or adapt and modify the cloud mask thresholds for their particular use case. In addition, longwave retrieved emissivity bands (e.g. 4 and 5) are usually good indicators of cloud contamination. e.g. band 4 emissivity values less than 0.9 should be regarded as suspect and possibly cloud contaminated in the presence of nearby cloud.

In Collection 2 and 3 (build 7/8*) the cloud mask algorithm was improved significantly, and instead of a thresholding approach, we employed a more rigorous single-channel Bayesian cloud thresholding look-up-table (LUT) approach (Bulgin et al. 2014). The LUT thresholds were derived from a climatology of atmospheric conditions and surface emissivity using Monte-Carlo radiative transfer simulations to model distributions of clear sky radiances for a given sensor’s spectral response function in the longwave infrared window region (10-12 µm). Specifically, we use an emissivity database (ASTER GEDv3) and reanalysis data on 3 hour intervals (from GEOS5-FP) to derive clear-sky ECOSTRESS brightness temperature for each global grid point (~1.0° resolution) as a function of time of day, month of year, and location. Please refer to the L2 cloud ATBD for details on how the cloud thresholds are computed and how the final mask is determined. The L2 cloud mask product will contain binary layers instead of a bitmask in Collection 2/3 for ease of use (see Table 7).

In collection 3 the cloud mask was further refined and improved upon by increasing the resolution of the input atmospheric modeling grids from GEOS5 from 1 to 0.5° resolution. The surface temperature data in the simulations were improved by directly using the 2d GEOS5 modeling grid data product (‘skin temperature’), instead of estimating the surface temperature from the near surface temperature. This results in higher accuracy along coastal/land interfaces and improved discrimination of low, warm clouds from land.

Table 7 details the data sets included in the L2G cloud output. Users can interpret the data sets as follows:

1. `Cloud_confidence` contains the results of the brightness temperature LUT test with confidence levels set according to different threshold levels: 0 = confident clear, 1 = probably clear, 2 = probably cloudy, and 3 = confident cloudy.
2. `Cloud_final` contains a final cloud mask (1=cloud, 0=clear) based on the following criteria:
a. For elevations < 2km, cloud = probably cloud + confident cloudy pixels
b. For elevations > 2km, cloud = confident cloudy pixels

While users can interpret the confidence levels and estimation of a final cloud mask as they wish, the generalized rationale for the above logic is that mountainous areas with higher elevation will have higher uncertainties in the cloud mask due to the presence of snow/ice and effects of shading and varying temperature lapse rates. As a result only confident cloudy pixels are classified as cloud. For low lying regions using both probably and confident cloudy pixels provided a good balance between false positive and false negative cloud errors. If the user has zero tolerance for any nearby or possible cloud, then only confident clear pixels should be used, similarly if some cloud can be tolerated, then only confident cloudy pixels can be used.

### 5.2 Cloud Scientific Data Sets (SDS)


**Table 7. The SDSs in the ECOSTRESS L2G Cloud product.**

| SDS | Long Name | Data type | Units | Valid Range | Fill Value | Scale Factor | Offset |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **Cloud_confidence** | Brightness temperature LUT test | uint8 | 3 = confident cloudy<br>2 = probably cloudy<br>1 = probably clear<br>0 = confident clear | 0-1 | 255 | 1 | 0 |
| **Cloud_final** | Final cloud mask | uint8 | 1 = cloud<br>0 = clear | 0-1 | 255 | 1 | 0 |

### 5.3 Attributes


**Table 8. The metadata definition in the ECOSTRESS L2G Cloud product.**

| Name | Type | Size | Example |
| --- | --- | --- | --- |
| **Group** |  |  | **L2 CLOUD Metadata** |
| QAPercentCloudCover | Int | 4 | 80 |
| CloudMeanTemperature | LongFloat | 8 | 231 |
| CloudMaxTemperature | LongFloat | 8 | 275 |
| CloudMinTemperature | LongFloat | 8 | 221 |
| CloudSDevTemperature | LongFloat | 8 | 0.45 |

---

## 6 References

Gillespie, A., Rokugawa, S., Matsunaga, T., Cothern, J.S., Hook, S., Kahle, A.B., 1998. A temperature and emissivity separation algorithm for Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER) images. IEEE Transactions on Geoscience and Remote Sensing 36, 1113-1126.

Hulley, G.C., Gottsche, F., Rivera, G., Hook, S., Freepartner, R., Radocinksi, R., Martin, M., Cawse-Nicholson, K., Johnson, W.R., 2021. Validation and quality assessment of the ECOSTRESS level-2 land surface temperature and emissivity product. IEEE Transactions on Geoscience and Remote Sensing 59, 2772-2785.

Hulley, G.C., Hook, S.J., 2011. Generating consistent land surface temperature and emissivity products between ASTER and MODIS data for earth science research. IEEE Transactions on Geoscience and Remote Sensing 49, 1304-1315.

Hulley, G.C., Hook, S.J., Abbott, E., Malakar, N., Islam, T., Abrams, M., 2015. The ASTER Global Emissivity Dataset (ASTER GED): Mapping Earth's emissivity at 100 meter spatial scale. Geophysical Research Letters 42, 7966-7976.

Hulley, G.C., Hughes, C.G., Hook, S.J., 2012. Quantifying uncertainties in land surface temperature and emissivity retrievals from ASTER and MODIS thermal infrared data. Journal of Geophysical Research: Atmospheres 117, D20105.

Hulley, G.C., Hook, S.J., Schneider, P., 2011. Optimized split-window coefficients for deriving surface temperatures from inland water bodies. Remote Sensing of Environment 115, 3758-3769.

Malakar, N., Hulley, G.C., 2016. A water vapor scaling model for improved land surface temperature and emissivity separation of MODIS thermal infrared data. Remote Sensing of Environment 182, 252-264.

Matricardi, M., 2008. The generation of RTTOV regression coefficients for IASI and AIRS using a new profile training set and a new line-by-line database. ECMWF Research Department Technical Memorandum 564.

Rienecker, M.M., Suarez, M.J., Gelaro, R., Todling, R., Bacmeister, J., Liu, E., Bosilovich, M.G., Schubert, S.D., Takacs, L., Kim, G.K., Bloom, S., Chen, J.Y., Collins, D., Conaty, A., da Silva, A., Gu, W., Joiner, J., Koster, R.D., Lucchesi, R., Molod, A., Owens, T., Pawson, S., Pegion, P., Redder, C.R., Reichle, R., Robertson, F.R., Ruddick, A.G., Sienkiewicz, M., Woollen, J., 2011. MERRA: NASA's Modern-Era Retrospective Analysis for Research and Applications. Journal of Climate 24, 3624-3648.

Saunders, R., Matricardi, M., Brunel, P., 1999. An improved fast radiative transfer model for assimilation of satellite radiance observations. Quarterly Journal of the Royal Meteorological Society 125, 1407-1425.

Tonooka, H., 2005. Accurate atmospheric correction of ASTER thermal infrared imagery using the WVS method. IEEE Transactions on Geoscience and Remote Sensing 43, 2778-2792.