# GEE-SedimentBalance

Script and example to identify and quantitatively analyze deposition and erosion patterns at lake shores and in lake beds.
Appropriate for lakes and reservoirs with varying water surface heights (WSH).
Example script based on Landsat Collection 1 Level-1 surface reflectance imagery.

The methodology is published here (preprint):
#### [https://esurf.copernicus.org/preprints/esurf-2021-99/](https://esurf.copernicus.org/preprints/esurf-2021-99/)

Main steps:
1. COREGISTRATION: Co-Registration Landsat 5/7/8 with orthoimage
2. ENSEMBLE OF LANDSAT SCENES
3. SHORELINES OF LANDSAT SCENES
4. MEDIAN ELEVATION OF SHORELINES
5. SHORELINE ANOMALIES PER YEAR
6. PERFORM ADDITIONAL ITERATIONS
7. CALCULATE SEDIMENT BALANCE PER PIXEL WITHIN AOI
8. OVERALL LAKE SEDIMENT BALANCE

Run the example script in Earth Engine:
#### [https://code.earthengine.google.com/9329201c7feedd54608db44d51032368](https://code.earthengine.google.com/9329201c7feedd54608db44d51032368)
