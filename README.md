# GEE-SedimentBalance

1. Co-Registration Landsat 5/7/8 and Sentinel-2

2. Water Edges Landsat 5/7/8 and Sentinel-2 ("edgeOtsu_func_optical")
3. Extract Elevation Information based on DEM ("get_elevation_DEM")
4. Extract Shoreline Anomalies of individual contours ("get_sedimentation_zone")
5. Calculate annual shoreline anomalies ("getannual_deltah") --> make image stack ("newCollectionToImage")
6. Pixel sediment balances based on annual anomalies and SensSlope ("stacktoimage")
7. Identify deposition and erosion zones based on a sediment balance threshold ("get_sedzone")

8. Iteratively Repeat steps 2-7
9. Calculate Lake sediment balance and extract final WSH Time Series
