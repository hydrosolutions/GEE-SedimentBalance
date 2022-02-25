// ****************************************************************************************************************** //
// IMPORT LAYERS AND CODES
// ****************************************************************************************************************** //
// Import external dependencies
var gee_codes = require('users/hydrosolutions/Sahel_Apps:SedimentBalance');

//Import DEM
var dem_sylvatrop = ee.Image("users/tabeadonauer/DEM_Abijata/LAC_WEGNIA_MNT");
//Import Ortho image
var Ortho = ee.Image("users/magosilvain/Ortho_Wegnia_5m_clip");
//rename Bands
Ortho = Ortho.select(['b1', 'b2', 'b3']).rename(['red', 'green', 'blue']);
//Import a shapefile of the area of interest
var aoi = ee.FeatureCollection('users/magosilvain/Lakes_Sahel_b500').filter(ee.Filter.eq('ID', 1009));
Map.centerObject(aoi);
//keep only the geometry, add a buffer
aoi=aoi.first().geometry().buffer(500);//500m buffer around polygon
Map.addLayer(ee.Image().paint(aoi, 0, 2),null, "area_of_interest");
Map.setOptions("HYBRID");
//alternatively: use the footprint of the available DEM
//var aoi=ee.Feature(ee.Geometry(Ortho.get('system:footprint'))).geometry();

// ****************************************************************************************************************** //
//COREGISTRATION
// ****************************************************************************************************************** //
/*Note: This example is for Landsat Level-1 products. However, several data processing, geometric, and radiometric improvements, 
along with a new data distribution process, define Collection 2 Level-1 data, that are now also available on GEE.
*/
//date for coregistration:
var selected_date=ee.Date.fromYMD(2019,5,11);

//list of band names
var bandIn_sentinel = ['B2','B3','B4','B8','B11','B12'];
var bandIn_landsat = ['B1','B2','B3','B4','B5','B7'];
var bandOut = ['blue','green','red','nir','swir1','swir2'];

//Load the Sentinel-2 image
var imgS2 = ee.Image(ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(aoi)
       .filterDate(selected_date,selected_date.advance(1,'day')).first()).select(bandIn_sentinel,bandOut);
// Choose to register using only the 'green' band.
var S2GreenBand = imgS2.select('green');
// Get the displacement between the Ortho and the Sentinel-2 image
var displacement_S2_ortho = S2GreenBand.displacement({
  referenceImage: Ortho.select('green'),
  maxOffset: 50.0,//The maximum offset allowed when attempting to align the input images, in meters
});

// Determine the displacement by matching the 'Green' bands.
var selected_date_S2=ee.Date.fromYMD(2020,4,25);
var imgS2 = ee.Image(ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(aoi)
       .filterDate(selected_date_S2,selected_date_S2.advance(1,'day')).first()).select(bandIn_sentinel,bandOut);
var displacement_L7 = gee_codes.coreg_getdisplacement(ee.Date.fromYMD(2020,4,25),"LANDSAT_7",imgS2.select('green').displace(displacement_S2_ortho),aoi);

var selected_date_L7=ee.Date.fromYMD(2007,6,9);
var selected_date_L5=ee.Date.fromYMD(2007,6,1);
var imgL7 = ee.Image(ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterBounds(aoi)
       .filterDate(selected_date_L7,selected_date_L7.advance(1,'day')).first()).select(bandIn_landsat,bandOut);
var displacement_L5 = gee_codes.coreg_getdisplacement(selected_date_L5,"LANDSAT_5",imgL7.select('green').displace(displacement_L7),aoi);

var selected_date_S2=ee.Date.fromYMD(2020,11,11);
var imgS2 = ee.Image(ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(aoi)
       .filterDate(selected_date_S2,selected_date_S2.advance(1,'day')).first()).select(bandIn_sentinel,bandOut);
var displacement_L8 = gee_codes.coreg_getdisplacement(ee.Date.fromYMD(2020,11,11),"LANDSAT_8",imgS2.select('green').displace(displacement_S2_ortho),aoi);

Map.addLayer(Ortho,  {bands: (['red', 'green', 'blue']), min: 30, max: 220},'Ortho RGB');   
var imgL5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterBounds(aoi)
    .filterDate(selected_date_L5,selected_date_L5.advance(1,'day')).select(bandIn_landsat,bandOut).mosaic().clip(aoi);
    
Map.addLayer(imgL5, {bands: (['red', 'green', 'blue']), min: 0, max: 3000}, 'RGB Landsat 5 Raw',false);
Map.addLayer(imgL5.displace(displacement_L5), {bands: (['red', 'green', 'blue']), min: 0, max: 3000}, 'RGB Landsat 5 coregistrated', false);

// ****************************************************************************************************************** //
//ENSEMBLE OF LANDSAT SCENES
// ****************************************************************************************************************** //

//define the study period
var start_date = ee.Date.fromYMD(1999,10,1);
var end_date = ee.Date.fromYMD(2021, 10, 1);
// number of years
var years = end_date.difference(start_date,'year').round();
print('number of years considered for analysis',years);

var merged_landsat = ee.ImageCollection(gee_codes.get_landsat_images(aoi,displacement_L5,displacement_L7,displacement_L8,start_date,end_date,
  {nodata_thresh: 30,
  cc_thresh: 30,
  nodata_thresh_l7: 50,
  }));
print('Number of scenes in collection:',merged_landsat.size());

// ****************************************************************************************************************** //
//SHORELINES OF ALL LANDSAT SCENES
// ****************************************************************************************************************** //
/*Note: Water detection currently can't handle floating vegetation/inundated vegetation. 
With the parameter 'th_veg' we can ignore pixels adjaccent to vegetation.*/ 

var edge_imgs_landsat = ee.ImageCollection(merged_landsat
    //returns the DEM clipped to the shorelines for each scene
      .map(function(img){return gee_codes.edgeOtsu_optical(img,dem_sylvatrop,
        {cannyThreshold: 0.7,
          cannySigma: 0.7,
          minimumValue: -0.3,
          maximumValue: 0.3,
          scale: 30,
          th_veg: 0.5 //threshold mean NDVI during the main wet season months (July to September)
        })}));

// ****************************************************************************************************************** //
//MEDIAN ELEVATION OF SHORELINES
// ****************************************************************************************************************** //

//map over the images to obtain the median elevation and also the standard deviation of the elevation ('dh_std')
var edge_img4DEM = ee.ImageCollection(gee_codes.get_corrected_elevation_DEM(edge_imgs_landsat, null, {
  filter_outliers: true, //optional: filter out the noisiest shorelines (dh_std >= two standard deviations of the entire ensemble)
  }));
print('Number of scenes with data available:',edge_img4DEM.size());

// ****************************************************************************************************************** //
//SHORELINE ANOMALIES OF EACH SCENE
// ****************************************************************************************************************** //

////Extract Shoreline Anomalies of individual contours and get the time information required for linear regression
var edge_sed_landsat=ee.ImageCollection(edge_img4DEM.map(function(img){return gee_codes.get_sedimentation_zone(img,start_date)}));

// ****************************************************************************************************************** //
//SHORELINE ANOMALIES PER YEAR
// ****************************************************************************************************************** //
//returns a stacked image (one image with one band per year)
var annual_deltah=ee.Image(gee_codes.getannual_deltah(edge_sed_landsat,start_date,end_date));

//define asset path and name
var imgname= 'SedZone_00to21';
var path= 'users/magosilvain/Wegnia_paper/example_stack/';

//export the stack as an asset
  Export.image.toAsset({
    image:annual_deltah.set('iterations',1),
    description: imgname + '_1st_iteration',
    assetId: path + imgname + '_1st_iteration',
    scale: 10,
    maxPixels: 1e13,
    region:aoi
  });

// ****************************************************************************************************************** //
//NOW PERFORM SEVERAL ADDITIONAL ITERATIONS
// ****************************************************************************************************************** //

//NUMBER OF ADDITIONAL ITERATIONS
var it= 4;

annual_deltah=ee.Image(ee.List.sequence(1,it).iterate(function(x,previous){
  previous=ee.Image(previous);
  //SEDIMENT BALANCE
  var sed_thislake=ee.Image(gee_codes.stacktoimage(ee.Image(previous)));
  //MEDIAN ELEVATION OF SHORELINES
  var edge_img4DEM = ee.ImageCollection(gee_codes.get_corrected_elevation_DEM(edge_imgs_landsat, sed_thislake, {
    filter_outliers: true, //optional: filter out the noisiest shorelines (dh_std >= two standard deviations of the entire ensemble)
    }));
  //SHORELINE ANOMALIES OF EACH SCENE
  var edge_sed_landsat=ee.ImageCollection(edge_img4DEM.map(function(img){return gee_codes.get_sedimentation_zone(img,start_date)}));
  //SHORELINE ANOMALIES PER YEAR
  var annual_deltah=ee.Image(gee_codes.getannual_deltah(edge_sed_landsat,start_date,end_date));

  return annual_deltah;
},annual_deltah));

//Unfortunately, for this example, the execution of several iterations at once will fail, due to out of memory. (Error code: 8)
//It is necessary to export after every 1-2 iterations
Export.image.toAsset({
  image:annual_deltah.set('iterations',5),
  description: imgname + '_5th_iteration',
  assetId: path +  imgname + '_5th_iteration',
  scale: 10,
  maxPixels: 1e13,
  region:aoi
});
  
// ****************************************************************************************************************** //
//CALCULATE SEDIMENT BALANCE PER PIXEL WITHIN AOI
// ****************************************************************************************************************** //

//if already exported, load the annual shoreline anomalies from the asset
var annual_deltah=ee.Image(path +  imgname +  '_10th_iteration');
//Sen's slope of shoreline anomalies per pixel --> sediment balance
var sed_thislake=ee.Image(gee_codes.stacktoimage(annual_deltah));

// ****************************************************************************************************************** //
//OVERALL LAKE SEDIMENT BALANCE
// ****************************************************************************************************************** //  
//convert the stack of shoreline anomalies back to an image collection
var sed_zone_collection=ee.ImageCollection(annual_deltah.bandNames().map(function(b){
  var nr=ee.Number.parse(ee.String(b).slice(1));
  return annual_deltah.select(ee.String(b)).rename('b1').set('t',nr);
}));
print('sed_zone_collection',sed_zone_collection);
//Count the number of available annual SEAs per pixel (2000 to 2021)
var count=sed_thislake.select('count');

var slope_final=
  sed_thislake.select('slope').updateMask(count.gte(6)).focal_mean({radius: 2, kernelType :'circle',units: 'pixels',iterations: 1})
  .blend(
    //focal_mean with radius 1 to smooth the image
    sed_thislake.select('slope').updateMask(count.gte(6)).focal_mean({radius: 1, kernelType :'circle',units: 'pixels',iterations: 1}).updateMask(count.gte(6))
  ).reproject({crs: count.projection().crs(), scale:10});
  
//Because surface water may not be correctly identified under canopies, mask all boundary pixels that are adjacent to lush vegetation
var th_veg=0.5;//threshold mean NDVI during the main wet season months (July to September)
var ndvi_wetseason=ee.Image(1).where(merged_landsat.filter(ee.Filter.calendarRange(7,9, 'month')).mean().select('NDVI').gte(th_veg),0);
slope_final=slope_final
  .updateMask(ndvi_wetseason.focal_min({radius: 10,kernelType:'square', units: 'meters'}))
  .updateMask(count.gte(1))//mask out all pixels that have never represented the shoreline
  .reproject({crs: count.projection().crs(), scale:10});

Map.addLayer(slope_final,{bands: (['slope']),palette:['red', 'white','blue'],min: -0.05, max: 0.05}, "erosion/deposition rates");

var dem_filled = ee.Image(-99).rename('b1').where(dem_sylvatrop.gt(0),dem_sylvatrop).clip(aoi);
var max_elevation_waterline=330.98;
var outflow_elevation=329.765;

print('average sediment balance in mm (negative: deposition, positive: erosion)',
  ee.Number(slope_final
  .reduceRegion(ee.Reducer.mean(),aoi).get('slope')).multiply(years).multiply(1000).round());

print('average sediment balance in mm, below outflow (negative: deposition, positive: erosion)',
  ee.Number(slope_final.updateMask(dem_sylvatrop.lt(outflow_elevation))
  .reduceRegion(ee.Reducer.mean(),aoi).get('slope')).multiply(years).multiply(1000).round());
  
print('average sediment balance in mm, above outflow (negative: deposition, positive: erosion)',
  ee.Number(slope_final.updateMask(dem_sylvatrop.gte(outflow_elevation))
  .reduceRegion(ee.Reducer.mean(),aoi).get('slope')).multiply(years).multiply(1000).round());
      
print('average sediment balance in mm at pixels with net deposition of sediments',
  ee.Number(slope_final.updateMask(slope_final.lte(-0.01))
.reduceRegion(ee.Reducer.mean(),aoi).get('slope')).multiply(years).multiply(1000).round());

print('average sediment balance in mm at pixels with net erosion',
  ee.Number(slope_final.updateMask(slope_final.gte(0.01))
.reduceRegion(ee.Reducer.mean(),aoi).get('slope')).multiply(years).multiply(1000).round());

print('Fraction of area with net erosion (%)',
  ee.Number(slope_final.updateMask(slope_final.gte(0.01))
  .reduceRegion(ee.Reducer.count(),aoi).get('slope'))
  .divide(ee.Number(slope_final.reduceRegion(ee.Reducer.count(),aoi).get('slope')))
  .multiply(100));
print('Fraction of area with net deposition (%)',
  ee.Number(slope_final.updateMask(slope_final.lte(-0.01))
  .reduceRegion(ee.Reducer.count(),aoi).get('slope'))
  .divide(ee.Number(slope_final.reduceRegion(ee.Reducer.count(),aoi).get('slope')))
  .multiply(100));
print('Fraction of area with stable ground (%)',
  ee.Number(slope_final.updateMask(slope_final.abs().lt(0.01))
  .reduceRegion(ee.Reducer.count(),aoi).get('slope')) 
  .divide(ee.Number(slope_final.reduceRegion(ee.Reducer.count(),aoi).get('slope')))
  .multiply(100));
