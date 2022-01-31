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

//date for coregistration:
var selected_date=ee.Date.fromYMD(2019,5,11);

//list of band names
var bandIn_s2 = ['B2','B3','B4','B8','B11','B12'];
var bandOut = ['blue','green','red','nir','swir1','swir2'];

//Load the Sentinel-2 iamge
var imgS2 = ee.Image(ee.ImageCollection('COPERNICUS/S2_SR').filterBounds(aoi)
       .filterDate(selected_date,selected_date.advance(1,'day')).first()).select(bandIn_s2,bandOut);
// Choose to register using only the 'green' band.
var S2GreenBand = imgS2.select('green');
// Get the displacement between the Ortho and the Sentinel-2 image
var displacement_S2_ortho = S2GreenBand.displacement({
  referenceImage: Ortho.select('green'),
  maxOffset: 50.0,//The maximum offset allowed when attempting to align the input images, in meters
});

// Determine the displacement by matching the 'Green' bands.
var displacement_L7 = gee_codes.coreg_getdisplacement (ee.Date.fromYMD(2020,4,25),"LANDSAT_7",S2GreenBand.displace(displacement_S2_ortho),aoi);
var displacement_L5 = displacement_L7;
var displacement_L8 = gee_codes.coreg_getdisplacement (ee.Date.fromYMD(2020,11,11),"LANDSAT_8",S2GreenBand.displace(displacement_S2_ortho),aoi);
print('displacement_L7',displacement_L7);

// ****************************************************************************************************************** //
//ENSEMBLE OF LANDSAT SCENES
// ****************************************************************************************************************** //

//define the study period
var start_date = ee.Date.fromYMD(1999,10,1);
var end_date = ee.Date.fromYMD(2021, 10, 1);
var merged_landsat = ee.ImageCollection(gee_codes.get_landsat_images(aoi,displacement_L5,displacement_L7,displacement_L8,start_date,end_date,
  {nodata_thresh: 30,
  cc_thresh: 30,
  nodata_thresh_l7: 50,
  }));
print('Number of scenes in collection:',merged_landsat.size());

// ****************************************************************************************************************** //
//SHORELINES OF ALL LANDSAT SCENES
// ****************************************************************************************************************** //

var edge_imgs_landsat = ee.ImageCollection(merged_landsat
    //returns the DEM clipped to the shorelines for each scene
      .map(function(img){return gee_codes.edgeOtsu_optical(img,dem_sylvatrop,
        {cannyThreshold: 0.7,
          cannySigma: 0.7,
          minimumValue: -0.3,
          maximumValue: 0.3,
          scale: 30,
          th_veg: 0.5
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
print('edge_sed_landsat',edge_sed_landsat.first());

// ****************************************************************************************************************** //
//SHORELINE ANOMALIES PER YEAR
// ****************************************************************************************************************** //
//returns a stacked image (one image with one band per year)
var annual_deltah=ee.Image(ee.ImageCollection(gee_codes.getannual_deltah(edge_sed_landsat,start_date,end_date)));
print('annual_deltah',annual_deltah);

//export the stack as an asset
  var thisname='ID_1009';
  Export.image.toAsset({
    image:annual_deltah.set('iterations',1),
    description: 'SedZone_L8_00to20_' + thisname + '_1st_iteration',
    assetId: 'users/magosilvain/Wegnia_paper/SedZones/stack_test/SedZone_L8_00to20_' + thisname + '_1st_iteration',
    scale: 10,
    maxPixels: 1e13,
    region:aoi
  });
// ****************************************************************************************************************** //
//CALCULATE SEDIMENT BALANCE PER PIXEL WITHIN AOI
// ****************************************************************************************************************** //
var sed_thislake=ee.Image(gee_codes.stacktoimage(annual_deltah));
print('sed_thislake',sed_thislake);

// ****************************************************************************************************************** //
//NOW PERFORM SEVERAL ITERATIONS
// ****************************************************************************************************************** //

//NUMBER OF ADDITIONAL ITERATIONS
var it= 8;

var sed_thislake=ee.Image(ee.List.sequence(1,it).iterate(function(x,previous){
  previous=ee.Image(previous);
  //MEDIAN ELEVATION OF SHORELINES
  var edge_img4DEM = ee.ImageCollection(gee_codes.get_corrected_elevation_DEM(edge_imgs_landsat, previous, {
    filter_outliers: true, //optional: filter out the noisiest shorelines (dh_std >= two standard deviations of the entire ensemble)
    }));
  //SHORELINE ANOMALIES OF EACH SCENE
  var edge_sed_landsat=ee.ImageCollection(edge_img4DEM.map(function(img){return gee_codes.get_sedimentation_zone(img,start_date)}));
  //SHORELINE ANOMALIES PER YEAR
  var annual_deltah=ee.Image(ee.ImageCollection(gee_codes.getannual_deltah(edge_sed_landsat,start_date,end_date)));
  //SEDIMENT BALANCE
  var sed_thislake=ee.Image(gee_codes.stacktoimage(annual_deltah));
  return sed_thislake;
},sed_thislake));

//one last iteration to obtain again the shoreline anomalies per year
//MEDIAN ELEVATION OF SHORELINES
var edge_img4DEM = ee.ImageCollection(gee_codes.get_corrected_elevation_DEM(edge_imgs_landsat, sed_thislake, {
  filter_outliers: true, //optional: filter out the noisiest shorelines (dh_std >= two standard deviations of the entire ensemble)
  }));
//SHORELINE ANOMALIES OF EACH SCENE
var edge_sed_landsat=ee.ImageCollection(edge_img4DEM.map(function(img){return gee_codes.get_sedimentation_zone(img,start_date)}));
//SHORELINE ANOMALIES PER YEAR
var annual_deltah=ee.Image(ee.ImageCollection(gee_codes.getannual_deltah(edge_sed_landsat,start_date,end_date)));
  
  Export.image.toAsset({
    image:annual_deltah.set('iterations',10),
    description: 'SedZone_L8_00to20_' + thisname + '_10th_iteration',
    assetId: 'users/magosilvain/Wegnia_paper/SedZones/stack_test/SedZone_L8_00to20_' + thisname + '_10th_iteration',
    scale: 10,
    maxPixels: 1e13,
    region:aoi
  });
  
  
