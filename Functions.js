var morefunctions = require('users/magosilvain/Tabea_METRIC:Otsu_thresholding');

/**
 * Otsu-Edge thresholding based on the algorithm used by Donchyts et al.2016
 * 
 * @param {Image} image: An image with three bands: 'MNDWI','NDVI' and 'cloud'
 * @param {Dictionary} options: The options consist of the following parameters:
 * - @param {Number} cannyThreshold: 
 * - @param {Number} cannySigma: 
 * - @param {Number} minimumValue: 
 * - @param {Number} maximumValue: Cap the highest MNDWI values. Such high values represent rarely the edge of water.
 *      A maximum value reduces the risk of detecting edges between free water and vegetation overgrown water 
 * - @param {Number} scale: 
 * 
 *  
*/
exports.edgeOtsu_optical = function(image,DEM,options){
  DEM=ee.Image(DEM);
  var cannyThreshold = options.cannyThreshold || 0.7;
  var cannySigma = options.cannySigma || 0.7;
  var minValue = options.minimumValue || -0.3;
  var maxValue = options.maximumValue || 0.3;
  var scale = options.scale || 30;
  var th_veg = options.th_veg || 0.5;

  var bounds=ee.Image(image).get('system:footprint');
  var mndwi = ee.Image(image).select('MNDWI');
  mndwi=mndwi.where(mndwi.gt(maxValue),maxValue);
  var imgcloud = ee.Image(image).select('cloud');

  var th = morefunctions.computeThresholdUsingOtsu(mndwi.updateMask(imgcloud.lt(100)), scale, bounds, cannyThreshold, cannySigma, minValue);
  
  function getEdge10m(mask) {
  //get a buffer around edge +- 10m
    return mask.focal_max({radius: 10, units: 'meters'}).subtract(mask.focal_min({radius: 10, units: 'meters'}));
  }
  var dem_edge=DEM.updateMask(getEdge10m(mndwi.gt(th))).updateMask(imgcloud.lt(100));

  var veg_edge=ee.Image(1).updateMask(
    ee.Image(image).select('NDVI').gt(th_veg).focal_max({radius: 10, units: 'meters'}));
  veg_edge=ee.Image(0).where(veg_edge.eq(1),1); 
  
  return dem_edge.rename('dem_edge').clip(bounds)
    .addBands(veg_edge.rename('veg_edge'))
    .set('otsu_th',th)
    .copyProperties(image, ['system:time_start','SENSING_TIME', 'SATELLITE','cloud_cover','nodata_cover','ID']);  
};

//Extract median elevation of contour lines
exports.get_corrected_elevation_DEM=function(ic,sed_thislake,options){
  ic=ee.ImageCollection(ic);
  var filter_outliers= options.filter_outliers || true;
  var edge_img4DEM=ee.ImageCollection(ic.map(function(image){
    var bounds=ee.Image(image).get('system:footprint');
    var dem_edge=ee.Image(image).select('dem_edge');
    var veg_edge=ee.Image(image).select('veg_edge');
    
    if (sed_thislake){
      var slope_filled=ee.Image(sed_thislake).select('sedzone_filled');
      dem_edge=dem_edge.updateMask(slope_filled.abs().lt(0.01));
    }
    var dem_elevation = 
    //Ignore pixels adjaccent to vegetation. Water detection currently can't handle floating vegetation/inundated vegetation.
      ee.Number(dem_edge.updateMask(veg_edge.neq(1)).reduceRegion({
        reducer: 'median',
        geometry: bounds,
        scale: 10,
        tileScale: 2,
      }).get('dem_edge'));  
    dem_elevation =ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(dem_elevation, null), -100, dem_elevation));
  
    //the standard deviation is an indicator of noise. Can be used to filter out noisy images (e.g., edges due to clouds)
    var dh_std = 
      ee.Number(dem_edge.subtract(ee.Image(dem_elevation)).updateMask(veg_edge.neq(1)).updateMask(dem_edge.gt(-100))        
      .reduceRegion({
        reducer: ee.Reducer.stdDev(),
        geometry: bounds,
        scale: 10,
        tileScale: 2,
      }).get('dem_edge'));    
      
    return ee.Image(image).set('dem_elev',dem_elevation).set('dh_std',dh_std);
  }));
  //returns -100 if the median elevation cannot be calculated (e.g., shorelines are outside the provided DEM)
  edge_img4DEM= edge_img4DEM.filter(ee.Filter.gt('dem_elev', -100));

  if (filter_outliers){
    //optional: filter out the noisiest shorelines (dh_std >= two standard deviations of the entire ensemble)
    var dh_std_list=edge_img4DEM.aggregate_array('dh_std');
    var dh_std_threshold=ee.Number(dh_std_list.reduce(ee.Reducer.mean())).add(ee.Number(2).multiply(ee.Number(dh_std_list.reduce(ee.Reducer.stdDev()))));
    edge_img4DEM=edge_img4DEM.filter(ee.Filter.lt('dh_std', dh_std_threshold));
  }
  return edge_img4DEM;
};

//Extract Shoreline Anomalies of individual contours 
exports.get_sedimentation_zone=function(image,start_date){
  start_date=ee.Date(start_date);
  var dem_edge=ee.Image(image).select('dem_edge');
  var veg_edge=ee.Image(image).select('veg_edge');

  var dem_elevation=ee.Number(dem_edge.get('dem_elev'));
  var deltat=ee.Number(ee.Date(ee.Image(image).get('system:time_start')).difference(start_date,'year')).floor();
  dem_edge=dem_edge.updateMask(veg_edge.neq(1)).subtract(ee.Image(dem_elevation))
  .reduceNeighborhood(ee.Reducer.mean(),ee.Kernel.square(30, 'meters')).updateMask(dem_edge.gt(-100));
  return dem_edge.rename('sed_diff')
    .copyProperties(ee.Image(image)).set('t',deltat).set('system:time_start',ee.Image(image).get('system:time_start'));
};

//Calculate annual shoreline anomalies 
exports.getannual_deltah=function(daily_deltah,start_date,end_date){    
  start_date=ee.Date(start_date);
  end_date=ee.Date(end_date); 
  daily_deltah=ee.ImageCollection(daily_deltah);
  var start_year=start_date.get('Year');
  var start_month=start_date.get('Month');
  var start_day=start_date.get('Day');
  var end_year=end_date.get('Year');
  var annual_deltah=ee.ImageCollection(ee.List.sequence(ee.Number(start_year).add(1), end_year).map(function(year) {
    var col=daily_deltah.filterDate(ee.Date.fromYMD(ee.Number(year).subtract(1),start_month,start_day), ee.Date.fromYMD(ee.Number(year),start_month,start_day));
    var col_mean=col.reduce(ee.Reducer.mean()).rename('sed_diff')
      .set('system:time_start',ee.Date.fromYMD(ee.Number(year), 1, 1));
    return col_mean.set('count',col.size()).set('delta_t',ee.Number(year).subtract(ee.Number(start_year).add(1)));
  })).filter(ee.Filter.gt('count', 0));

  //stacking is useful for exporting (keep only one image with one band per year)
  var stack = ee.Image(annual_deltah.select('sed_diff').iterate(function(img, prev) {
    return ee.Image(prev).addBands(img.rename(ee.String('t').cat(ee.String(ee.Number(ee.Image(img).get('delta_t')).int()))));
  }, ee.Image(1))).slice(1);
  return stack;
};

//Pixel sediment balances based on annual anomalies and SensSlope.
//The output of this function can be used as an input to function 'get_corrected_elevation_DEM' above (sed_thislake)
exports.stacktoimage=function(img){
  var bands=ee.Image(img).bandNames();
  var annual_deltah=
    bands.map(function(str){
    var x=ee.Number.parse(ee.String(str).slice(1));
    return ee.Image(img).select(ee.String('t').cat(x.int())).rename('sed_diff')
    .addBands(ee.Image(x).rename('t').float())
    .set('system:time_start',ee.Date.fromYMD(x.add(2000), 1, 1));
  });
  var correlation = ee.ImageCollection(annual_deltah).select('t', 'sed_diff').reduce(ee.Reducer.pearsonsCorrelation());
  var reducerSensSlope=ee.Reducer.sensSlope();
  var im_reduced_sens=ee.ImageCollection(annual_deltah).select('t', 'sed_diff').reduce(reducerSensSlope)
    .addBands(ee.ImageCollection(annual_deltah).select('sed_diff').count().rename('count'))
    .addBands(correlation.select('p-value'));
  var slope=im_reduced_sens.select('slope');
  var count=im_reduced_sens.select('count');
  
  var sed_zone_collection=ee.ImageCollection(bands.map(function(b){
    var nr=ee.Number.parse(ee.String(b).slice(1));
    return ee.Image(img).select(ee.String(b)).rename('b1').set('t',nr);
  }));
  //Count the number of available annual SEAs per pixel (2000-2021)
  var count2=sed_zone_collection.filter(ee.Filter.lte('t', 10)).count().rename('count');
  slope=slope.updateMask(count.gte(6)).updateMask(count2.gt(1));
  
  //Count the number of available annual SEAs per pixel (recent years; 2011-present)
  var count3=sed_zone_collection.filter(ee.Filter.gt('t', 10)).count().rename('count');
  count3=count3.where(count2.gte(2),0);
  //calculate also the slope of those pixels which represented water edge mainly in recent years
  var slope11to20=im_reduced_sens.select('slope').updateMask(count3.gte(5));

  var sedzones_tmp=slope.where(im_reduced_sens.select('p-value').gt(0.1),0);
  var sedzones_tmp11to20=slope11to20.where(im_reduced_sens.select('p-value').gt(0.1),0);  
  sedzones_tmp=sedzones_tmp11to20.blend(sedzones_tmp);
  var sedzone_filled=sedzones_tmp.focal_mean({radius: 20, kernelType :'circle',units: 'meters',iterations: 1}).blend(sedzones_tmp);
  sedzone_filled=
    im_reduced_sens.select('slope').updateMask(count.gte(6)).focal_mean({radius: 2, kernelType :'circle',units: 'pixels',iterations: 1})
    .blend(
      im_reduced_sens.select('slope').updateMask(count.gte(6))//wo smoothing
    )
    //pixels that represent the shorelines only in less than 3 out of 22 years will be masked
    .updateMask(count.gte(3)).reproject({crs: count.projection().crs(), scale:10}).rename('sedzone_filled');  
  return im_reduced_sens.addBands(sedzone_filled).copyProperties(ee.Image(img));
};

//////////////////////////////
exports.get_landsat_images=function(aoi,displacement_L5,displacement_L7,displacement_L8,start_date,end_date, options){
  aoi=ee.Geometry(aoi);
  start_date=ee.Date(start_date);
  end_date=ee.Date(end_date);
  var cc_thresh = options.cc_thresh || 30;
  var nodata_thresh = options.nodata_thresh_l8 || 30;
  var nodata_thresh_l7 = options.nodata_thresh_l7 || 50;
  var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterBounds(aoi)
   .filterDate(start_date, end_date);//
  var withMNDWI_l8 = ee.ImageCollection(l8).map(function(img){
    return ee.Image(img).set('SENSING_TIME',ee.Date(ee.Image(img).get('system:time_start')).format('YYYY-MM-dd'));
  });
  var img_names_l8=withMNDWI_l8.sort('system:time_start').aggregate_array('SENSING_TIME');
  var img_list_l8=ee.Dictionary.fromLists(img_names_l8.distinct(), img_names_l8.distinct()).keys();
  withMNDWI_l8=ee.ImageCollection(img_list_l8.map(function(x){
    var img1=ee.Image(withMNDWI_l8.filter(ee.Filter.eq('SENSING_TIME', x)).first());
    //Mosaic the available images
    return ee.Image(withMNDWI_l8.filter(ee.Filter.eq('SENSING_TIME', x)).mosaic())
    .clip(aoi).copyProperties(img1)
    .set('system:time_start',img1.get('system:time_start'));
  })).map(function(img){ return morefunctions.cloudscore_L8_pro(ee.Image(img),aoi)}).map(morefunctions.addbands_l8); 
  withMNDWI_l8=withMNDWI_l8.filter(ee.Filter.lt('cloud_cover', cc_thresh))
    .filter(ee.Filter.lt('nodata_cover', nodata_thresh));        
  img_names_l8=ee.ImageCollection(withMNDWI_l8.sort('SENSING_TIME',false).sort('nodata_cover')).aggregate_array('SENSING_TIME');
  withMNDWI_l8=withMNDWI_l8.map(function(img){
    return ee.Image(img).displace(displacement_L8).clip(aoi);
    });
  //due to coregistration nodata fraction might increase - recalculate
  withMNDWI_l8=withMNDWI_l8.map(function(img){ return morefunctions.cloudscore_L8_pro(ee.Image(img),aoi)})
    .filter(ee.Filter.lt('nodata_cover', nodata_thresh));

  ////////NOW LANDSAT 7//////////////
  var l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterBounds(aoi)
   .filterDate(start_date, end_date);
  var withMNDWI_l7=l7.map(function(img){
    return ee.Image(img).set('SENSING_TIME',ee.Date(ee.Image(img).get('system:time_start')).format('YYYY-MM-dd'));
  });
  var img_names_l7=withMNDWI_l7.sort('system:time_start').aggregate_array('SENSING_TIME');
  var img_list_l7=ee.Dictionary.fromLists(img_names_l7.distinct(), img_names_l7.distinct()).keys();
  withMNDWI_l7=ee.ImageCollection(img_list_l7.map(function(x){
    var img1=ee.Image(withMNDWI_l7.filter(ee.Filter.eq('SENSING_TIME', x)).first());
    //Mosaic the available Landsat 7 images
    return ee.Image(withMNDWI_l7.filter(ee.Filter.eq('SENSING_TIME', x)).mosaic())
    .clip(aoi).copyProperties(img1)
    .set('system:time_start',img1.get('system:time_start'));
  })).map(function(img){ return morefunctions.cloudscore_L7_pro(ee.Image(img),aoi)}).map(morefunctions.addbands_l7); 
  withMNDWI_l7=withMNDWI_l7.filter(ee.Filter.lt('cloud_cover', cc_thresh))
    .filter(ee.Filter.lt('nodata_cover', 50));
  withMNDWI_l7=withMNDWI_l7.map(function(img){
    return ee.Image(img).displace(displacement_L7).clip(aoi);
    });
  //due to coregistration nodata fraction might increase - recalculate
  withMNDWI_l7=withMNDWI_l7.map(function(img){ return morefunctions.cloudscore_L7_pro(ee.Image(img),aoi)})
    .filter(ee.Filter.lt('nodata_cover', nodata_thresh_l7));

  ////////NOW LANDSAT 5//////////////
  var l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterBounds(aoi)
   .filterDate(start_date, end_date);
  var withMNDWI_l5=l5.map(function(img){
    return ee.Image(img).set('SENSING_TIME',ee.Date(ee.Image(img).get('system:time_start')).format('YYYY-MM-dd'));
  });
  var img_names_l5=withMNDWI_l5.sort('system:time_start').aggregate_array('SENSING_TIME');
  var img_list_l5=ee.Dictionary.fromLists(img_names_l5.distinct(), img_names_l5.distinct()).keys();
  withMNDWI_l5=ee.ImageCollection(img_list_l5.map(function(x){
    var img1=ee.Image(withMNDWI_l5.filter(ee.Filter.eq('SENSING_TIME', x)).first());
    //Mosaic the available images
    return ee.Image(withMNDWI_l5.filter(ee.Filter.eq('SENSING_TIME', x)).mosaic())
    .clip(aoi).copyProperties(img1)
    .set('system:time_start',img1.get('system:time_start'));
  })).map(function(img){ return morefunctions.cloudscore_L7_pro(ee.Image(img),aoi)}).map(morefunctions.addbands_l7); 
  withMNDWI_l5=withMNDWI_l5.filter(ee.Filter.lt('cloud_cover', cc_thresh))
    .filter(ee.Filter.lt('nodata_cover', nodata_thresh));

  withMNDWI_l5=withMNDWI_l5.map(function(img){
    return ee.Image(img).displace(displacement_L5).clip(aoi);
    });
  //due to coregistration nodata fraction might increase - recalculate
  withMNDWI_l5=withMNDWI_l5.map(function(img){ return morefunctions.cloudscore_L7_pro(ee.Image(img),aoi)})
    .filter(ee.Filter.lt('nodata_cover', nodata_thresh));
  
  ///////////
   var merged_landsat=withMNDWI_l7.merge(withMNDWI_l8).merge(withMNDWI_l5).select(['B', 'G', 'R', 'NIR','SWIR1','MNDWI','cloud','NDVI'])
      .sort('system:time_start');
  return merged_landsat;
};

exports.coreg_getdisplacement = function(date,sat,img,bounds){
  bounds=ee.Geometry(bounds);
  date=ee.Date(date);
  var sensor_name=ee.String(sat);
  var imgS2=ee.Image(img);
  var bandIn_l8 = ['B2','B3','B4','B5','B6','B7'];
  var bandIn_l7 = ['B1','B2','B3','B4','B5','B7'];
  var bandIn_s2 = ['B2','B3','B4','B8','B11','B12'];
  var bandOut = ['blue','green','red','nir','swir1','swir2'];
  var sensor_list= ee.List(['LANDSAT_5','LANDSAT_7','LANDSAT_8']);
  var proj=imgS2.select('green').projection();

  var imgL5col = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterBounds(bounds)
         .filterDate(date,date.advance(1,'day')).select(bandIn_l7,bandOut);
  var imgl5a=ee.Image(imgL5col.mosaic().copyProperties(imgL5col.first()).set('system:time_start',imgL5col.first().get('system:time_start')))
    .clip(bounds).setDefaultProjection(proj);
  var imgl5b=imgL5col.first();
  var imgL5= ee.Image(ee.Algorithms.If(imgL5col.size().gt(1),imgl5a,imgl5b));

  var imgL7col = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterBounds(bounds)
         .filterDate(date,date.advance(1,'day')).select(bandIn_l7,bandOut);
  var imgL7a=ee.Image(imgL7col.mosaic().copyProperties(imgL7col.first()).set('system:time_start',imgL7col.first().get('system:time_start')))
    .clip(bounds).setDefaultProjection(proj);
  var imgL7b=imgL7col.first();
  var imgL7= ee.Image(ee.Algorithms.If(imgL7col.size().gt(1),imgL7a,imgL7b));
  
  var imgL8col = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterBounds(bounds)
         .filterDate(date,date.advance(1,'day')).select(bandIn_l8,bandOut);
  var imgL8a=ee.Image(imgL8col.mosaic().copyProperties(imgL8col.first()).set('system:time_start',imgL8col.first().get('system:time_start')))
    .clip(bounds).setDefaultProjection(proj);
  var imgL8b=imgL8col.first();
  var imgL8= ee.Image(ee.Algorithms.If(imgL8col.size().gt(1),imgL8a,imgL8b));
  
  var imgL=ee.Image(ee.List([imgL5,imgL7,imgL8]).get(sensor_list.indexOf(sensor_name)));
  
  var LRedBand = imgL.select('green');
  var S2RedBand = imgS2.select('green');

  var displacement = LRedBand.displacement({
    referenceImage: S2RedBand,
    maxOffset: 50.0,//The maximum offset allowed when attempting to align the input images, in meters
    patchWidth: 100 // Small enough to capture texture and large enough that ignorable 
  }).reproject(ee.Projection('EPSG:4326'), null, 30);
  return displacement;
};
