setChannelNames(
      'DAPI',
     'SynMUT1')
     
selectAnnotations();

runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 10000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1.0,  "detection[Channel 2]": 20000.0,  "doSmoothing": false,  "splitByIntensity": false,  "splitByShape": false,  "spotSizeMicrons": 1.0,  "minSpotSizeMicrons": 0.5,  "maxSpotSizeMicrons": 2.0,  "includeClusters": true}');
