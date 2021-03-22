setChannelNames(
      'DAPI',
     'SynMUT1')
     
selectAnnotations();

//Cytoplasm expansion
//runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.2,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 5000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

//Nucleus only
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.2,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 5000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.01,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1.0,  "detection[Channel 2]": 30000.0,  "doSmoothing": true,  "splitByIntensity": true,  "splitByShape": true,  "spotSizeMicrons": 1.0,  "minSpotSizeMicrons": 0.2,  "maxSpotSizeMicrons": 6.0,  "includeClusters": true}');

setCellIntensityClassifications("Subcellular: Channel 2: Num spots estimated", 1, 10, 20)

//Using Positive Cell Detection
//"MEAN" signal
//runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 5000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true,  "thresholdCompartment": "Cell: SynMUT1 mean",  "thresholdPositive1": 8000.0,  "thresholdPositive2": 12000.0,  "thresholdPositive3": 15000.0,  "singleThreshold": false}');

//Using Positive Cell Detection
//"MAX" signal
//runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 6000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true,  "thresholdCompartment": "Cell: SynMUT1 max",  "thresholdPositive1": 30000.0,  "thresholdPositive2": 40000.0,  "thresholdPositive3": 50000.0,  "singleThreshold": false}');
