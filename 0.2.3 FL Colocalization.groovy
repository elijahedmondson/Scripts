//setChannelNames(
//      'DAPI',
//     'GFP-488',
//     'NeuN-594')
     
selectAnnotations();

runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 0.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 100.0,  "threshold": 3000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

runObjectClassifier("CD8.IFN")