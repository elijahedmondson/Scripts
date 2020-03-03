setChannelNames(
      'DAPI',
     'CD3',
     'Cyotkeratin',
     'IFNg'
)
selectAnnotations();

runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 100.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 5.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

//Create "single measurement classifier"

runClassifier('C://Users//edmondsonef//Desktop//QuPath//Dong Zhang//ER IF//classifiers//object_classifiers//CD3.json');


// Load and run a trained classifier
def pathClassifier = buildFilePath(PROJECT_BASE_DIR , 'classifiers', 'CD3.json')
runClassifier(pathClassifier);