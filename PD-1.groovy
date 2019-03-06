setImageType('BRIGHTFIELD_H_DAB');
selectAnnotations();
runPlugin('qupath.imagej.detect.nuclei.PositiveCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 0.08,  "maxBackground": 2.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 1.5,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true,  "thresholdCompartment": "Cell: DAB OD mean",  "thresholdPositive1": 0.2,  "thresholdPositive2": 0.3,  "thresholdPositive3": 0.6,  "singleThreshold": true}');



//SAVE ANNOTATIONS //

def name = getProjectEntry().getImageName() + '.txt'
def path = buildFilePath('C:/Users/edmondsonef/Desktop/QuPath/Dong Zhang/ER PD1/', 'PD1')
mkdirs(path)
path = buildFilePath(path, name)
saveAnnotationMeasurements(path)