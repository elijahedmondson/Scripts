import qupath.ext.stardist.StarDist2D
//setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759", "Background" : " 228 223 226"}');

// Specify the model file (you will need to change this!)
var pathModel = 'F:/QuPath/Stardist/he_heavy_augment.pb'

// Get current image - assumed to have color deconvolution stains set
println '1'
//runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2', '{"threshold": 230,  "requestedPixelSizeMicrons": 20.0,  "minAreaMicrons": 1000000.0,  "maxHoleAreaMicrons": 5000.0,  "darkBackground": false,  "smoothImage": true,  "medianCleanup": true,  "dilateBoundaries": false,  "smoothCoordinates": true,  "excludeOnBoundary": false,  "singleAnnotation": true}');

selectAnnotations();

var stardist = StarDist2D.builder(pathModel)
      .ignoreCellOverlaps(false)   // Set to true if you don't care if cells expand into one another
      .threshold(0.2)              // Prediction threshold
      .normalizePercentiles(2, 99) // Percentile normalization
      .pixelSize(0.9)              // Resolution for detection
      //.includeProbability(true)    // Include prediction probability as measurement
      .cellExpansion(7.0)          // Approximate cells based upon nucleus expansion
      .cellConstrainScale(3)       // Constrain cell expansion using nucleus size
      //.measureShape()              // Add shape measurements
      //.measureIntensity()          // Add cell measurements (in all compartments)
      //.doLog()                     // Use this to log a bit more information while running the script
      .build()
println '2'
// Run detection for the selected objects
var imageData = getCurrentImageData()
var pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println '3'


//IHC//
//IHC//
//IHC//
//IHC//
//IHC//
//setDetectionIntensityClassifications("DAB: Mean", 0.13, 0.5, 0.9)



//RNASCOPE//
//RNASCOPE//
//RNASCOPE//
//RNASCOPE//
runPlugin('qupath.imagej.detect.nuclei.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 20.0,  "maxAreaMicrons": 400.0,  "threshold": 0.05,  "maxBackground": 2.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 10.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[DAB]": 0.3,  "doSmoothing": false,  "splitByIntensity": true,  "splitByShape": false,  "spotSizeMicrons": .2,  "minSpotSizeMicrons": 0.01,  "maxSpotSizeMicrons": 1.5,  "includeClusters": false}');
setCellIntensityClassifications("Subcellular: DAB: Num spots estimated", 1, 4, 10)

