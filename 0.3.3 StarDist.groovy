import qupath.ext.stardist.StarDist2D
var pathModel = 'F:/QuPath/Stardist/he_heavy_augment.pb'

selectAnnotations();

var stardist = StarDist2D.builder(pathModel)
      .ignoreCellOverlaps(false)   // Set to true if you don't care if cells expand into one another
      .threshold(0.5)              // Prediction threshold
      .normalizePercentiles(1, 99) // Percentile normalization
      .pixelSize(0.5)              // Resolution for detection
      .includeProbability(true)    // Include prediction probability as measurement
      .cellExpansion(0.0)          // Approximate cells based upon nucleus expansion
      .cellConstrainScale(3)       // Constrain cell expansion using nucleus size
      .measureShape()              // Add shape measurements
      .measureIntensity()          // Add cell measurements (in all compartments)
      .doLog()                     // Use this to log a bit more information while running the script
      .build()

// Run detection for the selected objects
var imageData = getCurrentImageData()
var pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println 'Done!'

setDetectionIntensityClassifications("DAB: Mean", 0.13, 0.5, 0.9)