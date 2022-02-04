import qupath.ext.stardist.StarDist2D
//setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759", "Background" : " 228 223 226"}');

// Specify the model file (you will need to change this!)
var pathModel = 'F:/QuPath/Stardist/he_heavy_augment.pb'

// Get current image - assumed to have color deconvolution stains set
println '1'
selectAnnotations();

var stardist = StarDist2D.builder(pathModel)
      .ignoreCellOverlaps(false)   // Set to true if you don't care if cells expand into one another
      .threshold(0.4)              // Prediction threshold
      .normalizePercentiles(1, 99) // Percentile normalization
      .pixelSize(0.2)              // Resolution for detection
      .includeProbability(true)    // Include prediction probability as measurement
      .cellExpansion(0.0)          // Approximate cells based upon nucleus expansion
      .cellConstrainScale(3)       // Constrain cell expansion using nucleus size
      .measureShape()              // Add shape measurements
      //.measureIntensity()          // Add cell measurements (in all compartments)
      .doLog()                     // Use this to log a bit more information while running the script
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

//setDetectionIntensityClassifications("DAB: Mean", 0.13, 0.5, 0.9)

println 'DONE'

//Additional options
//var stardist = StarDist2D.builder(pathModel)
//        .threshold(0.5)              // Probability (detection) threshold
//        .channels('DAPI')            // Select detection channel
//        .normalizePercentiles(1, 99) // Percentile normalization
//        .pixelSize(0.5)              // Resolution for detection
//        .tileSize(1024)              // Specify width & height of the tile used for prediction
//        .cellExpansion(5.0)          // Approximate cells based upon nucleus expansion
//        .cellConstrainScale(1.5)     // Constrain cell expansion using nucleus size
//        .ignoreCellOverlaps(false)   // Set to true if you don't care if cells expand into one another
//        .measureShape()              // Add shape measurements
//        .measureIntensity()          // Add cell measurements (in all compartments)
//        .includeProbability(true)    // Add probability as a measurement (enables later filtering)
//        .nThreads(4)                 // Limit the number of threads used for (possibly parallel) processing
//        .simplify(1)                 // Control how polygons are 'simplified' to remove unnecessary vertices
//        .doLog()                     // Use this to log a bit more information while running the script
//        .createAnnotations()         // Generate annotation objects using StarDist, rather than detection objects
//        .constrainToParent(false)    // Prevent nuclei/cells expanding beyond any parent annotations (default is true)
//        .classify("Tumor")           // Automatically assign all created objects as 'Tumor'
//        .build()


//var imageData = getCurrentImageData()
//var stains = imageData.getColorDeconvolutionStains()
//      .preprocess( // Extra preprocessing steps, applied sequentially
//            ImageOps.Channels.deconvolve(stains),
//            ImageOps.Channels.extract(0),
//            ImageOps.Filters.median(2),
//            ImageOps.Core.divide(1.5)
//       )