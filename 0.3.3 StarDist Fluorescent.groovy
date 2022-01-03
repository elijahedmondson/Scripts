
//setChannelNames('DAPI', 'GFP','CD45R')
import qupath.ext.stardist.StarDist2D

// Specify the model file (you will need to change this!)
println '1'
selectAnnotations();
var pathModel = 'C:/Users/edmondsonef/QuPath/Stardist Trained Models/dsb2018_heavy_augment.pb'

var stardist = StarDist2D.builder(pathModel)
        .threshold(0.2)              // Probability (detection) threshold
        .channels('DAPI')            // Specify detection channel
        .normalizePercentiles(1, 99) // Percentile normalization
        .pixelSize(0.3)              // Resolution for detection
        .cellExpansion(6.0)          // Approximate cells based upon nucleus expansion
        .cellConstrainScale(1.5)     // Constrain cell expansion using nucleus size
        .measureShape()              // Add shape measurements
        .measureIntensity()          // Add cell measurements (in all compartments)
        .includeProbability(true)    // Include prediction probability as measurement
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


//runObjectClassifier("GFP.CD45R")