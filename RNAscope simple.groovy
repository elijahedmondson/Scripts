setImageType('BRIGHTFIELD_H_DAB');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759 ", "Background" : " 255 255 255 "}');


runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2', '{"threshold": 230,  "requestedPixelSizeMicrons": 20.0,  "minAreaMicrons": 10000.0,  "maxHoleAreaMicrons": 1000000.0,  "darkBackground": false,  "smoothImage": true,  "medianCleanup": true,  "dilateBoundaries": false,  "smoothCoordinates": true,  "excludeOnBoundary": false,  "singleAnnotation": true}');
selectAnnotations();

runPlugin('qupath.imagej.detect.nuclei.WatershedCellDetection', '{"detectionImageBrightfield": "Hematoxylin OD",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 20.0,  "maxAreaMicrons": 400.0,  "threshold": 0.05,  "maxBackground": 2.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 10.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[DAB]": 0.4,  "doSmoothing": false,  "splitByIntensity": true,  "splitByShape": false,  "spotSizeMicrons": .3,  "minSpotSizeMicrons": 0.01,  "maxSpotSizeMicrons": 1.5,  "includeClusters": false}');
setCellIntensityClassifications("Subcellular: DAB: Num spots estimated", 1, 4, 10)


/**
 * Export a thumbnail image, with and without an overlay, using QuPath.
 *
 * For tissue microarrays, the scripting code written by the 'File -> Export TMA data'
 * command is probably more appropriate.
 *
 * However, for all other kinds of images where batch export is needed this script can be used.
 *
 * @author Pete Bankhead
 */


import qupath.lib.gui.ImageWriterTools
import qupath.lib.gui.QuPathGUI
import qupath.lib.gui.viewer.OverlayOptions
import qupath.lib.regions.RegionRequest
import qupath.lib.scripting.QPEx

// Aim for an output resolution of approx 20 ?m/pixel
double requestedPixelSize = 10

// Create the output directory, if required
def path = QPEx.buildFilePath(QPEx.PROJECT_BASE_DIR, "Image masks")
QPEx.mkdirs(path)

// Get the imageData & server
def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()

// Get the file name from the current server
def name = server.getShortServerName()

// We need to get the display settings (colors, line thicknesses, opacity etc.) from the current viewer, if available
def overlayOptions = QuPathGUI.getInstance() == null ? new OverlayOptions() : QuPathGUI.getInstance().getViewer().getOverlayOptions()

// Calculate downsample factor depending on the requested pixel size
double downsample = requestedPixelSize / server.getAveragedPixelSizeMicrons()
def request = RegionRequest.createInstance(imageData.getServerPath(), downsample, 0, 0, server.getWidth(), server.getHeight())

// Write output image, with and without overlay
def dir = new File(path)
def fileImage = new File(dir, name + ".jpg")
def img = ImageWriterTools.writeImageRegion(server, request, fileImage.getAbsolutePath())
def fileImageWithOverlay = new File(dir, name + "-overlay.jpg")
ImageWriterTools.writeImageRegionWithOverlay(img, imageData, overlayOptions, request, fileImageWithOverlay.getAbsolutePath())
print("Done")


//SAVE ANNOTATIONS //

def name = getProjectEntry().getImageName() + '.txt'
//def path = buildFilePath('C:/Users/edmondsonef/Desktop/QuPath/Venditti/', 'RNAscope quantification')
def path = buildFilePath(PROJECT_BASE_DIR, 'RNAscope quantification')
mkdirs(path)
path = buildFilePath(path, name)
saveAnnotationMeasurements(path)