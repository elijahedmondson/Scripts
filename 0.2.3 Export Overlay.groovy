import qupath.lib.images.servers.LabeledImageServer  
def imageData = getCurrentImageData()

// Define output path (relative to project)
def outputDir = buildFilePath(PROJECT_BASE_DIR, 'export')
mkdirs(outputDir)
print '01. Directory = ' + (outputDir)
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def Path = buildFilePath(outputDir, name + ".tif")
print '02. Filepath = ' + (Path)
def maskPath = buildFilePath(outputDir, name + "-mask.tif")
print '03. Masked Filepath = ' + (maskPath)



// Define output resolution
double requestedPixelSize = 20.0

// Convert to downsample
double downsample = requestedPixelSize / imageData.getServer().getPixelCalibration().getAveragedPixelSize()
print (downsample)

// Create an ImageServer where the pixels are derived from annotations
def labelServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.WHITE) // Specify background label (usually 0 or 255)
    .downsample(downsample)    // Choose server resolution; this should match the resolution at which tiles are exported
    .useUniqueLabels()  
    .useFilter({p -> p.getPathClass() == getPathClass('Tumor')})
    .lineThickness(2)          // Optionally export annotation boundaries with another label
    //.setBoundaryLabel('Boundary*', 4) // Define annotation boundary label
    .multichannelOutput(false) // If true, each label refers to the channel of a multichannel binary image (required for multiclass probability)
    .build()

writeImage(labelServer, Path)