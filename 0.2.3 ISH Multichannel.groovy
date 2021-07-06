setChannelNames(
      'DAPI',
     'TLR2',
     'CD206')
     
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.1,  "backgroundRadiusMicrons": 0.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 8000.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 7,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1,  "detection[Channel 2]": 21000.0,  "detection[Channel 3]": 8000.0,  "doSmoothing": false,  "splitByIntensity": true,  "splitByShape": true,  "spotSizeMicrons": 1.0,  "minSpotSizeMicrons": 0.1,  "maxSpotSizeMicrons": 100.0,  "includeClusters": true}');

setCellIntensityClassifications("Subcellular: Channel 2: Num spots estimated", 1, 10, 20)




// From Pete's post on Gitter, another way of applying cell classifications
//0.1.2 and 0.2.0 (though channel names have changed in 0.2.0)
// Get cells & reset all the classifications
def cells = getCellObjects()
resetDetectionClassifications()

cells.each {it.setPathClass(getPathClass('Negative'))}

// Get channel 1 & 2 positives
def ch1Pos = cells.findAll {measurement(it, "Subcellular: Channel 2: Num spots estimated") > 2}
ch1Pos.each {it.setPathClass(getPathClass('TLR2'))}

def ch2Pos = cells.findAll {measurement(it, "Subcellular: Channel 3: Num spots estimated") > 2}
ch2Pos.each {it.setPathClass(getPathClass('CD206'))}

// Overwrite classifications for double positives
def doublePos = ch1Pos.intersect(ch2Pos)
doublePos.each {it.setPathClass(getPathClass('Double positive'))}

fireHierarchyUpdate()