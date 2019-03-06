setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049 ", "Stain 2" : "Eosin", "Values 2" : "0.2159 0.8012 0.5581 ", "Background" : " 255 255 255 "}');
selectDetections();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "colorOD": false,  "colorStain1": false,  "colorStain2": false,  "colorStain3": false,  "colorRed": false,  "colorGreen": false,  "colorBlue": false,  "colorHue": false,  "colorSaturation": false,  "colorBrightness": false,  "doMean": true,  "doStdDev": true,  "doMinMax": true,  "doMedian": true,  "doHaralick": true,  "haralickDistance": 1,  "haralickBins": 32}');
selectCells();
runPlugin('qupath.lib.plugins.objects.ShapeFeaturesPlugin', '{"area": true,  "perimeter": true,  "circularity": true,  "useMicrons": true}');
runClassifier('C:\\Users\\edmondsonef\\Desktop\\QuPath\\OMAL\\190482\\classifiers\\necrosis.qpclassifier');