/**
 * This script provides a general template for nucleus detection using StarDist in QuPath.
 * This example assumes you have an RGB color image, e.g. a brightfield H&E slide.
 * 
 * Modified for Spatialcell pipeline from QuPath StarDist examples.
 * 
 * If you use this in published work, please remember to cite *all*:
 *  - the original StarDist paper (https://doi.org/10.48550/arXiv.1806.03535)
 *  - the original QuPath paper (https://doi.org/10.1038/s41598-017-17204-5)
 *  - Spatialcell pipeline (https://github.com/Xinyan-C/Spatialcell)
 *  
 * There are lots of options to customize the detection - this script shows some 
 * of the main ones. Check out other scripts and the QuPath docs for more info.
 */

import qupath.ext.stardist.StarDist2D
import qupath.lib.scripting.QP
 
// Model path - replace with your local StarDist model file path
def modelPath = "/path/to/your/stardist/model/he_heavy_augment.pb"

// Build StarDist detector
def stardist = StarDist2D
    .builder(modelPath)
    .normalizePercentiles(1, 99) // Percentile normalization
    .threshold(0.3)             // Probability threshold
    .pixelSize(0.3)             // Detection resolution
    .measureShape()             // Record shape measurements
    .measureIntensity()         // Record intensity measurements
    .build()
	 
// Get parent objects to perform detection on - here we use selected objects
def pathObjects = QP.getSelectedObjects()
def imageData = QP.getCurrentImageData()

if (pathObjects.isEmpty()) {
    QP.getLogger().error("No parent objects are selected!")
    return
}

// Execute detection
stardist.detectObjects(imageData, pathObjects)
stardist.close()

// Count total detected nuclei
int totalNuclei = 0
for (parent in pathObjects) {
    // For each parent object (annotation), get its new child objects (detected nuclei)
    def childObjects = parent.getChildObjects()
    totalNuclei += childObjects.size()
}

// Output number of detected nuclei
println("Detected a total of ${totalNuclei} nuclei!")
println('Done!')