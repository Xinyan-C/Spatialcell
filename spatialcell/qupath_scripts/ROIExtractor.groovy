/**
 * This script reads rectangular ROI coordinates from a text file, creates corresponding annotations in QuPath, and adds them to the image hierarchy, with support for up to three rectangles
 */

import qupath.lib.roi.RectangleROI
import qupath.lib.objects.PathAnnotationObject
import java.nio.file.Files
import java.nio.file.Paths

// 📌 1. Read ROI file
def filePath = "/path/to/E18.5_ranges.txt"
def lines = Files.readAllLines(Paths.get(filePath))

// Store parsed rectangle coordinates
def rectangles = []
def currentLabel = null

// 📌 2. Parse coordinate information
for (line in lines) {
    line = line.trim()
    // If the line ends with "Rectangle coordinate range:", it indicates a new label
    if (line.endsWith("Rectangle coordinate range:")) {
        currentLabel = line.replace("Rectangle coordinate range:", "").trim()
    } 
    // If there is a current label and the line starts with "X:" or "Y:" for coordinate information
    else if (currentLabel && (line.startsWith("X:") || line.startsWith("Y:"))) {
        // Skip if the name contains "Entire image"
        if (currentLabel.contains("Entire image")) {
            continue
        }

        def values = line.replace("X:", "")
                         .replace("Y:", "")
                         .trim()
                         .split(" - ")
        if (values.size() == 2) {
            def minVal = values[0].toDouble()
            def maxVal = values[1].toDouble()
            if (line.startsWith("X:")) {
                rectangles << [label: currentLabel, x1: minVal, x2: maxVal]
            } else if (line.startsWith("Y:")) {
                rectangles.last().y1 = minVal
                rectangles.last().y2 = maxVal
            }
        }
    }
}

// Keep only the first 3 rectangles
if (rectangles.size() > 3) {
    rectangles = rectangles.take(3)
}

// 📌 3. Create annotations and output coordinates
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()

for (rect in rectangles) {
    def x = rect.x1
    def y = rect.y1
    def width = rect.x2 - rect.x1
    def height = rect.y2 - rect.y1
    
    def roi = new RectangleROI(x, y, width, height)
    def annotation = new PathAnnotationObject(roi)
    annotation.setName(rect.label)
    hierarchy.addPathObject(annotation)
    
    // Output annotation coordinates and dimensions
    println("Annotation ${rect.label}:")
    println("  File coordinates: x=${rect.x1}, y=${rect.y1}, width=${width}, height=${height}")
    println("  QuPath coordinates: x=${roi.getBoundsX()}, y=${roi.getBoundsY()}, width=${roi.getBoundsWidth()}, height=${roi.getBoundsHeight()}")
}

// 📌 4. Refresh the view
fireHierarchyUpdate()
println("✅ Successfully added ${rectangles.size()} rectangle annotations!")

// Get image dimensions
def server = imageData.getServer()
def imageWidth = server.getWidth()
def imageHeight = server.getHeight()
println("Image dimensions: ${imageWidth} x ${imageHeight}")