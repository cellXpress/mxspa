////////////////////////////////////////////////////////////////////////////////////////////
/// Perform StartDist segmentation and save the results as tiff
/// Written by Lit-Hsin Loo (2022)
////////////////////////////////////////////////////////////////////////////////////////////
import qupath.ext.stardist.StarDist2D
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane

boolean touchingBoundary(server, roi) {
    return roi.getBoundsX() <= 0 ||
           roi.getBoundsY() <= 0 ||
           roi.getBoundsX() + roi.getBoundsWidth() >= server.getWidth()-1 ||
           roi.getBoundsY() + roi.getBoundsHeight() >= server.getHeight()-1
}

// Specify the model .pb file (you will need to change this!)
def pathModel = buildFilePath(PROJECT_BASE_DIR, "models/dsb2018_heavy_augment.pb")

def stardist = StarDist2D.builder(pathModel)
        .threshold(0.5)              // Probability (detection) threshold
        .channels('reg001_cyc025_ch001_DAPI-25.tif')             // Select detection channel
        .normalizePercentiles(1, 99.9) // Percentile normalization
        .pixelSize(0.5)                // Resolution for detection
        .tileSize(1024)              // Specify width & height of the tile used for prediction
        .cellExpansion(5.0)          // Approximate cells based upon nucleus expansion
        .cellConstrainScale(1.5)     // Constrain cell expansion using nucleus size
        .measureShape()              // Add shape measurements
        //.measureIntensity()          // Add cell measurements (in all compartments)
        .includeProbability(true)    // Add probability as a measurement (enables later filtering)
        .nThreads(8)
        .doLog()                     // Use this to log a bit more information while running the script
        //.createAnnotations()         // Generate annotation objects using StarDist, rather than detection objects
        .build()

// Run detection for the selected objects
def imageData = getCurrentImageData()
def server = getCurrentServer()

// Remove previous annotations
def curAnnotation = getAnnotationObjects().findAll{it.getLevel() == 1}
removeObjects(curAnnotation, true)

createSelectAllObject(true);

def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}

stardist.detectObjects(imageData, pathObjects)

// Remove big cells
int sizeBefore = getCellObjects().size()
def toDelete = getCellObjects().findAll{
    measurement(it, 'Cell: Area µm^2') < 15 ||
    measurement(it, 'Cell: Area µm^2') > 500 ||
    touchingBoundary(server, it.getROI())
}
removeObjects(toDelete, true)

int sizeAfter = getCellObjects().size()
int sizeRemoved = sizeBefore - sizeAfter
println(sprintf("Before: %d, After: %d, Removed: %d", sizeBefore, sizeAfter, sizeRemoved))

// Define output path (relative to project)
def outputDir = buildFilePath(PROJECT_BASE_DIR, 'rois')
mkdirs(outputDir)

def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def path = buildFilePath(outputDir, name + "-labels.ome.tif")

// Define how much to downsample during export (may be required for large images)
double downsample = 1

// Create an ImageServer where the pixels are derived from annotations
def labelServer = new LabeledImageServer.Builder(imageData)
  .backgroundLabel(0, ColorTools.WHITE) // Specify background label (usually 0 or 255)
  .downsample(downsample)    // Choose server resolution; this should match the resolution at which tiles are exported
  .useCells()
  .useInstanceLabels()
  .multichannelOutput(false) // If true, each label refers to the channel of a multichannel binary image (required for multiclass probability)
  .build()

// Write the image
writeImage(labelServer, path)



println 'Done!'