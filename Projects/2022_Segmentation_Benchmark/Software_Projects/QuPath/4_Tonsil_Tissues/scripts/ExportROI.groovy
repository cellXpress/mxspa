// Export results as ROI
import ij.plugin.frame.RoiManager

def imageData = getCurrentImageData()

def img_name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def output_name = img_name + "_Cell.zip"
    
def path = buildFilePath(PROJECT_BASE_DIR, "rois", output_name)

def annotations = getAnnotationObjects()
def roiMan = new RoiManager(false)
double x = 0
double y = 0
double downsample = 1 // Increase if you want to export to work at a lower resolution
annotations.each {
  def roi = IJTools.convertToIJRoi(it.getROI(), x, y, downsample)
  roiMan.addRoi(roi)
}
roiMan.runCommand("Save", path)