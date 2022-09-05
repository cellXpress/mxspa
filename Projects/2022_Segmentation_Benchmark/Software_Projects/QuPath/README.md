QuPath Benchmark <!-- omit in toc -->
======================

Setup and requirements
======================
1. Download and install the software from https://qupath.github.io/.
   We tested version 0.3.2 on MS Windows 11.

2. Download and install QuPath StarDist extension follwing the instructions
   from https://github.com/qupath/qupath-extension-stardist.

3. Update `Preferences > General > Number of parallel threads` to 
   the number of physical cores in the local processor.

Running benchmark
=================
1. Start QuPath.

2. Open one of the QuPath projects:
   * QuPath\BEAS2B_Lung_Cells\project.qpproj

3. If QuPath complains about missing image files, click on `Search` and select
   the appropriate local raw image folder:
   * 2022_07_cX2_BEAS2B_Lung_Cells\BEAS2B_Lung_Cell_Images-JY15135
   * 2022_07_cX2_COVID_Lung_Tissues\COVID19_Lung_Tissue_Images-JY22001
   * 2022_07_cX2_Lung_Tumors\Lung_Tumors_Images-JY22010
   * 2022_07_cX2_Tonsil_Tissues\Tonsil_Tissues-JY22002
   
   Press `Apply Changes` to continue. Please wait patiently while QuPath is
   reading the image file metadata.

4. If QuPath fails to find the files, click `Cancel` and add the images manually.

5. Click on `Annotations` and make sure the channel names are correct. If not,
   click on `...` > `Populate from image channels` and don't keep existing 
   classes.

6. Open the script editor via `Automate` > `Project Scripts` >
   `StarDistSegmentation`

7. Update the number of threads to the the number of physical cores 
   in the local processor.
   ```
   .nThreads(8)
   ``` 

8. Click on the `Brightness & contrast` button, and make sure that
   the correct DNA channel name is used for segmentation:
   ```
   .channels('DNA') 
   ```

9.  Select images to analyze by clicking `Run` > `Run for project`.
   All images in the project should be selected.

11. Press `OK` to start the analysis. The script starting and ending times are
    logged in the editor console.

12. The results are saved under `rois`. After all the images have been analyzed,
    move the final results for the selected images to `rois_final`.