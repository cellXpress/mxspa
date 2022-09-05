CellProfiler Benchmark <!-- omit in toc -->
======================

- [Setup and requirements](#setup-and-requirements)
- [Disabling Windows Search](#disabling-windows-search)
- [Running benchmark](#running-benchmark)

Setup and requirements
======================
1. Download and install the software from https://cellprofiler.org.
   We tested version 4.2.4 on MS Windows 11.

2. Start CellProfiler.

3. Under `Output Settings`,
   * Update `Default Input Folder` to the cX2 Benchmark
   code main folder (e.g., H:\Research\2022_07_cX2_benchmark).
   * Update `Default Output Folder` to the cX2 Benchmark code
   CellProfiler sub-folder (e.g., H:\Research\2022_07_cX2_benchmark\CellProfiler\results).

4. Under `File > Preferences`,  update the maximum number of workers to the
   number of physical cores in the current CPU.

5. Download and install all four cX2 reference projects and images directly 
   using cX2.

Disabling Windows Search
========================
1. For Windows, `Windows Search` indexer may affect the performance 
   during benchmarking, and thus need to be disabled.

2. Run the `Services` management console as an administrator.

3. Find `Windows Search`, righ-click on it, and select `Properties`.

4. Change `Startup type` to 'Disabled', and click `Apply`.

5. Press `Stop` to stop the service.

6. The `Service status` should be 'Stopped'.

7. Disable all the anti-virus realtime protection.


Running benchmark
=================
1. Start CellProfiler.

2. Open one of the CellProfiler projects:
   * CellProfiler\1_BEAS2B_Lung_Cells.cpproj
   * CellProfiler\2_COVID_Lung_Tissues.cpproj
   * CellProfiler\3_Lung_Tumors.cpproj
   * CellProfiler\4_Tonsil_Tissues.cpproj

3. Click on the `Images` module and select `Clear File List`.

4. Drag the correponding raw image folder into the 
   `Images` module:
   * 2022_07_cX2_BEAS2B_Lung_Cells\BEAS2B_Lung_Cell_Images-JY15135
   * 2022_07_cX2_COVID_Lung_Tissues\COVID19_Lung_Tissue_Images-JY22001
   * 2022_07_cX2_Lung_Tumors\Lung_Tumors_Images-JY22010
   * 2022_07_cX2_Tonsil_Tissues\Tonsil_Tissues-JY22002

5. and make sure 140 rows have been found.

6. There is a bug in CellPrfoiler (v4.2.4) that prevents it from loading
   the image metadata properly. To overcome this, click on the `Metadata` module, 
   remove the `Extract from image file headers` method by clicking on the 
   `Remove this extraction method` button.

7. Click on the `Update` button to refresh the metadata.
   
8. Click on `Add another extraction method`, select `Extract from image   
   file headers`, and click `Extract metadata`.

9.  Click on the `Update` button again to refresh the metadata. Now, the number
   of rows found should be correct.

11. Click on the `NamesAndTypes` module, and click on the `Update` button
    to detect all the markers. The number of image sets found should be 140.

12. Click on the `Start Analysis` button to start the analysis.
    
13. Once the analysis is done, copy the final results to `results_final`.