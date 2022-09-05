cellXpress2 Benchmark <!-- omit in toc -->
=====================


Setup and requirements
======================
1. Download the reference projects direct from cX2 via `Project` >
   `Download a reference project`.

2. Set the number of threads via `Options` > `Parallel Processing` to the
   number of physical cores in the local processor.

Running benchmark for the normal cX2 projects
=============================================
1. BEAS-2B Lung cells and COVID Lung tissues are normal projects.

2. Both projects have two marker sets. One is using CellShape AI (CSAI), and
   another one is using watershed algorithm for cell segmentation.

3. Select the intended marker set.

4. Run cell segmentation for the marker set by clicking on `Profiling` > 
   `Segment cells for the current marker set`.

5. The results are saved as HDF5 under cX2 project folder
   ```
   2022_07_cX2_BEAS2B_Lung_Cells/220726_BEAS2B_Lung_Cells_Project/markerset0/segmentation/0/seg28.cXd
   ```

Running benchmark for the merged cX2 projects
=============================================
1. Lung tumors and Tonsil tissues are merged projects.

2. Each tissue type has two cX2 projects, one is using CellShape AI (CSAI) and
   another one is using watershed algorithm for cell segmentation.

   | Project name                     | Segmentation  |Marker set type |
   |----------------------------------|---------------|----------------|   
   |Lung tumors - 7 markers (CSAI)    | CellShape AI  |Merged          |
   |Lung tumors - 7 markers (WS)      | Watershed     |Single          | 
   |Tonsil tissues - 49 markers (WS)  | CellShape AI  |Merged          |
   |Tonsil tissues - 49 markers (CSAI)| Watershed     |Single          |

3. Open the intended cX2 project via `Project` > `Open an existing project`.

4. Run cell segmentation for the marker set by clicking on `Profiling` > 
   `Segment cells for all marker sets`.

5. To export the multi-market-set cell masks, open the image in Cell Viewer,
   and then under `Image` > `Save merged cell masks as...`
   