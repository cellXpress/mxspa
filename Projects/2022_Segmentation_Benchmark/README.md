# 2022 Tissue segmentation benchmark

Python scripts to compare the segmentation results from cX2 and various other
tools.


# Installing the reference data sets

1. Assuming the project root directory is `/home/user/projects`
   (You may substitute it with your own root directory location).

2. Clone the repo using Git under the project root directory.
   ```
   git clone https://github.com/ccpa/mxspa /home/user/projects
   ```

3. Download and install cellXpress 2 (cX2) from http://www.cellXpress.org

4. Run cX2 and make sure the project root path is point to
   `/home/user/projects`.

5. Download the four reference cX2 projects and images.


# Runnning the benchmark

1. Install required library

   ```
   pip3 install --user matplotlib scipy h5py numpy roifile opencv-python tiffile cmapy pandas imagecodecs
   ```

2. Update the location of benchmark data at segmentation_benchmark.py

3. Run the benchmark comparison code
   ```
   python segmentation_benchmark.py
   ```

4. The segmentation boundaries are saved under `figures`, and the raw
   performance values are saved under `results`