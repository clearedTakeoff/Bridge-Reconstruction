# Reconstruction of bridges in aerial LIDAR point cloud data
**Advanced computer graphics**

**Faculty of Computer and Information Science, University of Ljubljana (2019 / 2020)**

Program that takes a point cloud (.laz file) and shapefile with information about bridges 
as an input and tries to reconstruct the missing parts of bridges inside the point cloud. 

Lidar data can be accessed at [eVode](http://gis.arso.gov.si/evode/profile.aspx?id=atlas_voda_Lidar@Arso&culture=en-US)
while shapefile data is available with a free account at [EGP: e-Geodetski podatki](https://egp.gu.gov.si/egp/).

The final implementation contains 4 files: one for [manipulating shapefiles](https://github.com/nzupancic/Bridge-Reconstruction/blob/master/read_shapefile.py), another 
one for [manipulation of laz files](https://github.com/nzupancic/Bridge-Reconstruction/blob/master/read_laz.py),
[bridge.py](https://github.com/nzupancic/Bridge-Reconstruction/blob/master/bridge.py) which defines a new Bridge class
and the most important one ([bridge_reconstructor.py](https://github.com/nzupancic/Bridge-Reconstruction/blob/master/bridge_reconstructor.py)) which
processes the data and outputs a new laz file with reconstructed bridges (if there are any).

## Usage
Reconstruction can be started from command line:

`python bridge_reconstructor.py -l "lazfile.laz" -s "shapefile" -o "outputfile.laz" -i hermite`

The switches have the following roles:
- `-l` gives the path to the input laz file
- `-s` path to input shapefile (only file name but don't add the extension)
- `-o` arbitrary name for writing the output point cloud
- `-i hermite` optional switch that can be supplied to use Hermite (instead of default cosine) interpolation for 
connecting shoreline points but this also introduces a certain time penalty

Ideally the input files are located in the same directory as the script.

## Dependencies
The program was developed using Python version 3.6 and uses the following packages:
- pyshp (2.1.0)
- pylas (0.3.4)
- lazperf (1.3)
- numpy (1.17.2)

Additionally in case of issues opening laz files, following libraries and executables are included:
laszip.dll, laszip.exe, laszip64.dll, laszip64.exe.

## Results
Because of the size of the complete point cloud, the results coupled with input files are 
archived and available to download [here](https://www.dropbox.com/s/oc3if4qafiafzab/bridgeoutput.zip?dl=0). Since most of the runtime is used for 
reading the whole file I also provided point clouds (and shapefiles) with smaller area that only include one 
bridge for faster testing and execution (available in the same archive inside folders bridge1, 
bridge2, bridge3).ł