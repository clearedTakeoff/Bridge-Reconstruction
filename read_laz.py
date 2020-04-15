import pylas
import numpy as np
# Lazpef >= 1.3 required for writing .laz file

class LazManipulator:

    def __init__(self, filename):
        fh = pylas.open(filename)
        self.laz = fh.read()
        self.x_min = self.laz.x.min()
        self.x_max = self.laz.x.max()
        self.y_min = self.laz.y.min()
        self.y_max = self.laz.y.max()
        #self.points = list(zip(self.laz.x, self.laz.y))
        self.points = self.laz.points
        self.dt = self.points.dtype
        # print(self.points)
        # print("X min ", self.x_min, "X max", self.x_max, "Y min", self.y_min, "Y max", self.y_max)

    # Adds new points to the point cloud and writes it to an output file
    # Expects new_points to be a nested list where each entry/line consists of 10 parameters:
    # [x, y, z, intensity, bit_fiels, classification, scan_angle_rank, user_data, point_source_id, gps_time]
    def writeListToFile(self, new_points, name="output.laz"):
        #for entry in new_points:
            #new_point = np.array(tuple(entry), dtype=self.dt)
        #    self.points = np.append(self.points, entry)
        #    print("New point")
        numpy_array = np.array(new_points, dtype=self.dt)
        self.points = np.append(self.points, numpy_array)
        print("End for")
        self.laz.points = self.points
        self.laz.write(name)

if __name__ == "__main__":
    LazManipulator("GK_430_136.laz")