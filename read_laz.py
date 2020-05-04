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
        self.laz.points = self.points
        self.laz.write(name)

    def filterPoints(self, min_x, min_y, max_x, max_y):
        result = []
        for pnt in self.points:
            if min_x <= pnt[0] <= max_x and min_y <= pnt[1] <= max_y:
                result.append(pnt)
        return result

if __name__ == "__main__":
    l = LazManipulator("TM_462_101.laz")
    tmp = l.filterPoints(46200403, 10107873, 46206213, 10115060)
    l.points = tmp
    l.writeListToFile(tmp, "TM_462_101_short2.laz")
