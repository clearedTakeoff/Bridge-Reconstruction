from read_laz import LazManipulator
from read_shapefile import ShapefileReader
import math
import numpy as np


class BridgeReconstructor:

    def __init__(self, laz_file, shp_file):
        self.laz = LazManipulator(laz_file)
        shp = ShapefileReader(shp_file)
        self.bridges = shp.bridgesInsideCoords(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)

    def reconstructIt(self):
        step_size = 2
        added_points = []
        for bridge in self.bridges:
            print("Reconstructing new one")
            width = bridge.record["SIRCES"] * 100
            for i in range(len(bridge.shape.points) - 1):
                print("BLLALAAL")
                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                point_b = list(bridge.shape.points[i + 1]) + [bridge.shape.z[i + 1]]
                point_b = [x * 100 for x in point_b]
                magnitude = math.sqrt((point_a[0] - point_b[0]) * (point_a[0] - point_b[0]) +
                                      (point_a[1] - point_b[1]) * (point_a[1] - point_b[1]))
                z_step = (step_size / magnitude) * (point_b[2] - point_a[2])
                unit_vector = [(point_b[0] - point_a[0]) / magnitude, (point_b[1] - point_a[1]) / magnitude]
                for step in range(math.floor(magnitude / step_size)):
                    point_a[0] += (unit_vector[0] * step_size)
                    point_a[1] += (unit_vector[1] * step_size)
                    print("New point a is ", point_a[0], point_a[1])
                    point_a[2] += z_step
                    point_c, point_c2 = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                    for depth in range(8):
                        print("Adding point c", point_c[0], point_c[1])
                        added_points.append((point_c[0], point_c[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                        added_points.append((point_c2[0], point_c2[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
        print("Construction finished, writing...")
        self.laz.writeListToFile(added_points, "test.laz")

    def findThirdPoint(self, x0, y0, x1, y1, width):
        magnitude = math.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
        vec_x = ((x1 - x0) / magnitude) * width
        vec_y = ((y1 - y0) / magnitude) * width
        point1 = [x0 - vec_y, y0 + vec_x]
        point2 = [x0 + vec_y, y0 - vec_x]
        return point1, point2


if __name__ == "__main__":
    b = BridgeReconstructor("TM_462_101.laz", "TN_CESTE_L")
    b.reconstructIt()
