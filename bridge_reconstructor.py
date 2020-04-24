from read_laz import LazManipulator
from read_shapefile import ShapefileReader
import math
import numpy as np


class BridgeReconstructor:

    def __init__(self, laz_file, shp_file):
        self.laz = LazManipulator(laz_file)
        shp = ShapefileReader(shp_file)
        self.bridges = shp.bridgesInsideCoords(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)
        self.points = self.getRelevantLazPoints(self.bridges, self.laz.laz  )

    def getRelevantLazPoints(self, bridges, laz):
        result = []
        for i in bridges:
            x_min = i.shape.bbox[0]
            y_min = i.shape.bbox[1]
            x_max = i.shape.bbox[2]
            y_max = i.shape.bbox[3]
            current_points = laz.points[np.logical_and(np.logical_and(laz.x >= x_min, laz.x <= x_max),
                                                       np.logical_and(laz.y >= y_min, laz.y <= y_max))]
            print("For bridge found", len(current_points), "points")
            result.append(current_points)
        return result

    def reconstructIt(self):
        step_size = 20
        added_points = []
        for _i in range(len(self.bridges)):
            bridge = self.bridges[_i]
            print("Reconstructing new one")
            width = bridge.record["SIRCES"] * 100
            first_point = [i * 100 for i in bridge.shape.points[0]]
            last_point = [i * 100 for i in bridge.shape.points[-1]]
            length = math.sqrt((first_point[0] - last_point[0]) * (first_point[0] - last_point[0]) +
                               (first_point[1] - last_point[1]) * (first_point[1] - last_point[1]))
            thickness = math.ceil(length / 8)
            for i in range(len(bridge.shape.points) - 1):
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
                    point_a[2] += z_step
                    point_c, point_c2 = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                    # Find relevant points at current (x,y) position and check their z value
                    # z_limits = [pnt[2] for pnt in self.points if pnt[0] == point_c[0] and pnt[1] == point_c[1]]
                    for depth in range(thickness):
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
    #b.reconstructIt()
