from read_laz import LazManipulator
from read_shapefile import ShapefileReader
import math
import numpy as np


class BridgeReconstructor:

    def __init__(self, laz_file, shp_file):
        self.laz = LazManipulator(laz_file)
        shp = ShapefileReader(shp_file)
        self.bridges = shp.bridgesInsideCoords(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)
        #shp.writeBridges(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)
        # f = open("bridges.txt", "r")
        # self.bridges = eval(f.read())
        # f = open("bridges.txt", "w")
        # f.write(str(self.bridges))
        # f.close()
        self.points = self.getRelevantLazPoints(self.bridges, self.laz.laz)
        # f = open("points.txt", "r")
        # self.points = eval(f.read())
        # f = open("points.txt", "w")
        # f.write(str(self.points))
        # f.close()

    def getRelevantLazPoints(self, bridges, laz):
        result = []
        for i in bridges:
            current_points = []
            for j in range(len(i.shape.points) - 1):
                point_a = list(i.shape.points[j]) + [i.shape.z[j]]
                point_a = [x * 100 for x in point_a]
                point_b = list(i.shape.points[j + 1]) + [i.shape.z[j + 1]]
                point_b = [x * 100 for x in point_b]
                magnitude = math.sqrt((point_a[0] - point_b[0]) * (point_a[0] - point_b[0]) +
                                      (point_a[1] - point_b[1]) * (point_a[1] - point_b[1]))
                unit_vector = [(point_b[0] - point_a[0]) / magnitude, (point_b[1] - point_a[1]) / magnitude]
                point_a, point_c = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1],
                                                       (i.record["SIRCES"] + 5) * 100)
                point_b = [point_a[0] + unit_vector[0] * magnitude, point_a[1] + unit_vector[1] * magnitude]
                tmp = [pnt for pnt in laz.points if self.isInsideRect(point_a, point_b, point_c, [pnt[0], pnt[1]])]
                current_points.extend(tmp)
            print("For bridge found", len(current_points), "points")
            result.append(current_points)
        return result

    # Check if point test_point is inside a rectangle defined by points p1, p2, p3
    # https://math.stackexchange.com/a/190373
    def isInsideRect(self, p1, p2, p3, test_point):
        # "create" vectors from points
        ab = [p2[0] - p1[0], p2[1] - p1[1]]
        am = [test_point[0] - p1[0], test_point[1] - p1[1]]
        ad = [p3[0] - p1[0], p3[1] - p1[1]]

        # Dot product
        ab_ab = ab[0] * ab[0] + ab[1] * ab[1]
        am_ab = am[0] * ab[0] + am[1] * ab[1]
        am_ad = am[0] * ad[0] + am[1] * ad[1]
        ad_ad = ad[0] * ad[0] + ad[1] * ad[1]
        return 0 < am_ab < ab_ab and 0 < am_ad < ad_ad


    def reconstructIt(self):
        step_size = 40
        added_points = []
        for _i in range(len(self.bridges)):
            bridge = self.bridges[_i]
            print("Reconstructing new one")
            width = bridge.record["SIRCES"] * 100
            first_point = [i * 100 for i in bridge.shape.points[0]]
            last_point = [i * 100 for i in bridge.shape.points[-1]]
            length = math.sqrt((first_point[0] - last_point[0]) * (first_point[0] - last_point[0]) +
                               (first_point[1] - last_point[1]) * (first_point[1] - last_point[1]))
            thickness = math.ceil(length / 25)
            # Calculate thickness of the reconstruction but limit it to interval [100, 250]
            thickness = min(max(thickness, 100), 250)
            print("Thickness", thickness)
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
                    point_c, point_c2 = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1],
                                                            width / 2)
                    # Find relevant points at current (x,y) position and check their z value
                    z_limit = None
                    z_limits = [pnt[2] for pnt in self.points[_i] if math.sqrt((pnt[0] - point_c[0]) ** 2 +
                                                                               (pnt[1] - point_c[1]) ** 2) < 100
                                and pnt[2] < point_a[2] - 50]
                    if len(z_limits) > 0:
                        z_limit = sum(z_limits) / len(z_limits)
                    z_limit2 = None
                    z_limits2 = [pnt[2] for pnt in self.points[_i] if math.sqrt((pnt[0] - point_c2[0]) ** 2 +
                                                                                (pnt[1] - point_c2[1]) ** 2) < 100
                                 and pnt[2] < point_a[2] - 50]
                    if len(z_limits2) > 0:
                        z_limit2 = sum(z_limits2) / len(z_limits2)
                    for depth in range(thickness):
                        if z_limit is not None and point_a[2] - depth > z_limit:
                            added_points.append((point_c[0], point_c[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                 109863062.84541322))
                        elif z_limit is None:
                            added_points.append((point_c[0], point_c[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                 109863062.84541322))
                        if z_limit2 is not None and point_a[2] - depth > z_limit2:
                            added_points.append((point_c2[0], point_c2[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                109863062.84541322))
                        elif z_limit2 is None:
                            added_points.append((point_c2[0], point_c2[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                 109863062.84541322))
                # Add points to the bottom face of bridge (reuse x,y coords from already existing points)
                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                depth = point_a[2] - thickness
                point_a, point_c = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                point_b = [point_a[0] + unit_vector[0] * magnitude, point_a[1] + unit_vector[1] * magnitude]
                points_under_bridge = [pnt for pnt in self.points[_i] if self.isInsideRect(point_a, point_b, point_c, [pnt[0], pnt[1]])]
                for point in points_under_bridge:
                    added_points.append((point[0], point[1], depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))


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
    b = BridgeReconstructor("TM_462_101_short.laz", "CESTE_SHORT")
    b.reconstructIt()
