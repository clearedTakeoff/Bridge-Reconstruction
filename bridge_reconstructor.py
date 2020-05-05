from read_laz import LazManipulator
from read_shapefile import ShapefileReader
from bridge import Bridge
import math
import numpy as np

class BridgeReconstructor:

    def __init__(self, laz_file, shp_file):
        self.laz = LazManipulator(laz_file)
        shp = ShapefileReader(shp_file)
        self.bridges_shp = shp.bridgesInsideCoords(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)
        self.bridges = [Bridge(br) for br in self.bridges_shp]
        # shp.writeBridges(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max, "CESTE_SHORT2")
        # self.points = self.getRelevantLazPoints(self.bridges_shp, self.laz.laz)
        self.points = self.getRelevantLazPoint2(self.bridges, self.laz.laz)
        for br in self.bridges:
            print("Found", len(br.points), "points")

    def getRelevantLazPoint2(self, bridges, laz):
        result = []
        for pnt in laz.points:
            for br in bridges:
                if br.isRelevant(pnt):
                    result.append(pnt)
        return result

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

                point_a, point_c = self.findThirdPoint(point_a[0] - 150 * unit_vector[0], point_a[1] - 150 * unit_vector[1], point_b[0], point_b[1],
                                                       (i.record["SIRCES"] + 5) * 100)
                point_b = [point_a[0] + unit_vector[0] * magnitude, point_a[1] + unit_vector[1] * magnitude]
                point_b = [point_b[0] + 150 * unit_vector[0], point_b[1] + 150 * unit_vector[1]]
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
        return 0 <= am_ab <= ab_ab and 0 <= am_ad <= ad_ad

    def reconstructIt(self):
        step_size = 40
        added_points = []
        for _i in range(len(self.bridges_shp)):
            bridge = self.bridges_shp[_i]
            print("Reconstructing new one")
            width = bridge.record["SIRCES"] * 100
            first_point = [i * 100 for i in bridge.shape.points[0]]
            last_point = [i * 100 for i in bridge.shape.points[-1]]
            length = math.sqrt((first_point[0] - last_point[0]) * (first_point[0] - last_point[0]) +
                               (first_point[1] - last_point[1]) * (first_point[1] - last_point[1]))
            thickness = math.ceil(length / 25)
            # Calculate thickness of the reconstruction but limit it to interval [100, 250]
            thickness = min(max(thickness, 100), 250)
            for i in range(len(bridge.shape.points) - 1):
                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                point_b = list(bridge.shape.points[i + 1]) + [bridge.shape.z[i + 1]]
                point_b = [x * 100 for x in point_b]
                magnitude = math.sqrt((point_a[0] - point_b[0]) * (point_a[0] - point_b[0]) +
                                      (point_a[1] - point_b[1]) * (point_a[1] - point_b[1]))
                z_step = (step_size / magnitude) * (point_b[2] - point_a[2])
                unit_vector = [(point_b[0] - point_a[0]) / magnitude, (point_b[1] - point_a[1]) / magnitude]

                depth = point_a[2] - thickness
                point_a, point_c = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                point_b = [point_a[0] + unit_vector[0] * magnitude * 0.5, point_a[1] + unit_vector[1] * magnitude * 0.5]
                left_side, left_z1 = self.findPointsUnderBridge(point_a, point_b, self.bridges[_i].points, depth)
                left_side_1 = [pnt for pnt in left_side if pnt[2] > left_z1 + 20]
                point_a = [point_a[0] + unit_vector[0] * magnitude * 0.5, point_a[1] + unit_vector[1] * magnitude * 0.5]
                point_b = [point_a[0] + unit_vector[0] * magnitude * 0.5, point_a[1] + unit_vector[1] * magnitude * 0.5]
                # Search for points along left side of the bridge, on one shoreline
                left_side, left_z2 = self.findPointsUnderBridge(point_a, point_b, self.bridges[_i].points, depth)
                left_side_2 = [pnt for pnt in left_side if pnt[2] > left_z2 + 20]
                point_b = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                right_side, right_z1 = self.findPointsUnderBridge(point_c, point_b, self.bridges[_i].points, depth)
                right_side_1 = [pnt for pnt in right_side if pnt[2] > right_z1 + 20]
                point_c = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                point_b = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                # Search along the other side of bridge
                right_side, right_z2 = self.findPointsUnderBridge(point_c, point_b, self.bridges[_i].points, depth)
                right_side_2 = [pnt for pnt in right_side if pnt[2] > right_z2 + 20]
                # self.generateNewBoundaryPoints(right_side_2, right_z2, depth) WHEN???

                interpolated_points = []
                # TESTING THIS DELETE later (or not???)
                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                point_b = list(bridge.shape.points[i + 1]) + [bridge.shape.z[i + 1]]
                point_b = [x * 100 for x in point_b]
                point_a, point_c = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                point_b = [point_a[0] + unit_vector[0] * magnitude, point_a[1] + unit_vector[1] * magnitude]
                points_under_bridge = [pnt for pnt in self.bridges[_i].points if
                                       self.isInsideRect(point_a, point_b, point_c, [pnt[0], pnt[1]])]
                for point in left_side_1:
                    for point2 in right_side_1:
                        dist = self.distance([point[0], point[1]], [point2[0], point2[1]])
                        if abs(point[2] - point2[2]) < 100 and 0 < dist < width + 150:
                            point1 = None
                            #for pnt in test:
                            #    if self.distance(point, pnt) < 100 and pnt != point:
                            #        point1 = pnt
                            #        break
                            point3 = None
                            #for pnt in test2:
                            #    if self.distance(point2, pnt) < 100 and pnt != point2:
                            #        point3 = pnt
                            #        break
                            interpolated_points.extend(self.interpolate(point1, point, point2, point3))
                for point in left_side_2:
                    for point2 in right_side_2:
                        dist = self.distance([point[0], point[1]], [point2[0], point2[1]])
                        if abs(point[2] - point2[2]) < 100 and 0 < dist < width + 150:
                            point1 = None
                            #for pnt in test:
                            #    if self.distance(point, pnt) < 100 and pnt != point:
                            #        point1 = pnt
                            #        break
                            point3 = None
                            #for pnt in test2:
                            #    if self.distance(point2, pnt) < 100 and pnt != point2:
                            #        point3 = pnt
                            #        break
                            interpolated_points.extend(self.interpolate(point1, point, point2, point3))

                print("Interpolated")
                for point in points_under_bridge:
                    max_z = 0
                    # z_limits = [pnt[2] for pnt in interpolated_points if math.sqrt((pnt[0] - point[0]) ** 2 +
                    #                                                               (pnt[1] - point[1]) ** 2) < 100]
                    # if len(z_limits) > 0:
                    #    max_z = max(z_limits)
                    # Calculate z value (especially visible in ascending or descending bridge)
                    depth = bridge.shape.z[i] * 100 + (self.findProjection(point_a, point_b, point) / step_size) * z_step
                    if depth > max_z:
                        # Point on the bottom face of bridge
                        added_points.append((point[0], point[1], depth - thickness, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    # Adding a point on surface of the bridge
                    added_points.append((point[0], point[1], depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    # Existing points on the bridge re-added for testing with new classification
                    # added_points.append((point[0], point[1], point[2] + 1, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    # else:
                    #   added_points.append((point[0], point[1], depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))

                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                point_b = list(bridge.shape.points[i + 1]) + [bridge.shape.z[i + 1]]
                point_b = [x * 100 for x in point_b]
                # Generate sides of the bridge
                for step in range(math.floor(magnitude / step_size)):
                    point_a[0] += (unit_vector[0] * step_size)
                    point_a[1] += (unit_vector[1] * step_size)
                    point_a[2] += z_step
                    point_c, point_c2 = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1],
                                                            width / 2)
                    # Find relevant points at current (x,y) position and check their z value
                    z_limit = None
                    #z_limits = [pnt[2] for pnt in interpolated_points if math.sqrt((pnt[0] - point_c[0]) ** 2 +
                    #                                                           (pnt[1] - point_c[1]) ** 2) < 100
                    #            and pnt[2] < point_a[2] - 90]
                    #if len(z_limits) > 0:
                    #    z_limit = max(z_limits)
                    z_limit2 = None
                    #z_limits2 = [pnt[2] for pnt in interpolated_points if math.sqrt((pnt[0] - point_c2[0]) ** 2 +
                    #                                                            (pnt[1] - point_c2[1]) ** 2) < 100
                    #            and pnt[2] < point_a[2] - 90]
                    #if len(z_limits2) > 0:
                    #    z_limit2 = max(z_limits2)
                    for depth in range(0, thickness, 15):
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
                ##interpolated_points.extend(left_side)
                #interpolated_points.extend(right_side)
                print("Interpolated", len(interpolated_points))
                for point in interpolated_points:
                    added_points.append((point[0], point[1], point[2], 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    #added_points.append((point[0], point[1], point[2] + 1, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
        print(len(added_points))
        print("Construction finished, writing...")
        self.laz.writeListToFile(added_points, "test.laz")

    def interpolate(self, point_0, point_a, point_b, point_1):
        result = []
        new_point = [point_a[0], point_a[1], point_a[2]]
        vector = [point_b[0] - point_a[0], point_b[1] - point_a[1]]
        magnitude = math.sqrt(vector[0] * vector[0] + vector[1] * vector[1])
        vector = [vector[0] / magnitude, vector[1] / magnitude]
        step = magnitude / 10
        t = 0
        while t < magnitude:
            # new_point[0] += t * vector[0]
            # new_point[1] += t * vector[1]
            if point_0 is not None and point_1 is not None:
                tmp = self.hermiteInterpolate(point_0, point_a, point_b, point_1, t / magnitude)
            else:
                tmp = self.cosInterpolate(point_a, point_b, t / magnitude)
            result.append(tuple(tmp))
            #result.append((new_point[0] + t * vector[0], new_point[1] + t * vector[1], new_point[2]))
            t += step
        return result

    def cosInterpolate(self, p1, p2, t):
        t2 = (1 - math.cos(t * math.pi)) / 2
        new_px = (p1[0] * (1 - t2) + p2[0] * t2)
        new_py = (p1[1] * (1 - t2) + p2[1] * t2)
        new_pz = (p1[2] * (1 - t2) + p2[2] * t2)
        return [new_px, new_py, new_pz]

    def hermiteInterpolate(self, p1, p2, p3, p4, t, tension=0, bias=0):
        t3 = t * t * t
        t2 = t * t
        # p3 - p1
        p2_p1 = [np.float64(p2[0] - p1[0]), np.float64(p2[1] - p1[1]), np.float64(p2[2] - p1[2])]
        p3_p2 = [np.float64(p3[0] - p2[0]), np.float64(p3[1] - p2[1]), np.float64(p3[2] - p2[2])]  # p3 - p2
        p4_p3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]]  # p4 - p3
        print(p2_p1[0] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[0] * (1 + bias) * (1 - tension) * 0.5)
        m0 = [p2_p1[0] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[0] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[1] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[1] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[2] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[2] * (1 + bias) * (1 - tension) * 0.5]
        m1 = [p2_p1[0] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[0] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[1] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[1] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[2] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[2] * (1 + bias) * (1 - tension) * 0.5]

        f1 = 2 * t3 - 3 * t2 + 1
        f2 = -(2 * t3) + 3 * t2
        f3 = t3 - 2 * t2 + t
        f4 = t3 - t2
        print(m0)
        print(m1)
        result = [f1 * p2[0] + f2 * m0[0] + f3 * m1[0] + f4 * p3[0],
                  f1 * p2[1] + f2 * m0[1] + f3 * m1[1] + f4 * p3[1],
                  f1 * p2[2] + f2 * m0[2] + f3 * m1[2] + f4 * p3[2]]
        print("Result", result)
        print(type(result[0]))
        return result


    def findThirdPoint(self, x0, y0, x1, y1, width):
        magnitude = math.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
        vec_x = ((x1 - x0) / magnitude) * width
        vec_y = ((y1 - y0) / magnitude) * width
        point1 = [x0 - vec_y, y0 + vec_x]
        point2 = [x0 + vec_y, y0 - vec_x]
        return point1, point2

    # Returns true if point p is "close" to projection on vector AB (distance less than 10)
    def isCloseToVector(self, point_a, point_b, p, dist_threshold=150):
        vector = [point_b[0] - point_a[0], point_b[1] - point_a[1]]
        p_min_a = [p[0] - point_a[0], p[1] - point_a[1]]
        project = (p_min_a[0] * vector[0] + p_min_a[1] * vector[1]) / (vector[0] * vector[0] + vector[1] * vector[1])
        if not 0 < project < 1:
            return False
        vector[0] *= project
        vector[1] *= project
        projection = [point_a[0] + vector[0], point_a[1] + vector[1]]
        if self.distance(projection, [p[0], p[1]]) < dist_threshold:
            return True
        return False

    # Finds distance from point_a to orthogonal projection of pnt onto vector a->b (2D only)
    def findProjection(self, point_a, point_b, pnt):
        vector = [point_b[0] - point_a[0], point_b[1] - point_a[1]]
        magnitude = math.sqrt(vector[0] * vector[0] + vector[1] * vector[1])
        p_min_a = [pnt[0] - point_a[0], pnt[1] - point_a[1]]
        project = (p_min_a[0] * vector[0] + p_min_a[1] * vector[1]) / (vector[0] * vector[0] + vector[1] * vector[1])
        # Project should be between 0 and 1, just in case it's not limit it by 1
        return magnitude * min(1, project)

    # Calculates distance between 2 points in 2d or 3d
    def distance(self, point_a, point_b):
        if len(point_a) > 2:
            return math.sqrt((point_a[0] - point_b[0]) ** 2 + (point_a[1] - point_b[1]) ** 2 +
                             (point_a[2] - point_b[2]) ** 2)
        else:
            return math.sqrt((point_a[0] - point_b[0]) ** 2 + (point_a[1] - point_b[1]) ** 2)

    def generateNewBoundaryPoints(self, current_list, water_level, bridge_bottom, no_of_new_points=20):
        avg_point = [0, 0, water_level + 20]
        min_z_point = current_list[0]
        print(type(avg_point[0]))
        for pnt in current_list:
            avg_point[0] += pnt[0]
            avg_point[1] += pnt[1]
            if pnt[2] < min_z_point[2]:
                min_z_point = pnt
        avg_point[0] /= len(current_list)
        avg_point[1] /= len(current_list)
        length = self.distance(avg_point, min_z_point)
        step = length / no_of_new_points
        t = 0
        while t < length:
            current_list.append(self.cosInterpolate(avg_point, min_z_point, t/length))
            t += step
        return current_list

        print(avg_point)
        while water_level < bridge_bottom:
            water_level += 20
            current_list.append([avg_point[0], avg_point[1], water_level])
        return current_list

    def findPointsUnderBridge(self, point_a, point_b, laz_points, depth):
        tmp1 = []
        z = laz_points[0][2]
        for pnt in laz_points:
            if pnt[2] < depth and self.isCloseToVector(point_a, point_b, pnt, 250):
                tmp1.append(pnt)
                z = min(z, pnt[2])
            #tmp1 = [pnt for pnt in laz_points if pnt[2] < depth and
            #        self.isCloseToVector(point_a, point_b, pnt, 250)]
        print("Found", len(tmp1), "points")
        print("Minimum is", z)
        return tmp1, z


if __name__ == "__main__":
    #b = BridgeReconstructor("TM_462_101.laz", "TN_CESTE_L")
    b = BridgeReconstructor("TM_462_101_short.laz", "CESTE_SHORT")
    #b = BridgeReconstructor("TM_462_101_short2.laz", "CESTE_SHORT2")
    b.reconstructIt()
