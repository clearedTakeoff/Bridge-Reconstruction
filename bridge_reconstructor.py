from read_laz import LazManipulator
from read_shapefile import ShapefileReader
from bridge import Bridge
import math
import sys

class BridgeReconstructor:

    def __init__(self, laz_file, shp_file, output, interpolation=0):
        self.interpolation = interpolation  # if 0 cosine interpolation, if 1 hermite interpolation
        self.laz = LazManipulator(laz_file)
        self.output = output
        shp = ShapefileReader(shp_file)
        self.bridges_shp = shp.bridgesInsideCoords(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max)
        self.bridges = [Bridge(br) for br in self.bridges_shp]
        # shp.writeBridges(self.laz.x_min, self.laz.y_min, self.laz.x_max, self.laz.y_max, "CESTE_SHORT3")
        print("Reading points now")
        self.points = self.getRelevantLazPoints(self.bridges, self.laz.laz)
        # for br in self.bridges:
            # print("Found", len(br.points), "points")

    def getRelevantLazPoints(self, bridges, laz):
        result = []
        for pnt in laz.points:
            for br in bridges:
                if br.isRelevant(pnt):
                    result.append(pnt)
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
            print("Reconstructing new bridge")
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
                # Get points along one side of the bridge on one shoreline, returns list of points +
                # lowest z found (water surface)
                left_side, left_z1 = self.findPointsUnderBridge(point_a, point_b, self.bridges[_i].points, depth)
                left_side_1 = [pnt for pnt in left_side if pnt[2] > left_z1 + 40] # Remove water surface points
                if len(left_side_1) == 0:
                    left_side_1 = [point_b + [bridge.shape.z[i + 1]]]
                point_a = [point_a[0] + unit_vector[0] * magnitude * 0.5, point_a[1] + unit_vector[1] * magnitude * 0.5]
                point_b = [point_a[0] + unit_vector[0] * magnitude * 0.5, point_a[1] + unit_vector[1] * magnitude * 0.5]
                # Search for points along one side of the bridge, on other shoreline
                left_side, left_z2 = self.findPointsUnderBridge(point_a, point_b, self.bridges[_i].points, depth)
                left_side_2 = [pnt for pnt in left_side if pnt[2] > left_z2 + 40]
                point_b = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                # Search for points along other side of the bridge, on one shoreline
                right_side, right_z1 = self.findPointsUnderBridge(point_c, point_b, self.bridges[_i].points, depth)
                right_side_1 = [pnt for pnt in right_side if pnt[2] > right_z1 + 40]
                point_c = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                point_b = [point_c[0] + unit_vector[0] * magnitude * 0.5, point_c[1] + unit_vector[1] * magnitude * 0.5]
                # Search along the other side of bridge, on other shoreline
                right_side, right_z2 = self.findPointsUnderBridge(point_c, point_b, self.bridges[_i].points, depth)
                right_side_2 = [pnt for pnt in right_side if pnt[2] > right_z2 + 40]
                # Preprocess points (delete outliers and add new points)
                right_side_1, right_side_2, \
                left_side_1, left_side_2 = self.preprocessBoundaryPoints(right_side_1, right_z1, right_side_2,
                                                                         right_z2, left_side_1, left_z1, left_side_2,
                                                                         left_z2, depth, 0.5 * magnitude)

                interpolated_points = []
                point_a = list(bridge.shape.points[i]) + [bridge.shape.z[i]]
                point_a = [x * 100 for x in point_a]
                point_b = list(bridge.shape.points[i + 1]) + [bridge.shape.z[i + 1]]
                point_b = [x * 100 for x in point_b]
                point_a, point_c = self.findThirdPoint(point_a[0], point_a[1], point_b[0], point_b[1], width / 2)
                point_b = [point_a[0] + unit_vector[0] * magnitude, point_a[1] + unit_vector[1] * magnitude]
                points_under_bridge = [pnt for pnt in self.bridges[_i].points if
                                       self.isInsideRect(point_a, point_b, point_c, [pnt[0], pnt[1]])]
                bridge_center = [(point_a[0] + point_b[0]) * 0.5, (point_a[1] + point_b[1]) * 0.5]
                max_right_z1 = max_left_z1 = max_right_z2 = max_left_z2 = None
                # Calculate how many steps in the interpolation (more points = less steps)
                # But limit between 10 and 30, to not create too many points
                steps = 10000 / (len(left_side_1) + len(right_side_1))
                steps = max(10, min(steps, 30))
                # Only to this in the first bridge segment
                if len(left_side_1) > 0 and len(right_side_1) > 0 and i == 0:
                    # Finds the highest part of each shoreline, used for vector along the shoreline
                    # to fit bottom face of bridge to terrain
                    right_side_done = 0
                    max_left_z1 = left_side_1[0]
                    max_right_z1 = right_side_1[0]

                    for point in left_side_1:
                        if point[2] > max_left_z1[2]:
                            max_left_z1 = point
                        for point2 in right_side_1:
                            if right_side_done == 0:
                                if point2[2] > max_right_z1[2]:
                                    max_right_z1 = point2
                            dist = self.distance([point[0], point[1]], [point2[0], point2[1]])
                            if abs(point[2] - point2[2]) < 100 and 0 < dist < width + 200:
                                point1 = None
                                point3 = None
                                if self.interpolation == 1:
                                    dst1 = 999999999
                                    for pnt in left_side_1:
                                        dst = self.distance(point, pnt)
                                        if dst1 > dst > 0 and pnt != point:
                                            point1 = pnt
                                            dst1 = dst

                                    dst3 = 999999999
                                    for pnt in right_side_1:
                                        dst = self.distance(point2, pnt)
                                        if 0 < dst < dst3 and pnt != point2:
                                            point3 = pnt
                                            dst3 = dst
                                    if self.distance(point, point2) > self.distance(point1, point2):
                                        tmp = point
                                        point = point1
                                        point1 = tmp
                                interpolated_points.extend(self.interpolate(point1, point, point2, point3, steps))
                        right_side_done = 1
                # Only do this in the last section of the bridge
                if len(left_side_2) > 0 and len(right_side_2) > 0 and i == len(bridge.shape.points) - 2:
                    max_left_z2 = left_side_2[0]
                    max_right_z2 = right_side_2[0]
                    right_side_done = 0
                    steps = 10000 / (len(left_side_2) + len(right_side_2))
                    steps = max(10, min(steps, 30))
                    for point in left_side_2:
                        if point[2] > max_left_z2[2]:
                            max_left_z2 = point
                        for point2 in right_side_2:
                            if right_side_done == 0 and point2[2] > max_right_z2[2]:
                                max_right_z2 = point2
                            dist = self.distance([point[0], point[1]], [point2[0], point2[1]])
                            if abs(point[2] - point2[2]) < 100 and 0 < dist < width + 200:
                                point1 = None
                                point3 = None
                                if self.interpolation == 1:
                                    dst1 = 999999999
                                    for pnt in left_side_2:
                                        dst = self.distance(point, pnt)
                                        if dst1 > dst > 0 and pnt != point:
                                            point1 = pnt
                                            dst1 = dst

                                    dst3 = 999999999
                                    for pnt in right_side_2:
                                        dst = self.distance(point2, pnt)
                                        if 0 < dst < dst3 and pnt != point2:
                                            point3 = pnt
                                            dst3 = dst
                                    if self.distance(point, point2) > self.distance(point1, point2):
                                        tmp = point
                                        point = point1
                                        point1 = tmp
                                interpolated_points.extend(self.interpolate(point1, point, point2, point3, steps))
                        right_side_done = 1

                # Check on which side of the vector bridge center point lies and use that as a reference
                side1_center = self.whichSide(max_right_z1, max_left_z1, bridge_center)
                side2_center = self.whichSide(max_right_z2, max_left_z2, bridge_center)
                for point in points_under_bridge:
                    # Calculate z value (especially visible in ascending or descending bridge)
                    depth = bridge.shape.z[i] * 100 + (self.findProjection(point_a, point_b, point) / step_size) * z_step
                    testTerrain1 = self.whichSide(max_right_z1, max_left_z1, point)
                    testTerrain2 = self.whichSide(max_right_z2, max_left_z2, point)
                    # If point lies on the correct side of vector add it to the bottom face
                    if testTerrain1 == side1_center and testTerrain2 == side2_center:
                        # Point on the bottom face of bridge
                        added_points.append((point[0], point[1], depth - thickness, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    # Adding a point on surface of the bridge
                    added_points.append((point[0], point[1], depth, 20, 82, 65, 32, 0, 2418, 109863062.84541322))

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
                    for depth in range(0, thickness, 15):
                        testTerrain1 = self.whichSide(max_right_z1, max_left_z1, point_c)
                        testTerrain2 = self.whichSide(max_right_z2, max_left_z2, point_c)
                        if testTerrain1 == side1_center and testTerrain2 == side2_center:
                            added_points.append((point_c[0], point_c[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                109863062.84541322))

                        testTerrain1 = self.whichSide(max_right_z1, max_left_z1, point_c2)
                        testTerrain2 = self.whichSide(max_right_z2, max_left_z2, point_c2)
                        if testTerrain1 == side1_center and testTerrain2 == side2_center:
                            added_points.append((point_c2[0], point_c2[1], point_a[2] - depth, 20, 82, 65, 32, 0, 2418,
                                                 109863062.84541322))
                for point in interpolated_points:
                    added_points.append((point[0], point[1], point[2], 20, 82, 65, 32, 0, 2418, 109863062.84541322))
                    #added_points.append((point[0], point[1], point[2] + 1, 20, 82, 65, 32, 0, 2418, 109863062.84541322))
        print("Construction finished, writing...")
        self.laz.writeListToFile(added_points, self.output)

    # Function that generates new points between point_a and point_b. If point_0 and point_1 are also given it uses
    # Hermite interpolation, else only cosine interpolation
    def interpolate(self, point_0, point_a, point_b, point_1, steps=15):
        result = []
        new_point = [point_a[0], point_a[1], point_a[2]]
        vector = [point_b[0] - point_a[0], point_b[1] - point_a[1]]
        magnitude = math.sqrt(vector[0] * vector[0] + vector[1] * vector[1])
        vector = [vector[0] / magnitude, vector[1] / magnitude]
        step = magnitude / steps
        t = 0
        while t < magnitude:
            # new_point[0] += t * vector[0]
            # new_point[1] += t * vector[1]
            if point_0 is not None or point_1 is not None:
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
        p2_p1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]
        p3_p2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]]  # p3 - p2
        p4_p3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]]  # p4 - p3

        m0 = [p2_p1[0] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[0] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[1] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[1] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[2] * (1 + bias) * (1 - tension) * 0.5 + p3_p2[2] * (1 + bias) * (1 - tension) * 0.5]
        m1 = [p2_p1[0] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[0] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[1] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[1] * (1 + bias) * (1 - tension) * 0.5,
              p2_p1[2] * (1 + bias) * (1 - tension) * 0.5 + p4_p3[2] * (1 + bias) * (1 - tension) * 0.5]

        f1 = 2 * t3 - 3 * t2 + 1
        f2 = t3 - 2 * t2 + t
        f3 = t3 - t2
        f4 = -2 * t3 + 3 * t2

        result = [f1 * p2[0] + f2 * m0[0] + f3 * m1[0] + f4 * p3[0],
                  f1 * p2[1] + f2 * m0[1] + f3 * m1[1] + f4 * p3[1],
                  f1 * p2[2] + f2 * m0[2] + f3 * m1[2] + f4 * p3[2]]
        # print("Result", result)
        return result

    # Function that removes outliers and/or generates additional points if any are missing
    # Args are 4 lists with points for each side and each shore, accompanied by the level of water on each side
    # and depth of the bridge (z value for bottom of the bridge)
    def preprocessBoundaryPoints(self, right_side_1, right_z1, right_side_2, right_z2,
                                 left_side_1, left_z1, left_side_2, left_z2, depth, magnitude,
                                 generate_points=True, remove_outliers=True):
        density = 0.15
        if remove_outliers is True:
            right_side_1 = self.removeOutliers(right_side_1)
            right_side_2 = self.removeOutliers(right_side_2)
            left_side_1 = self.removeOutliers(left_side_1)
            left_side_2 = self.removeOutliers(left_side_2)
        if generate_points is True:
            no_of_points = math.ceil(magnitude * density)
            right_side_1 = self.generateNewBoundaryPoints(right_side_1, right_z1, depth, no_of_points - len(right_side_1))
            right_side_2 = self.generateNewBoundaryPoints(right_side_2, right_z2, depth, no_of_points - len(right_side_2))
            left_side_1 = self.generateNewBoundaryPoints(left_side_1, left_z1, depth, no_of_points - len(left_side_1))
            left_side_2 = self.generateNewBoundaryPoints(left_side_2, left_z2, depth, no_of_points - len(left_side_2))
        return right_side_1, right_side_2, left_side_1, left_side_2


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

    # Checks whether test point lies on left side, right side or on the vector defined by point_a->point_b
    # Returns 1 if on the right, -1 if left or 0 if on vector
    # https://stackoverflow.com/a/1560510
    def whichSide(self, point_a, point_b, test_point):
        if point_a is None and point_b is None:
            return 1
        test = (point_b[0] - point_a[0]) * (test_point[1] - point_a[1]) - \
               (point_b[1] - point_a[1]) * (test_point[0] - point_a[0])
        if test > 0:
            return 1
        if test < 0:
            return -1
        return 0

    # Calculates distance between 2 points in 2d or 3d
    def distance(self, point_a, point_b):
        if len(point_a) > 2:
            return math.sqrt((point_a[0] - point_b[0]) ** 2 + (point_a[1] - point_b[1]) ** 2 +
                             (point_a[2] - point_b[2]) ** 2)
        else:
            return math.sqrt((point_a[0] - point_b[0]) ** 2 + (point_a[1] - point_b[1]) ** 2)

    def generateNewBoundaryPoints(self, current_list, water_level, bridge_bottom, no_of_new_points=20):
        if len(current_list) == 0:
            return current_list
        min_z_point = current_list[0]
        max_z_point = current_list[0]
        # Calculate complete distance traveled from first to last point in list
        # Used to calculate how many new points to insert based on distance and no. of new points needed
        complete_distance = 0
        for i in range(1, len(current_list)):
            pnt = current_list[i]
            complete_distance += self.distance(pnt, current_list[i - 1])
            if pnt[2] < min_z_point[2]:
                min_z_point = pnt
            if pnt[2] > max_z_point[2]:
                max_z_point = pnt
        try:
            points_per_unit = no_of_new_points / complete_distance
        except ZeroDivisionError:
            points_per_unit = 1
        max_z_point[2] = bridge_bottom - 10
        min_z_point[2] = water_level
        current_list.insert(0, max_z_point)
        current_list.append(min_z_point)
        min_z_point[2] -= 10
        current_list.append(min_z_point)
        min_z_point[2] -= 10
        current_list.append(min_z_point)
        min_z_point[2] -= 10
        current_list.append(min_z_point)
        interpolated = []
        for i in range(1, len(current_list) - 2):
            p1 = current_list[i - 1]
            p2 = current_list[i]
            p3 = current_list[i + 1]
            p4 = current_list[i + 2]
            length = self.distance(p2, p3)
            no_of_points = length * points_per_unit
            if no_of_points >= 1:
                step = length / no_of_points
                t = 0
                while t < length:
                    interpolated.append(self.hermiteInterpolate(p1, p2, p3, p4, t / length))
                    t += step
        current_list.extend(interpolated)
        return current_list

    def removeOutliers(self, cloud):
        if len(cloud) == 0:
            return cloud
        to_remove = []
        dist_to_neighbors = []
        global_avg = 0
        for pnt in cloud:
            top5 = 10 * [(9999999, None)]
            for pnt2 in cloud:
                if pnt2 != pnt and pnt2:
                    dist = self.distance(pnt, pnt2)
                    if dist < top5[-1][0]:
                        i = 0
                        while dist > top5[i][0]:
                            i += 1
                        top5.insert(i, (dist, pnt2))
                        top5 = top5[:-1]
            avg_dist = 0
            for entry in top5:
                avg_dist += entry[0]
            avg_dist /= len(top5)
            dist_to_neighbors.append((avg_dist, pnt))
            global_avg += avg_dist
            #for entry in top5:
             #   if entry[0] > 2 * avg_dist:
              #      to_remove.append(entry[1])
        global_avg /= len(cloud)
        for pnt in dist_to_neighbors:
            if pnt[0] > 1.7 * global_avg:
                to_remove.append(pnt[1])
        return [point for point in cloud if point not in to_remove]

    def findPointsUnderBridge(self, point_a, point_b, laz_points, depth):
        tmp1 = []
        z = laz_points[0][2]
        for pnt in laz_points:
            if pnt[2] < depth and self.isCloseToVector(point_a, point_b, pnt, 250):
                tmp1.append(pnt)
                z = min(z, pnt[2])
        return tmp1, z


if __name__ == "__main__":
    # Default  values
    interpolation = 0
    output = "test.laz"
    input_laz = "TM_462_101.laz"
    input_shp = "TN_CESTE_L"
    #input_laz = "TM_462_101_short.laz"
    #input_shp = "CESTE_SHORT"
    #input_laz = "TM_462_101_short2.laz"
    #input_shp = "CESTE_SHORT2"
    #input_laz = "TM_462_101_short3.laz"
    #input_shp = "CESTE_SHORT3"
    # Parse command line arguments
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv) - 1, 2):
            opt = sys.argv[i]
            arg = sys.argv[i + 1]
            if opt == "-l":
                input_laz = arg
            elif opt == "-s":
                input_shp = arg
            elif opt == "-i":
                if arg == "hermite":
                    interpolation = 1
                else:
                    interpolation = 0
            elif opt == "-o":
                output = arg
    b = BridgeReconstructor(input_laz, input_shp, output, interpolation)
    b.reconstructIt()
