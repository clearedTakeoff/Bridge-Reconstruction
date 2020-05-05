import math

class Bridge:

    def __init__(self, bridge_shp):
        point_a = list(bridge_shp.shape.points[0]) + [bridge_shp.shape.z[0]]
        point_a = [x * 100 for x in point_a]
        point_b = list(bridge_shp.shape.points[-1]) + [bridge_shp.shape.z[-1]]
        point_b = [x * 100 for x in point_b]
        magnitude = math.sqrt((point_a[0] - point_b[0]) * (point_a[0] - point_b[0]) +
                              (point_a[1] - point_b[1]) * (point_a[1] - point_b[1]))
        unit_vector = [(point_b[0] - point_a[0]) / magnitude, (point_b[1] - point_a[1]) / magnitude]
        self.point_a, self.point_c = self.findThirdPoint(point_a[0] - 150 * unit_vector[0],
                                                         point_a[1] - 150 * unit_vector[1],
                                                         point_b[0], point_b[1],
                                                         (bridge_shp.record["SIRCES"] + 5) * 100)
        point_b = [self.point_a[0] + unit_vector[0] * magnitude * 1, self.point_a[1] + unit_vector[1] * magnitude * 1]
        self.point_b = [point_b[0] + 150 * unit_vector[0], point_b[1] + 150 * unit_vector[1]]
        self.points = []

    # Returns true if point is inside the extended rectangle that belongs to the bridge
    def isRelevant(self, point):
        if self.isInsideRect(self.point_a, self.point_b, self.point_c, [point[0], point[1]]):
            self.points.append(point)
            return True
        else:
            return False

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

    def findThirdPoint(self, x0, y0, x1, y1, width):
        magnitude = math.sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1))
        vec_x = ((x1 - x0) / magnitude) * width
        vec_y = ((y1 - y0) / magnitude) * width
        point1 = [x0 - vec_y, y0 + vec_x]
        point2 = [x0 + vec_y, y0 - vec_x]
        return point1, point2
