import pylas


class LazReader:

    def __init__(self, filename):
        fh = pylas.open(filename)
        laz = fh.read()
        self.x_min = laz.x.min()
        self.x_max = laz.x.max()
        self.y_min = laz.y.min()
        self.y_max = laz.y.max()
        self.points = list(zip(laz.x, laz.y))
        # print(self.points)
        # print("X min ", self.x_min, "X max", self.x_max, "Y min", self.y_min, "Y max", self.y_max)


if __name__ == "__main__":
    LazReader("GK_430_136.laz")