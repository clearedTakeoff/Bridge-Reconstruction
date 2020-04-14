from read_laz import LazManipulator
from read_shapefile import ShapefileReader

class BridgeReconstructor:

    def __init__(self, laz_file, shp_file):
        laz = LazManipulator(laz_file)
        shp = ShapefileReader(shp_file)
        bridges = shp.bridgesInsideCoords(laz.x_min, laz.y_min, laz.x_max, laz.y_max)
        print("Found", len(bridges), "relevant bridges")


if __name__ == "__main__":
    BridgeReconstructor("TM_462_101.laz", "TN_CESTE_L")