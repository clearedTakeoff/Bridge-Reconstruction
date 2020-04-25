import shapefile


class ShapefileReader:

    def __init__(self, filename):
        # Open and read the shapefile (file extension not required, only filename)
        self.sf = shapefile.Reader(filename)
        # Save all entries in the file, each entry containing both shape and record object
        self.entries = self.sf.shapeRecords()
        self.bridges = []
        for entry in self.entries:
            # Only extract shapeRecord entries of bridges (tipobj 3), limit also by coordinates??
            if entry.record["TIPOBJ_CES"] == 3:
                self.bridges.append(entry)
        print(len(self.bridges))

    # Returns bridges inside coordinates bound between points (lowX, lowY) and (highX, highY)
    def bridgesInsideCoords(self, lowX, lowY, highX, highY):
        bridges = []
        for bridge in self.bridges:
            if lowX <= bridge.shape.bbox[0] <= highX and lowY <= bridge.shape.bbox[1] <= highY\
                    and lowX <= bridge.shape.bbox[2] <= highX and lowY <= bridge.shape.bbox[3] <= highY:
                bridges.append(bridge)
        print("Found", len(bridges), "bridges")
        return bridges

    def writeBridges(self, lowX, lowY, highX, highY):
        sf = shapefile.Writer("CESTE_SHORT")
        sf.fields = self.sf.fields[1:]
        for bridge in self.sf.iterShapeRecords():
            if lowX <= bridge.shape.bbox[0] <= highX and lowY <= bridge.shape.bbox[1] <= highY\
                    and lowX <= bridge.shape.bbox[2] <= highX and lowY <= bridge.shape.bbox[3] <= highY:
                sf.record(*bridge.record)
                sf.shape(bridge.shape)
        sf.close()

if __name__ == "__main__":
    s = ShapefileReader("TN_CESTE_L")
