from numpy import array
from slopegeometry import NaturalSlope
from circlespacking import pckCirclesInPolygon
surfaceCoords = array([[-2.4900, 18.1614],
                       [0.1022, 17.8824],
                       [1.6975, 17.2845],
                       [3.8909, 15.7301],
                       [5.8963, 14.3090],
                       [8.1183, 13.5779],
                       [9.8663, 13.0027],
                       [13.2865, 3.6058],
                       [20.2865, 3.6058],
                       [21.4347, 3.3231],
                       [22.2823, 2.7114],
                       [23.4751, 2.2252],
                       [24.6522, 1.2056],
                       [25.1701, 0.2488]])
slopeGeometry = NaturalSlope(surfaceCoords)
boundCoords = slopeGeometry.boundCoords
pckCircles = pckCirclesInPolygon(boundCoords)
pckCircles.plot(plotTriMesh=True)