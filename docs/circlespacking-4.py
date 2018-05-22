from numpy import array
from basegeometry import Polygon
from circlespacking import pckCirclesInPolygon
coordinates = array([[1, 1], [2, 5], [4.5, 6], [6, 4], [8, 3],
                     [7, 1], [4.5, 1], [4, 0]])
polygon = Polygon(coordinates)
boundCoords = polygon.boundCoords
pckCircles = pckCirclesInPolygon(boundCoords, 10)
pckCircles.loglogDiagram()