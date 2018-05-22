from numpy import array
from basegeometry import Polygon
coords = array([[1, 1], [2, 5], [4.5, 6], [8, 3],
                [7, 1], [4, 0]])
polygon = Polygon(coords)
polygon.plot()