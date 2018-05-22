from numpy import array
from slopegeometry import NaturalSlope
surfaceCoords = array([[  0.        ,  16.57142857],
                       [ 10.        ,  16.57142857],
                       [ 18.        ,   4.57142857],
                       [ 28.        ,   4.57142857]])
slopeGeometry = NaturalSlope(surfaceCoords)
slopeGeometry.plotSlope()