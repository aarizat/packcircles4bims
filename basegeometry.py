# -*- coding: utf-8 -*-
'''Module to define the ``Triangle``, ``Circle`` and ``Polygon`` classes. Some
properties related to the geometry of these classes are determined.
These classes are the basic inputs to pack circular particles in a
closed polygon in :math:`\\mathbb{R}^2`.
'''


# %%
class Circle:
    '''Creates an instance of an object that defines a Circle once the
    cartesian coordinates of its center and the radius are given.

    Attributes:
        center (`tuple` or `list`): (x, y)-cartesian coordinates of\
            circle center.
        radius (`float` or `int`): Length of the segment that joins the center\
            with any point of the circumference.

    Examples:
        >>> center, radius = (0, 0), 1
        >>> circle = Circle(center, radius)
        >>> circle.__dict__
        {'area': 3.141592653589793,
         'center': array([0, 0]),
         'curvature': 1.0,
         'diameter': 2,
         'perimeter': 6.283185307179586,
         'radius': 1}

        >>> center, radius = (2, 5), 2.5
        >>> circle = Circle(center, radius)
        >>> circle.__dict__
        {'area': 19.634954084936208,
        'center': array([2, 5]),
        'curvature': 0.4,
        'diameter': 5.0,
        'perimeter': 15.707963267948966,
        'radius': 2.5}

    Note:
        The class ``Circle`` requires `NumPy <http://www.numpy.org/>`_
    '''

    def __init__(self, center, radius):
        '''Method for initializing the attributes of the class.'''
        import numpy as np
        self.center = np.array(center)
        self.radius = radius
        self.curvature = 1 / radius
        self.diameter = 2 * radius
        self.area = np.pi * radius**2
        self.perimeter = 2 * radius * np.pi

    def descartesTheorem(self, circle1, circle2=None):
        '''
        Method to determine the tangent circles of the `Descartes theorem\
        <https://en.wikipedia.org/wiki/Descartes%27_theorem#Special_cases>`_.

        To find centers of these circles, it calculates the
        intersection points of two circles by using the construction of
        triangles, proposed by `Paul Bourke, 1997\
        <http://paulbourke.net/geometry/circlesphere/>`_.

        Parameters:
            circle1 (`circle` object): Tangent circle to the circle object\
                intantiated.
            circle2 (`circle` object): Tangent circle to the circle object\
                intantiated and to the `circle1`.

       Returns:
           circles (`tuple`): Each element of the tuple is a circle object.

       Examples:
           >>> import matplotlib.pyplot as plt
           >>> from basegeometry import Circle
           >>> # Special case Descartes' Theorem (cicle with infite radius)
           >>> circle = Circle((4.405957, 2.67671461), 0.8692056336001268)
           >>> circle1 = Circle((3.22694724, 2.10008003), 0.4432620600509628)
           >>> c2, c3 = circle.descartesTheorem(circle1)
           >>> # plotting
           >>> plt.axes()
           >>> plt.gca().add_patch(plt.Circle(circle.center,
                                   circle.radius, fill=False))
           >>> plt.gca().add_patch(plt.Circle(circle1.center,
                                   circle1.radius, fill=False))
           >>> plt.gca().add_patch(plt.Circle(c2.center,
                                   c2.radius, fc='r'))
           >>> plt.gca().add_patch(plt.Circle(c3.center,
                                   c3.radius, fc='r'))
           >>> plt.axis('equal')
           >>> plt.show()

           >>> import matplotlib.pyplot as plt
           >>> from basegeometry import Circle
           >>> # General case Descartes Theorem (three circle tangent mutually)
           >>> circle = Circle((4.405957, 2.67671461), 0.8692056336001268)
           >>> circle1 = Circle((3.22694724, 2.10008003), 0.4432620600509628)
           >>> circle2 = Circle((3.77641134, 1.87408749), 0.1508620255299397)
           >>> c3, c4 = circle.descartesTheorem(circle1, circle2)
           >>> # plotting
           >>> plt.axes()
           >>> plt.gca().add_patch(plt.Circle(circle.center,
                                   circle.radius, fill=False))
           >>> plt.gca().add_patch(plt.Circle(circle1.center,
                                   circle1.radius, fill=False))
           >>> plt.gca().add_patch(plt.Circle(circle2.center,
                                   circle2.radius, fill=False))
           >>> plt.gca().add_patch(plt.Circle(c3.center,
                                   c3.radius, fc='r'))
           >>> plt.axis('equal')
           >>> plt.show()
        '''

        if circle2 is None:
            # Special case Descartes' theorem
            radius = (self.curvature + circle1.curvature +
                      2*(self.curvature*circle1.curvature)**0.5)**-1
        else:
            # General case Descartes Theorem
            radius = (self.curvature + circle1.curvature +
                      circle2.radius + circle2.curvature +
                      2*((self.curvature*circle1.curvature) +
                         (circle1.curvature*circle2.curvature) +
                         (circle2.curvature*self.curvature))**0.5)**-1
        # Distances between centers and intersection points
        R1, R2 = self.radius + radius, circle1.radius + radius
        # Distance between the centers of the intersected circles
        dist = self.radius + circle1.radius
        cos, sin = (circle1.center - self.center) / dist
        # Distance to the chord
        chordDist = (R1**2 - R2**2 + dist**2) / (2*dist)
        # Half-length of the chord
        halfChord = (R1**2 - chordDist**2)**0.5
        center3 = (self.center[0] + chordDist*cos - halfChord*sin,
                   self.center[1] + chordDist*sin + halfChord*cos)
        center4 = (self.center[0] + chordDist*cos + halfChord*sin,
                   self.center[1] + chordDist*sin - halfChord*cos)
        circles = Circle(center3, radius), Circle(center4, radius)
        return circles


# %%
class Triangle:
    '''Creates an instance of an object that defines a Triangle once the
    coordinates of its three vertices in cartesian :math:`\\mathbb{R}^2` space
    are given.

    It considers the usual notation for the triangle ``ABC`` in which `A`, `B`
    and `C` represent the vertices and ``a``, ``b``, ``c`` are the
    lengths of `the segments BC`, `CA` and `AB` respectively.

    Attributes:
        coordinates ((3, 2) `numpy.ndarray`): Coordinates of three vertices\
            of the triangle.

    Note:
        The class ``Triangle`` requires `NumPy <http://www.numpy.org/>`_,
        `SciPy <https://www.scipy.org/>`_\
        and `Matplotlib <https://matplotlib.org/>`_.

    Examples:
        >>> from numpy import array
        >>> from basegeometry import Triangle
        >>> coords = array([(2, 1.5), (4.5, 4), (6, 2)])
        >>> triangle = Triangle(coords)
        >>> triangle.__dict__.keys()
        dict_keys(['vertices', 'area', 'sides', 'perimeter', 'distToIncenter',
                   'incircle'])

        >>> from numpy import array
        >>> from basegeometry import Triangle
        >>> coords = array([[2, 1], [6, 1], [4, 5.5]])
        >>> triangle = Triangle(coords)
        >>> triangle.__dict__.keys()
        dict_keys(['vertices', 'area', 'sides', 'perimeter', 'distToIncenter',
                   'incircle'])
    '''

    def __init__(self, coordinates):
        '''Method for initializing the attributes of the class.'''
        self.vertices = dict(zip('ABC', coordinates))
        # same geometric properties of the triangle
        self.getGeomProperties()

    def getGeomProperties(self):
        '''Method to set the attributes to the instanced object

        The established geometric attributes are the following:
            * Area.
            * Lenght of its three sides.
            * Perimeter.
            * Incircle (``Circle`` Object)
            * Distance of each vertice to incenter.
        '''

        from scipy.spatial.distance import euclidean
        import numpy as np

        vertsArray = np.array([*self.vertices.values()])
        v = np.vstack((vertsArray, vertsArray[0]))
        # Area (Gauss Equation)
        area = 0.5*abs(sum(v[:-1, 0] * v[1:, 1] - v[:-1, 1] * v[1:, 0]))
        # length three sides triangle.
        sides = {'a': euclidean(self.vertices['B'], self.vertices['C']),
                 'b': euclidean(self.vertices['A'], self.vertices['C']),
                 'c': euclidean(self.vertices['A'], self.vertices['B'])}
        perimeter = sides['a'] + sides['b'] + sides['c']
        # Inscribed circle radius
        radius = 2*area / perimeter
        # Inscribed circle center
        center = (sides['a'] * self.vertices['A'] +
                  sides['b'] * self.vertices['B'] +
                  sides['c'] * self.vertices['C']) / perimeter
        # Distances of each vertice to incenter
        distToIncenter = [euclidean(center, v) for v in self.vertices.values()]
        # Set the attribute to the instanced object.
        setattr(self, 'area', area)
        setattr(self, 'sides', sides)
        setattr(self, 'perimeter', perimeter)
        setattr(self, 'distToIncenter', distToIncenter)
        setattr(self, 'incircle', Circle(center, radius))
        return

    def packCircles(self, depth=None, want2plot=False):
        '''Method to pack circular particles within of a triangle. It apply
        the Descartes theorem (special and general case) to generate mutually
        tangent circles in a fractal way in the triangle.

        Parameters:
            depth (`int`): Fractal depth. Number that indicate how many\
                circles are fractally generated from the `incirle` to each\
                vertice of the triangle. If this number is not given, then,\
                the fractal generation of circles is done up to a circle\
                reachs a radius to lower than the five percent of the\
                incircle radius.
            want2plot (`bool`): Variable to check if a plot is wanted.\
                The default value is ``False``.

        Returns:
            listCircles (`list`): `list` that contains all the circular\
                particles packed in the triangle.

        Note:
            Large values of `depth` might produce internal variables that tend
            to infinte, then a ``ValueError`` is produced with a warning
            message ``array must not contain infs or NaNs``.

        Examples:
            >>> from numpy import array
            >>> from basegeometry import Triangle
            >>> coords = array([(2, 1.5), (4.5, 4), (6, 2)])
            >>> triangle = Triangle(coords)
            >>> cirsInTri = triangle.packCircles(depth=2, want2plot=True)

            >>> from numpy import array
            >>> from basegeometry import Triangle
            >>> coords = array([[2, 1], [6, 1], [4, 5.5]])
            >>> triangle = Triangle(coords)
            >>> cirsInTri = triangle.packCircles(depth=5, want2plot=True)

            .. plot::

                from numpy import array
                from basegeometry import Triangle
                coords = array([(2, 1), (2, 8), (7, 1)])
                triangle = Triangle(coords)
                triangle.packCircles(want2plot=True)
        '''

        import numpy as np
        from scipy.spatial.distance import euclidean
        import matplotlib.pyplot as plt

        listCircles = [self.incircle]
        for vert, distance in zip(self.vertices.values(), self.distToIncenter):
            auxCircle = self.incircle
            auxDist = distance
            if not depth:
                while True:
                    radius = (auxCircle.radius*auxDist -
                              auxCircle.radius**2) / (auxCircle.radius+auxDist)
                    # Checking condition stop
                    # (circle with center on the bisectrix)
                    if radius < 0.050*self.incircle.radius:
                        break
                    dist = auxCircle.radius + radius  # dist between centers
                    center = np.array(auxCircle.center) + (dist/auxDist) * (
                            np.array(vert) - np.array(auxCircle.center))
                    circle = Circle(center, radius)
                    # Genereting circles within triangle (Descartes circles)
                    c31, c32 = auxCircle.descartesTheorem(circle)
                    c41, c42 = auxCircle.descartesTheorem(circle, c31)
                    c51, c52 = circle.descartesTheorem(c31)
                    c61, c62 = circle.descartesTheorem(c32)
                    c71, c72 = auxCircle.descartesTheorem(c31)
                    c81, c82 = auxCircle.descartesTheorem(c32)
                    listCircles.extend((circle, c31, c32, c41, c42, c52, c61,
                                        c71, c82))
                    # Updating variables
                    auxDist = euclidean(vert, circle.center)
                    auxCircle = circle
            else:
                for i in range(1, depth+1):
                    radius = (auxCircle.radius*auxDist -
                              auxCircle.radius**2) / (auxCircle.radius+auxDist)
                    dist = auxCircle.radius + radius  # dist between centers
                    center = np.array(auxCircle.center) + (dist/auxDist) * (
                            np.array(vert) - np.array(auxCircle.center))
                    circle = Circle(center, radius)
                    # Genereting circles within triangle (Descartes circles)
                    c31, c32 = auxCircle.descartesTheorem(circle)
                    c41, c42 = auxCircle.descartesTheorem(circle, c31)
                    c51, c52 = circle.descartesTheorem(c31)
                    c61, c62 = circle.descartesTheorem(c32)
                    c71, c72 = auxCircle.descartesTheorem(c31)
                    c81, c82 = auxCircle.descartesTheorem(c32)
                    listCircles.extend((circle, c31, c32, c41, c42, c52, c61,
                                        c71, c82))
                    # Updating variables
                    auxDist = euclidean(vert, circle.center)
                    auxCircle = circle
        # Set the attribute to the instanced object.
        setattr(self, 'listCircles', listCircles)
        # plotting
        if want2plot:
            figure = self.plotTriangle()
            for circle in listCircles:
                figure.add_patch(plt.Circle(circle.center, circle.radius,
                                            fill=False, lw=1, ec='k'))
        else:
            return listCircles

    def plotTriangle(self):
        '''Method for show a graphic of the triangle object.

        Returns:
            ax (`matplotlib.axes._subplots.AxesSubplot`): object associated\
                with a matplotlib graphic.

        Examples:
            >>> from numpy import array
            >>> from basegeometry import Triangle
            >>> coords = array([(1, 1), (4, 8), (8, 5)])
            >>> triangle = Triangle(coords)
            >>> triangle.plotTriangle()

            .. plot::

                from numpy import array
                from basegeometry import Triangle
                coords = array([(2, 1), (2, 8), (7, 1)])
                triangle = Triangle(coords)
                triangle.plotTriangle()
        '''
        import matplotlib.pyplot as plt
        import numpy as np

        vert = np.array([*self.vertices.values()])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # plotting
        ax.plot(np.hstack((vert[:, 0], vert[0, 0])),
                np.hstack((vert[:, 1], vert[0, 1])),
                'k-', lw=2, label='Triangle')
        ax.axis('equal')
        ax.grid(ls='--', lw=0.6)
        return ax


# %%
class Polygon:
    '''Creates an instance of an object that defines a Polygon once the
    cartesian coordinates of its vertices are given.

    Attributes:
        coordinates ((n, 2) `numpy.ndarray`): Coordinates of the vertices\
            of the polygon.
    '''

    def __init__(self, coordinates):
        '''Method for initializing the attributes of the class.'''
        import numpy as np
        self.coordinates = coordinates
        self.boundCoords = np.vstack((coordinates, coordinates[0]))
        self.area()

    def area(self):
        '''Method for determine the area of the polygon.

        Returns:
            area (`float`): area of the polygon surface.

        Examples:
            >>> from numpy import array
            >>> from basegeometry import Polygon
            >>> coords = array([(1, 1), (4, 8), (8, 5)])
            >>> polygon = Polygon(coords)
            >>> polygon.area
            18.5

            >>> from numpy import array
            >>> from basegeometry import Polygon
            >>> coords = array([[1, 1], [2, 5], [4.5, 6], [8, 3],
                                [7, 1], [4, 0]])
            >>> polygon = Polygon(coords)
            >>> polygon.area
            27.5
        '''
        # polygon area by applying the gauss equation
        area = 0.5*abs(sum(self.boundCoords[:-1, 0] * self.boundCoords[1:, 1] -
                           self.boundCoords[:-1, 1] * self.boundCoords[1:, 0]))
        setattr(self, 'area', area)
        return area

    def plot(self):
        '''Method for show the graph of the polygon.

        Examples:

            .. plot::

                from numpy import array
                from basegeometry import Polygon
                coords = array([[1, 1], [2, 5], [4.5, 6], [8, 3],
                                [7, 1], [4, 0]])
                polygon = Polygon(coords)
                polygon.plot()
        '''

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.boundCoords[:, 0], self.boundCoords[:, 1], '-k', lw=2)
        ax.plot(self.coordinates[:, 0], self.coordinates[:, 1], 'ok', ms=5)
        ax.grid(ls='--', lw=0.5)
        ax.axis('equal')
        return


# %%
'''
BSD 2 license.

Copyright (c) 2018, Universidad Nacional de Colombia, Andres Ariza-Triana
and Ludger O. Suarez-Burgoa.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
