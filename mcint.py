from random import uniform

class Vect():
    """an N-d vector."""

    def __init__(self, coords):
        self.coords = coords
        self.N = len(coords)

    def add(self, x):
        if x.N == self.N:
            return Vect([self.coords[i] + x.coords[i] for i in range(self.N)])
        else:
            print "Error: cannot add vectors of different dimension."

    def dot(self, x):
        return sum([self.coords[i]*x.coords[i] for i in range(self.N)])

    def scale(self, k):
        return Vect([self.coords[i]*k for i in range(self.N)])

    def project(self, x):
        return x.scale(self.dot(x)/float(x.dot(x)))

class Domain():
    """a rectangular domain for N-d monte carlo integration. initializes to an empty 1x1x1 box with default args."""

    def __init__(self, minimum=[0,0,0], maximum=[1,1,1]):
        if len(minimum) == len(maximum):
            self.regions = []
            self.min = minimum
            self.max = maximum
            self.N = len(minimum)
        else:
            print "Error: input args must have same dimension."

    def addRegion(self, region):
        if self.N == region.N:
            self.regions.append(region)
        else:
            print "Error: region dimension does not match domain dimension."

    def containsPoint(self, point):
        for region in self.regions:
            if region.containsPoint(point):
                return True
        return False

    def generatePoint(self):
        return [uniform(self.min[i], self.max[i]) for i in range(self.N)]

    def integrate(self, trials, func=(lambda point: 1)):
        """this will be updated to use the VEGAS algorithm for greater efficiency / reduction of variance."""
        integral = 0
        i = 0
        V = 1
        while i < self.N:
            V = V*(self.max[i] - self.min[i])
            i += 1
        count = 0
        dV = V/float(trials)
        while count < trials:
            point = self.generatePoint()
            if self.containsPoint(point):
                integral += func(point)
            count += 1
        return integral*dV

class Box():
    """an N-d box region. initializes under default args to 1x1x1."""

    def __init__(self, vertex=[0,0,0], oppositeVertex=[1,1,1]):
        if len(vertex) == len(oppositeVertex):
            self.N = len(vertex)
            self.min = [min(vertex[i], oppositeVertex[i]) for i in range(self.N)]
            self.max = [max(vertex[i], oppositeVertex[i]) for i in range(self.N)]
        else:
            print "Error: input args must have same dimension."

    def containsPoint(self, point):
        i = 0
        while i < self.N:
            if not (self.min[i] < point[i] < self.max[i]):
                return False
            i += 1
        return True

class Sphere():
    """an N-d spherical region. initializes under default args to 3d sphere centered at the origin with radius 1."""

    def __init__(self, center=[0,0,0], radius=1):
        self.center = center
        self.radius = radius
        self.N = len(self.center)

    def containsPoint(self, point):
        delta = Vect([point[i] - self.center[i] for i in range(self.N)])
        if delta.dot(delta) < self.radius**2:
            return True
        return False

class Cone():
    """an N-d conical region. initializes under default args to 3d cone with radius 1, centered at origin with tip at [0,0,1]."""

    def __init__(self, center=[0,0,0], tip=[0,0,1], radius=1):
        if len(center) == len(tip):
            self.center = Vect(center)
            self.tip = Vect(tip)
            self.axis = self.tip.add(self.center.scale(-1))
            self.N = len(center)
            self.R = radius
        else:
            print "Error: center and tip must have same dimension."

    def containsPoint(self, point):
        X = Vect(point)
        X_0 = X.add(self.center.scale(-1))                         #vector from center to point
        X_a = X_0.project(self.axis)                               #axial projection of X_0
        X_r = X_0.add(X_a.scale(-1))                               #radial projection of X_0
        s = X_a.dot(self.axis)/float(self.axis.dot(self.axis))
        r = (1 - s)*self.R
        if not (0 < s < 1):                                              #axial condition
            return False
        elif not (X_r.dot(X_r) < r**2):                                  #radial condition
            return False
        return True

if __name__ == '__main__':
    """compute the volume of a 4-d cone."""
    D = Domain([0,0,0,0], [1,1,1,1])
    C = Cone([0.5,0.5,0.5,0], [0.5,0.5,0.5,1], 0.5)
    D.addRegion(C)
    print D.integrate(100000)
