import pysal as ps



def mbr(points):

    """
    Minimum bounding rectangle for a two dimensional point set


    Arguments
    =========

    points: array (nx2)
            x,y coordinates for n points

    Returns
    =======
    rect: array (4x1)
          [left,bottom, right, top]
    """

    #pl = ps.cg.locators.PointLocator([ps.cg.shapes.Point(p) for p in points])
    maxc = points.max(axis=0)
    minc = points.min(axis=0)

    return [minc[0],minc[1], maxc[0], maxc[1]] 


class PointsCollection(object):
    """
    Container class for point pattern analysis
    
    """
    def __init__(self, points, window=None):

        self._n = len(points)
        self._area = None
        self._density = None
        self._mbr = None
        self._mtd = None
        window = self.set_window(window)
        self.points = points
        self.locator = ps.cg.locators.PointLocator([ps.cg.shapes.Point(p) for
            p in points])

    def get_window(self):
        return self._window

    def set_window(self, window):
        if window is None:
            window = ps.cg.shapes.Rectangle(*mbr(points))
        self._window = window
        self._area = window.area # update area and density if window changes
        self._density = self._n / self._area

    window = property(get_window, set_window)

    def get_n(self):
        return self._n

    n = property(get_n)

    def get_area(self):
        return self._area

    area = property(get_area)

    def get_density(self):
        return self._density

    density = property(get_density)

    def get_mbr(self):
        if self._mbr is None:
            self._mbr = ps.cg.shapes.Rectangle(*mbr(points))
        return self._mbr

    mbr = property(get_mbr)

    def get_max_theoretical_distance(self):
        if self._mtd is None:
            if self._mbr is None:
                self._mbr = ps.cg.shapes.Rectangle(*mbr(points))
            dx = self._mbr.right - self._mbr.left
            dy = self._mbr.upper - self._mbr.lower
            self._mtd = np.sqrt(dx*dx + dy*dy)
        return self._mtd

    mtd = property(get_max_theoretical_distance)



def k(points, n_bins=10, delta=1.00001):
    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)
    width = points.mtd / n_bins
    print points.mtd
    bins = np.arange(0, pc.mtd+width, width)
    bins[-1] *= delta # to have max contained
    counts = np.zeros(len(bins)-1,)
    print len(counts)

    pairs = 0
    maxd = 0
    # brute force distances
    for i in xrange(points.n-1):
        p_i = points.points[i]
        dx = points.points[i,0] - points.points[:,0]
        dy = points.points[i,1] - points.points[:,1]
        d_i = np.sqrt(dx*dx + dy*dy)
        print len(d_i)
        pairs += len(d_i)
        counts_i = np.histogram(d_i, bins)[0]
        counts = counts + counts_i
        maxd = max(maxd, max(d_i))
        #print len(counts_i)

    print pairs
    print maxd



    return bins, points.mtd, counts




if __name__ == '__main__':

    import numpy as np

    points = np.random.random((50,2))*10

    points_mbr = mbr(points)

    pc = PointsCollection(points)




