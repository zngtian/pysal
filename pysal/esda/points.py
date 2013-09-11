import pysal as ps
import scipy.spatial.distance as DISTANCE
from  scipy.spatial import cKDTree as CKDTREE



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

def csr_rect(n, rectangle):
    # generate n csr points in a rectangle
    left, lower, right, upper = rectangle
    x = np.random.random(n) * (right - left)  + left 
    y = np.random.random(n) * (upper - lower) + lower
    return np.vstack((x,y)).T




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


def g(points, n_bins=10, delta=1.00001, permutations=99, pct=0.05,
        max_d=None):
    """
    Minimum Nearest Neighbor distance test
    """

    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)
    tree = CKDTREE(points.points)
    d, ids = tree.query(tree.data, k=2)

    if max_d is None:
        max_d = d.max()

    width = max_d / n_bins
    bins = np.arange(0, max_d, width)
    bins[-1] *= delta
    counts = np.zeros(len(bins)-1,)
    gd = np.histogram(d[:,1], bins)[0]

    return d, ids, gd, bins





def k(points, n_bins=10, delta=1.00001, permutations=99, pct=0.05, max_d =
        None):
    """
    k-function
    """

    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)
    if max_d is None:
        max_d = (points.area / 2.)**(1/2.)
    width = max_d / n_bins
    bins = np.arange(0, max_d, width)
    bins[-1] *= delta # to have max contained
    counts = np.zeros(len(bins)-1,)

    pairs = 0
    maxd = 0
    pd = DISTANCE.pdist(points.points) #upper triangle of distance matrix
    counts = np.histogram(pd, bins)[0] * 2

    counts = counts.cumsum()
    results = {}
    results['max_d'] = max_d
    results['width'] = width
    results['counts'] = counts
    results['bins'] = bins
    kd = counts * points.area  /  points.n**2
    results['k'] = kd
    results['l'] = np.sqrt(kd/np.pi) - bins[1:]


    if permutations:
        pCounts = np.zeros((len(bins)-1, permutations))
        left = points.mbr.left
        lower = points.mbr.lower
        right = points.mbr.right
        upper = points.mbr.upper

        for permutation in xrange(permutations):
            pCount = np.zeros(len(bins)-1,)

            # sample within collection window
            # for now just in the mbr, later for the window
            rpoints = csr_rect(points.n, (left, lower, right, upper)) 
            pd = DISTANCE.pdist(rpoints)
            pCount = np.histogram(pd, bins)[0] * 2
            pCounts[:,permutation] = pCount.cumsum()

        counts.shape = (len(counts),1)
        pCounts = np.hstack((pCounts,counts))
        pCounts.sort(axis=1)

        # lower and upper pct envelopes
        lp = np.int(pct * (permutations+1))
        up = np.int((1-pct) * (permutations+1))
        results['pCounts'] = pCounts[:, [lp,up]]

    return results



if __name__ == '__main__':

    import numpy as np

    points = np.random.random((50,2))*10

    points_mbr = mbr(points)

    pc = PointsCollection(points)

    res = k(pc)



    r1 = g(points)

    

