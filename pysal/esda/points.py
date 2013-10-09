import pysal as ps
import scipy.spatial.distance as DISTANCE
from scipy.spatial import cKDTree as CKDTREE
from pysal.common import stats


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

    return [minc[0], minc[1], maxc[0], maxc[1]]


def csr_rect(n, rectangle):
    # generate n csr points in a rectangle
    left, lower, right, upper = rectangle
    x = np.random.random(n) * (right - left) + left
    y = np.random.random(n) * (upper - lower) + lower
    return np.vstack((x, y)).T


class PointsCollection(object):
    """
    Container class for point pattern analysis

    """
    def __init__(self, points, window=None):

        if type(points) != 'numpy.ndarray':
            points = np.array(points)

        self._n = len(points)
        self._area = None
        self._density = None
        self._mbr = None
        self._mtd = None
        self._kdtree = None
        window = self.set_window(window)
        self.points = points
        self.locator = ps.cg.locators.PointLocator([ps.cg.shapes.Point(p)
                                                    for p in points])

    def get_window(self):
        return self._window

    def set_window(self, window):
        if window is None:
            window = ps.cg.shapes.Rectangle(*mbr(points))
        self._window = window
        self._area = window.area  # update area and density if window changes
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
            self._mtd = np.sqrt(dx * dx + dy * dy)
        return self._mtd

    mtd = property(get_max_theoretical_distance)

    def get_kdtree(self):
        if self._kdtree is None:
            self._kdtree = CKDTREE(self.points)
        return self._kdtree

    kdtree = property(get_kdtree)


def G(points, n_bins=10, delta=1.00001, permutations=99, pct=0.05,
        max_d=None):
    """
    Minimum Nearest Neighbor distance test
    """

    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)
    d, ids = points.kdtree.query(points.kdtree.data, k=2)

    if max_d is None:
        max_d = d.max()

    width = max_d / n_bins
    bins = np.arange(0, max_d, width)
    bins[-1] *= delta
    gd = np.histogram(d[:, 1], bins)[0]

    return d, ids, gd, bins


def R(points, sampling_rate=0.10, two_tailed=True, n_samples=99):
    """
    Clark and Evan's R statistic


    Arguments
    =========

    points: nx2 array
            point data

    sampling_rate: float
                   percentage of points to take for each random sample to
                   avoid full intensive sampling
    two_tailed: Boolean
                True: null is CSR
                False: deviations from 0 indicate which tail of the
                distribution is used to calculate p-vaues
    n_samples: int
               Number of samples to take to avoid full intensive sampling

    """
    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)

    e_d = 1. / (2 * np.sqrt(points.density))  # expected value
    nn_query = points.kdtree.query(points.points, k=2)
    nn_distances = nn_query[0][:, 1]
    mean_d = nn_distances.mean()
    # Bailey and Gatrell, 1995, p 100)
    var_d = (4 - np.pi) / (4 * np.pi * points.density * points.n)

    value = mean_d / e_d
    z = mean_d - e_d
    z /= np.sqrt(var_d)
    if two_tailed:
        p = 1 - stats.norm.cdf(np.abs(z))
        p *= 2.
    else:
        if z <= 0:
            p = stats.norm.cdf(z)
        else:
            p = 1 - stats.norm.cdf(z)

    results = {}
    results['mean[nnd]'] = mean_d
    results['R'] = value
    results['E[nnd]'] = e_d
    results['V[nnd]'] = var_d
    results['z[nnd]'] = z
    results['p-value[nnd]'] = p
    results['two_tailed'] = two_tailed
    results['nn_d'] = nn_distances

    # random sampling to minimize dependence of nnd
    # See http://www.seas.upenn.edu/~ese502/#notebook Chapter 3

    m = np.int(sampling_rate * points.n)
    results['m'] = m
    var_md = (4 - np.pi) / (4 * np.pi * points.density * m)
    z_values = np.zeros((n_samples, 1))
    # random sample of nn_distances
    for sample in xrange(n_samples):
        np.random.shuffle(nn_distances)
        mean_md = nn_distances[0:m].mean()
        z_values[sample] = mean_md
    z_values = (z_values - e_d) / np.sqrt(var_md)
    results['z_values'] = z_values
    z_values_mean = z_values.mean()
    results['z_mean'] = z_values_mean
    if two_tailed:
        p_r = 1 - stats.norm.cdf(np.abs(z_values_mean))
        p_r *= 2.
    else:
        if z <= 0:
            p_r = stats.norm.cdf(z_values_mean)
        else:
            p_r = 1 - stats.norm.cdf(z_values_mean)
    results['p-value[sampled]'] = p_r
    results['m'] = m
    results['sampling_rate'] = sampling_rate

    return results


def k(points, n_bins=10, delta=1.00001, permutations=99, pct=0.05,
      max_d=None):
    """
    k-function
    """

    if not isinstance(points, PointsCollection):
        points = PointsCollection(points)
    if max_d is None:
        max_d = (points.area / 2.) ** (1 / 2.)
    width = max_d / n_bins
    bins = np.arange(0, max_d, width)
    bins[-1] *= delta  # to have max contained
    counts = np.zeros(len(bins) - 1,)
    pd = DISTANCE.pdist(points.points)  # upper triangle of distance matrix
    counts = np.histogram(pd, bins)[0] * 2
    counts = counts.cumsum()
    results = {}
    results['max_d'] = max_d
    results['width'] = width
    results['counts'] = counts
    results['bins'] = bins
    kd = counts * points.area / points.n ** 2
    results['k'] = kd
    results['l'] = np.sqrt(kd / np.pi) - bins[1:]

    if permutations:
        pCounts = np.zeros((len(bins) - 1, permutations))
        left = points.mbr.left
        lower = points.mbr.lower
        right = points.mbr.right
        upper = points.mbr.upper

        for permutation in xrange(permutations):
            pCount = np.zeros(len(bins) - 1,)

            # sample within collection window
            # for now just in the mbr, later for the window
            rpoints = csr_rect(points.n, (left, lower, right, upper))
            pd = DISTANCE.pdist(rpoints)
            pCount = np.histogram(pd, bins)[0] * 2
            pCounts[:, permutation] = pCount.cumsum()

        counts.shape = (len(counts), 1)
        pCounts = np.hstack((pCounts, counts))
        pCounts.sort(axis=1)

        # lower and upper pct envelopes
        lp = np.int(pct * (permutations + 1))
        up = np.int((1 - pct) * (permutations + 1))
        results['pCounts'] = pCounts[:, [lp, up]]

    return results


if __name__ == '__main__':

    import numpy as np

    points = np.random.random((50, 2)) * 10
    points_mbr = mbr(points)
    pc = PointsCollection(points)
    res = k(pc)
    r1 = G(points)
    # point pattern from O'Sullivan and Unwin (2003) pg 90
    points = np.array([
                      [66.22, 32.54],
                      [22.52, 22.39],
                      [31.01, 81.21],
                      [9.47, 31.02],
                      [30.78, 60.10],
                      [75.21, 58.93],
                      [79.26,  7.68],
                      [8.23, 39.93],
                      [98.73, 42.53],
                      [89.78, 42.53],
                      [65.19, 92.08],
                      [54.46, 8.48]])
    rres = R(points)
