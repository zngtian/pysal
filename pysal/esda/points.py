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
    maxc = points.max(axis=1)
    minc = points.min(axis=1)

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
        window = self.set_window(window)
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




    


if __name__ == '__main__':

    import numpy as np

    points = np.random.random((50,2))

    points_mbr = mbr(points)

    pc = PointsCollection(points)




