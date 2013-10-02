import math
import scipy
import numpy
from pysal.cg.shapes import Rectangle, Point, LineSegment
from pysal.cg.standalone import get_segment_point_dist, get_bounding_box
import random
import time

__all__ = ["SegmentGrid", "SegmentLocator",
           "Polyline_Shapefile_SegmentLocator", "random_segments"]
DEBUG = False


class BruteSegmentLocator(object):
    def __init__(self, segments):
        self.data = segments
        self.n = len(segments)

    def nearest(self, pt):
        d = self.data
        distances = [get_segment_point_dist(
            d[i], pt)[0] for i in xrange(self.n)]
        return numpy.argmin(distances)


class SegmentLocator(object):
    def __init__(self, segments, nbins=500):
        self.data = segments
        if hasattr(segments, 'bounding_box'):
            bbox = segment.bounding_box
        else:
            bbox = get_bounding_box(segments)
        self.bbox = bbox
        res = max((bbox.right - bbox.left), (bbox.upper -
                                             bbox.lower)) / float(nbins)
        self.grid = SegmentGrid(bbox, res)
        for i, seg in enumerate(segments):
            self.grid.add(seg, i)

    def nearest(self, pt):
        d = self.data
        possibles = self.grid.nearest(pt)
        distances = [get_segment_point_dist(d[i], pt)[0] for i in possibles]
        #print "possibles",possibles
        #print "distances",distances
        #print "argmin", numpy.argmin(distances)
        return possibles[numpy.argmin(distances)]


    def intersections(self):
        """
        Find all segment intersections
 
        Example
        =======
        >>> segments = [[12,7,13,11],
        ...            [3,3,8,9],
        ...         [3,11, 9,10],
        ...         [4,10, 5,12],
        ...         [6,9, 9,7]]
        >>> segs = []
        >>> for segment in segments:
        ...     head = Point((segment[0], segment[1]))
        ...     tail = Point((segment[2], segment[3]))
        ...     segs.append(LineSegment(head, tail))
        >>> locator = SegmentLocator(segs)
        >>> locator.intersections()
        [(1, 4), (2, 3)]
        >>>
        """
        return self._intersections_bf()



    def _intersections_bf(self):
        """
        Brute force segment intersection

        """
        n_segments = len(self.data)
        hits = []
        for i in xrange(n_segments-1):
            seg_i = self.data[i]
            for j in xrange(i+1, n_segments):
                if seg_i.intersect(self.data[j]):
                    hits.append((i,j))
        return hits



    def intersections_ps(self):
        """
        Plane sweep intersection detection

        XXX INCOMPLETE
        """

        points2Segments = {}
        eventPoints = set()
        for i, segment in enuemrate(self.data):
            coordsHead = segment.p1
            coordsTail = segment.p2
            if coordsHead not in points2Segments:
                points2Segments[coordsHead] = []
            if coordsTail not in points2Segments:
                points2Segments[coordsHead] = []
            points2Segments[coordsHead].append(i) 
            points2Segments[coordsTail].append(i) 
            eventPoints.add(coordsHead)
            eventPoints.add(coordsTail)
        eventPoints = list(s).sort()

        segments = self.data

        # first event
        que = Node(eventPoints[0])
        firstSegments = points2Segments[eventPoints[0]]
        status = Node(firstSegments[0])

        # populate the rest of the que
        for pnt in eventPoints[1:]:
            que.insert(pnt)

        results = {}
        results['que'] = que
        results['status'] = status
        results['segments' ] = segments
        return results










class Polyline_Shapefile_SegmentLocator(object):
    def __init__(self, shpfile, nbins=500):
        self.data = shpfile
        bbox = Rectangle(*shpfile.bbox)
        res = max((bbox.right - bbox.left), (bbox.upper -
                                             bbox.lower)) / float(nbins)
        self.grid = SegmentGrid(bbox, res)
        for i, polyline in enumerate(shpfile):
            for p, part in enumerate(polyline.segments):
                for j, seg in enumerate(part):
                    self.grid.add(seg, (i, p, j))

    def nearest(self, pt):
        d = self.data
        possibles = self.grid.nearest(pt)
        distances = [get_segment_point_dist(
            d[i].segments[p][j], pt)[0] for (i, p, j) in possibles]
        #print "possibles",possibles
        #print "distances",distances
        #print "argmin", numpy.argmin(distances)
        return possibles[numpy.argmin(distances)]


class SegmentGrid(object):
    """
    Notes:
        SegmentGrid is a low level Grid class.
        This class does not maintain a copy of the geometry in the grid.
        It returns only approx. Solutions.
        This Grid should be wrapped by a locator.
    """
    def __init__(self, bounds, resolution):
        """
        Returns a grid with specified properties.

        __init__(Rectangle, number) -> SegmentGrid

        Parameters
        ----------
        bounds      : the area for the grid to encompass
        resolution  : the diameter of each bin

        Examples
        --------
        TODO: complete this doctest
        >>> g = SegmentGrid(Rectangle(0, 0, 10, 10), 1)
        """
        if resolution == 0:
            raise Exception('Cannot create grid with resolution 0')
        self.res = resolution
        self.hash = {}
        self._kd = None
        self._kd2 = None
        self._hashKeys = None
        self.x_range = (bounds.left, bounds.right)
        self.y_range = (bounds.lower, bounds.upper)
        try:
            self.i_range = int(math.ceil((self.x_range[1] -
                                          self.x_range[0]) / self.res)) + 1
            self.j_range = int(math.ceil((self.y_range[1] -
                                          self.y_range[0]) / self.res)) + 1
            self.mask = numpy.zeros((self.i_range, self.j_range), bool)
            self.endMask = numpy.zeros((self.i_range, self.j_range), bool)
        except Exception:
            raise Exception('Invalid arguments for SegmentGrid(): (' + str(self.x_range) + ', ' + str(self.y_range) + ', ' + str(self.res) + ')')
    @property
    def hashKeys(self):
        if self._hashKeys == None:
            self._hashKeys = numpy.array(self.hash.keys(),dtype=float)
        return self._hashKeys

    @property
    def kd(self):
        if self._kd == None:
            self._kd = scipy.spatial.cKDTree(self.hashKeys)
        return self._kd

    @property
    def kd2(self):
        if self._kd2 == None:
            self._kd2 = scipy.spatial.KDTree(self.hashKeys)
        return self._kd2

    def in_grid(self, loc):
        """
        Returns whether a 2-tuple location _loc_ lies inside the grid bounds.
        """
        return (self.x_range[0] <= loc[0] <= self.x_range[1] and
                self.y_range[0] <= loc[1] <= self.y_range[1])

    def _grid_loc(self, loc):
        i = int((loc[0] - self.x_range[0]) / self.res)  # floored
        j = int((loc[1] - self.y_range[0]) / self.res)  # floored
        #i = min(self.i_range-1, max(int((loc[0] - self.x_range[0])/self.res), 0))
        #j = min(self.j_range-1, max(int((loc[1] - self.y_range[0])/self.res), 0))
        #print "bin:", loc, " -> ", (i,j)
        return (i, j)

    def _real_loc(self, grid_loc):
        x = (grid_loc[0] * self.res) + self.x_range[0]
        y = (grid_loc[1] * self.res) + self.y_range[0]
        return x, y

    def bin_loc(self, loc, id):
        grid_loc = self._grid_loc(loc)
        if grid_loc not in self.hash:
            self.hash[grid_loc] = set()
            self.mask[grid_loc] = True
        self.hash[grid_loc].add(id)
        return grid_loc

    def add(self, segment, id):
        """
        Adds segment to the grid.

        add(segment, id) -> bool

        Parameters
        ----------
        id -- id to be stored int he grid.
        segment -- the segment which identifies where to store 'id' in the grid.

        Examples
        --------
        >>> g = SegmentGrid(Rectangle(0, 0, 10, 10), 1)
        >>> g.add(LineSegment(Point((0.2, 0.7)), Point((4.2, 8.7))), 0)
        True
        """
        if not (self.in_grid(segment.p1) and self.in_grid(segment.p2)):
            raise Exception('Attempt to insert item at location outside grid bounds: ' + str(segment))
        i, j = self.bin_loc(segment.p1, id)
        I, J = self.bin_loc(segment.p2, id)
        self.endMask[i, j] = True
        self.endMask[I, J] = True

        bbox = segment.bounding_box
        left = bbox.left
        lower = bbox.lower
        res = self.res
        line = segment.line
        tiny = res / 1000.
        for i in xrange(1 + min(i, I), max(i, I)):
            #print 'i',i
            x = self.x_range[0] + (i * res)
            y = line.y(x)
            self.bin_loc((x - tiny, y), id)
            self.bin_loc((x + tiny, y), id)
        for j in xrange(1 + min(j, J), max(j, J)):
            #print 'j',j
            y = self.y_range[0] + (j * res)
            x = line.x(y)
            self.bin_loc((x, y - tiny), id)
            self.bin_loc((x, y + tiny), id)
        self._kd = None
        self._kd2 = None
        return True

    def remove(self, segment):
        self._kd = None
        self._kd2 = None
        pass

    def nearest(self, pt):
        """
        Return a set of ids.

        The ids identify line segments within a radius of the query point.
        The true nearest segment is guaranteed to be within the set.

        Filtering possibles is the responsibility of the locator not the grid.
        This means the Grid doesn't need to keep a reference to the underlying segments,
        which in turn means the Locator can keep the segments on disk.

        Locators can be customized to different data stores (shape files, SQL, etc.)
        """
        grid_loc = numpy.array(self._grid_loc(pt))
        possibles = set()

        if DEBUG:
            print "in_grid:", self.in_grid(pt)
            i = pylab.matshow(self.mask, origin='lower',
                              extent=self.x_range + self.y_range, fignum=1)
        # Use KD tree to search out the nearest filled bin.
        # it may be faster to not use kdtree, or at least check grid_loc first
        # The KD tree is build on the keys of self.hash, a dictionary of stored bins.
        dist, i = self.kd.query(grid_loc, 1)

        ### Find non-empty bins within a radius of the query point.
        # Location of Q point
        row, col = grid_loc
        # distance to nearest filled cell +2.
        # +1 returns inconsistent results (compared to BruteSegmentLocator)
        # +2 seems to do the trick.
        radius = int(math.ceil(dist)) + 2
        if radius < 30:
            a, b = numpy.ogrid[-radius:radius + 1, -radius:radius +
                               1]   # build square index arrays centered at 0,0
            index = a ** 2 + b ** 2 <= radius ** 2                        # create a boolean mask to filter indicies outside radius
            a, b = index.nonzero()
                # grad the (i,j)'s of the elements within radius.
            rows, cols = row + a - radius, col + b - radius                   # recenter the (i,j)'s over the Q point
            #### Filter indicies by bounds of the grid.
            ### filters must be applied one at a time
            ### I havn't figure out a way to group these
            filter = rows >= 0
            rows = rows[filter]
            cols = cols[filter]  # i >= 0
            filter = rows < self.i_range
            rows = rows[filter]
            cols = cols[filter]  # i < i_range
            filter = cols >= 0
            rows = rows[
                filter]
            cols = cols[filter]  # j >= 0
            filter = cols < self.j_range
            rows = rows[
                filter]
            cols = cols[filter]  # j < j_range
            if DEBUG:
                maskCopy = self.mask.copy().astype(float)
                maskCopy += self.endMask.astype(float)
                maskCopy[rows, cols] += 1
                maskCopy[row, col] += 3
                i = pylab.matshow(maskCopy, origin='lower', extent=self.x_range + self.y_range, fignum=1)
                #raw_input('pause')
            ### All that was just setup for this one line...
            idx = self.mask[rows, cols].nonzero()[0] # Filter out empty bins.
            rows, cols = rows[idx], cols[idx]        # (i,j)'s of the filled grid cells within radius.

            for t in zip(rows, cols):
                possibles.update(self.hash[t])

            if DEBUG:
                print "possibles", possibles
        else:
        ### The old way...
        ### previously I was using kd.query_ball_point on, but the performance was terrible.
            I = self.kd2.query_ball_point(grid_loc, radius)
            for i in I:
                t = tuple(self.kd.data[i])
                possibles.update(self.hash[t])
        return list(possibles)


def random_segments(n):
    segs = []
    for i in xrange(n):
        a, b, c, d = [random.random() for x in [1, 2, 3, 4]]
        seg = LineSegment(Point((a, b)), Point((c, d)))
        segs.append(seg)
    return segs


def random_points(n):
    return [Point((random.random(), random.random())) for x in xrange(n)]


def test_combo(bins, segments, qpoints):
    G = SegmentLocator(segments, bins)
    G2 = BruteSegmentLocator(segs)
    for pt in qpoints:
        a = G.nearest(pt)
        b = G2.nearest(pt)
        if a != b:
            print a, b, a == b
            global DEBUG
            DEBUG = True
            a = G.nearest(pt)
            print a
            a = segments[a]
            b = segments[b]
            print "pt to a (grid)", get_segment_point_dist(a, pt)
            print "pt to b (brut)", get_segment_point_dist(b, pt)
            raw_input()
            pylab.clf()
            DEBUG = False


def test_brute(segments, qpoints):
    t0 = time.time()
    G2 = BruteSegmentLocator(segs)
    t1 = time.time()
    print "Created Brute in %0.4f seconds" % (t1 - t0)
    t2 = time.time()
    q = map(G2.nearest, qpoints)
    t3 = time.time()
    print "Brute Found %d matches in %0.4f seconds" % (len(qpoints), t3 - t2)
    print "Total Brute Time:", t3 - t0
    print
    return q


def test_grid(bins, segments, qpoints, visualize=False):
    t0 = time.time()
    G = SegmentLocator(segments, bins)
    t1 = time.time()
    G.grid.kd
    t2 = time.time()
    print "Created Grid in %0.4f seconds" % (t1 - t0)
    print "Created KDTree in %0.4f seconds" % (t2 - t1)
    if visualize:
        i = pylab.matshow(G.grid.mask, origin='lower',
                          extent=G.grid.x_range + G.grid.y_range)

    t2 = time.time()
    q = map(G.nearest, qpoints)
    t3 = time.time()
    print "Grid Found %d matches in %0.4f seconds" % (len(qpoints), t3 - t2)
    print "Total Grid Time:", t3 - t0
    qps = len(qpoints) / (t3 - t2)
    print "q/s:", qps
    #print
    return qps


def binSizeTest():
    q = 100
    minN = 1000
    maxN = 10000
    stepN = 1000
    minB = 250
    maxB = 2000
    stepB = 250
    sizes = range(minN, maxN, stepN)
    binSizes = range(minB, maxB, stepB)
    results = numpy.zeros((len(sizes), len(binSizes)))
    for row, n in enumerate(sizes):
        segs = random_segments(n)
        qpts = random_points(q)
        for col, bins in enumerate(binSizes):
            print "N, Bins:", n, bins
            qps = test_grid(bins, segs, qpts)
            results[row, col] = qps
    return results


# prototype of binary tree for segment intersection tests

class Node:
    """ """
    def __init__(self, data):
        self.left = None
        self.right = None
        self.data = data

    def insert(self, data):

        if data < self.data:
            if self.left is None:
                self.left = Node(data)
            else:
                self.left.insert(data)
        else:
            if self.right is None:
                self.right = Node(data)
            else:
                self.right.insert(data)

    def lookup(self, data, parent=None):
        if data < self.data:
            if self.left is None:
                return None, None
            return self.left.lookup(data, self)
        else:
            return self, parent

    def children_count(self):
        cnt = 0
        if self.left:
            cnt += 1
        if self.right:
            cnt += 1
        return cnt

    def delete(self, data):
        node, parent = self.lookup(data)
        if node is not None:
            children_count = node.children_count()
            if children_count == 0:
                if parent.left is node:
                    parent.left = None
                else:
                    parent.right = NOne
                del node
            elif children_count == 1:
                if node.left:
                    n = node.left
                else:
                    n = node.right
                if parent:
                    if parent.left is node:
                        parent.left = n
                    else:
                        parent.right = n
                del node
            else:
                parent = node
                successor = node.right
                while successor.left:
                    parent = successor
                    successor = successor.left
                node.data = successor.data
                if parent.left == successor:
                    parent.left = successor.right
                else:
                    parent.right = successor.right

    def print_tree(self):
        if self.left:
            self.left.print_tree()
        print self.data
        if self.right:
            self.right.print_tree()


def ccw(p0,p1,p2):
    """
    tests if in going from p0->p1->p2 direction is counter-clocwise


    Parameters
    ==========
    p0: point (tuple or list)

    p1: point (tuple or list)

    p2: point (tuple or list)





    Returns
    =======
    
    0: when p2 is on the line between p0 and p1
    1: when p1 is on the line between p0 and p2 or p2 is to the left of (p0,p1)
    -1: when p0 is between p1 and p2


    Notes
    =====
    From Sedgwick, R. (1992) Algorithms in C++. Addison Wesley.  p 350.

    """
    x=0
    y=1
    dx1 = p1[x] - p0[x]
    dy1 = p1[y] - p0[y]
    dx2 = p2[x] - p0[x]
    dy2 = p2[y] - p0[y]

    if dx1*dy2 > dy1*dx2:
        return 1
    if dx1*dy2 < dy1*dx2:
        return -1
    if (( dx1*dx2 <0) or (dy1*dy2 < 0)):
        return -1
    if ((dx1*dx1+dy1*dy1) < (dx2*dx2+dy2*dy2)):
        return 1

    return 0


def intersect(seg1, seg2):
    """
    Test if two line segments intersect

    Parameters
    ==========

    seg1: segment (pair of points)

    seg2: segment (pair of points)


    Returns
    =======

    1: if segments intersect

    0: if segments do not intersect

    >>> seg1 = [ (0,0), (5,4) ]
    >>> seg2 = [ (0,6), (7,0) ]
    >>> seg3 = [ (0,6), (0,10) ]
    >>> intersect(seg1, seg2)
    1
    >>> intersect(seg2, seg3)
    1
    >>> intersect(seg1, seg3)
    0

    """
    a = ccw(seg1[0], seg1[1], seg2[0])  
    b = ccw(seg1[0], seg1[1], seg2[1])
    c = ccw(seg2[0], seg2[1], seg1[0])
    d = ccw(seg2[0], seg2[1], seg1[1])
    return (a*b <= 0) * (c*d <= 0)



def pointOnSegment(point, segment):
    a,b = segment
    if ccw(a,b,point) == 0:
        return 1
    return 0


def intersectionPoint(seg1, seg2):
    p1,p2 = seg1
    p3,p4 = seg2
    x=0
    y=1
    d = 1.0 * (p1[x]-p2[x])*(p3[y]-p4[y]) - \
            (p1[y]-p2[y])*(p3[x]-p4[x])

    if d == 0:
        return None
    
    xi = ((p3[x]-p4[x])*(p1[x]*p2[y]-p1[y]*p2[x]) - \
            (p1[x]-p2[x]) * (p3[x]*p4[y]-p3[y]*p4[x])) / d
    yi = ((p3[y]-p4[y])*(p1[x]*p2[y]-p1[y]*p2[x]) - \
            (p1[y]-p2[y]) * (p3[x]*p4[y]-p3[y]*p4[x])) / d
    return (xi,yi)




from heapq import heappush, heappop, heapify


class BinaryHeap:
    """ """
    def __init__(self):
        self.heap_list = []

    def insert(self, element):
        self.heap_list.append(element)
        heapify(self.heap_list)

    def remove(self, element):
        if element in self.heap_list:
            self.heap_list.remove(element)
            heapify(self.heap_list)

    def popMin(self):
        return heappop(self.heap_list)

    def order(self):
        tmp = [ self.popMin() for i in range(len(self.heap_list)) ]
        self.heap_list = tmp[:]
        heapify(self.heap_list)
        return tmp



if __name__ == '__main__':
    import pylab
    pylab.ion()

    n = 100
    q = 1000

    t0 = time.time()
    segs = random_segments(n)
    t1 = time.time()
    qpts = random_points(q)
    t2 = time.time()
    print "segments:", t1 - t0
    print "points:", t2 - t1
    #test_brute(segs,qpts)
    #test_grid(50, segs, qpts)

    SG = SegmentLocator(segs)
    grid = SG.grid


    #bh = BinaryHeap()
    #alist = [9, 6, 5, 2, 3]
    #bh.buildHeap(alist)

    heap = []

    # line segments
    segs = [
            [(3,11), (9,10)],
            [(4,10), (5,12) ],
            [(6,9), (9,7)],
            [(3,3), (8,9)],
            [(12,7), (12,11)] ]
    pnts = []
    s2p ={}
    p2s = {}
    l2s = {} # left to segment bridge
    r2s = {} # right to segment bridge
    i2segs = {}
    for i,seg in enumerate(segs):
        seg.sort() # sort points by x coord
        h,t = seg
        i2segs[i] = seg
        pnts.extend(seg)
        if h not in p2s:
            p2s[h] = []
        if t not in p2s:
            p2s[t] = []
        if h not in l2s:
            l2s[h] = []
        l2s[h].append(i)
        if t not in r2s:
            r2s[t] = []
        r2s[t].append(i)
        p2s[h].append(i)
        p2s[t].append(i)
        s2p[i] = h,t

    heapify(pnts)
    #bh = BinaryHeap()
    #bh.buildHeap(pnts)

    status = []

    intersections =  {} 
    left = {} 
    right = {} 
    y2p = {}


    onstatus = set()
    s = BinaryHeap() # status
    q = BinaryHeap() # event que
    for pnt in pnts:
        q.insert(pnt)

    while q.heap_list:
        event_point = q.popMin()
        print 'event_point: ', event_point
        entering = set() 
        if event_point in l2s:
            for segment in l2s[event_point]:
                entering.add(segment)
        leaving = set()
        if event_point in r2s:
            for segment in r2s[event_point]:
                leaving.add(segment)
        print 'leaving: ', leaving
        print 'entering: ', entering
        # get segments on status that contain event point
        ep_s = [ i for i in onstatus if pointOnSegment(event_point, segs[i]) ]
        l_u_s = leaving.union(ep_s)
        e_u_s = entering.union(ep_s)
        uels = l_u_s.union(e_u_s)
        print 'uels: ', uels
        print 'lus: ', l_u_s
        print 'eus: ', e_u_s
        if len(uels) > 1:
            print event_point, ' is an intersecting point'
            print uels
        for segment in l_u_s:
            onstatus.remove(segment)
        for segment in e_u_s:
            onstatus.add(segment)

        status = BinaryHeap()
        for i in onstatus:
            segment = i2segs[i]
            x,y = segment[0]
            key = (y,x,i)
            status.insert(key)

        
        if e_u_s:
            print e_u_s
            # get order of segments of status
            order = [i[-1] for i in status.order()]
            # find bottom most segment of those entering on p or containing p
            # on status
            bottom_most = [ i for i in order if i in e_u_s][0] 
            if bottom_most:
                # check if bottom_most has a neighbor on the status below it
                bm_index = order.index(bottom_most)
                if bm_index > 0:
                    neighbor = order[bm_index-1]
                    si = i2segs[bottom_most]
                    sj = i2segs[neighbor]
                    if intersect(si, sj):
                        pair = [bottom_most,neighbor]
                        pair.sort()
                        ip = intersectionPoint(si,sj)
                        intersections[tuple(pair)] = ip
                        if ip[1] > event_point[1] or (ip[1]==event_point[1]
                                and ip[0] > event_point[0]):
                            q.remove(ip) # no dupes
                            q.insert(ip)


            # find top most segment of those entering on the status or contain
            # p on the status
            top_most = [ i for i in order if i in e_u_s][-1]
            if top_most:
                # check if top_most has a neighbor on the status above it
                tm_index = order.index(top_most)
                if order[-1] != order[tm_index]:
                    neighbor = order[tm_index + 1]
                    si = i2segs[top_most]
                    sj = i2segs[neighbor]
                    if intersect(si, sj):
                        pair = [top_most, neighbor]
                        pair.sort()
                        ip = intersectionPoint(si,sj)
                        intersections[tuple(pair)] = ip
                        if ip[1] > event_point[1] or (ip[1]==event_point[1]
                                and ip[0] > event_point[0]):
                            q.remove(ip) # no dupes
                            q.insert(ip)
        else:
            # find new event using bottom and top neighbors of p on status
            order = [ i[-1] for i in status.order() ]
            ys = numpy.array([i2segs[i][0][1] for i in order])
            bottom = numpy.nonzero(ys<=event_point[1])[0]
            top = numpy.nonzero(ys > event_point[1])[0]
            if bottom and top:
                si = i2segs[bottom[-1]]
                sj = i2segs[top[0]]
                if si != sj and intersect(si,sj):
                    ip = intersectionPoint(si,sj)
                    if ip[1] > event_point[1] or (ip[1]==event_point[1]
                                and ip[0] > event_point[0]):
                            q.remove(ip) # no dupes
                            q.insert(ip)


        print 'status: ', onstatus
        print '\n\n'

        #raw_input('here')


