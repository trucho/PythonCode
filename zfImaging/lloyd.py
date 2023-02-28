import numpy as np
import scipy.spatial
import pandas as pd
import os, warnings, time , sys
from skimage.io import imread
from tqdm.notebook import tqdm, trange
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from skimage.measure import label, regionprops, regionprops_table
from plotParams import *

baseColor = '#ffffff'

class Field():
    '''
    Create a Voronoi map that can be used to run Lloyd
    relaxation on an array of 2D points. For background,
    see: https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
    Created by https://github.com/duhaime/lloyd/
    Modified by Angueyra (Jan 2023):
        modified bounding box to be hard coded for 1024x1024 images, instead of self-bounding by the data itself
        added clipping of regions during relaxation by bounding box using Sutherland-Hodgman algorithm
        (see: https://github.com/mhdadk/sutherland-hodgman)
    '''

    def __init__(self, *args, **kwargs):
        '''
        Store the points and bounding box of the points to which
        Lloyd relaxation will be applied.
        @param np.array `arr`: a numpy array with shape n, 2, where n
          is the number of 2D points to be moved
        @param float `epsilon`: the delta between the input point
          domain and the pseudo-points used to constrain the points
        '''
        arr = args[0]
        if not isinstance(arr, np.ndarray) or arr.shape[1] != 2:
            raise Exception('Please provide a numpy array with shape n,2')
        self.points = arr
            # find the bounding box of the input data
            # self.domains = self.get_domains(arr) # from original code: bounds are set by data itself
        # make bounding box
        self.buffer = kwargs.get('buffer', 0)
        self.domains = self.get_domains1024(self.buffer) # use for 1024 x 1024 images from confocal imaging
        # ensure no two points have the exact same coords
        self.jitter_points()
        self.bb_points = self.get_bb_points(arr)
        self.constrain = kwargs.get('constrain', True)
        self.build_voronoi()

    def jitter_points(self, scalar=.000000001):
        '''
        Ensure no two points have the same coords or else the number
        of regions will be less than the number of input points
        '''
        while self.points_contain_duplicates():
            positive = np.random.rand( len(self.points), 2 ) * scalar
            negative = np.random.rand( len(self.points), 2 ) * scalar
            self.points = self.points + positive - negative
            self.constrain_points()

    def constrain_points(self):
        '''
        Update any points that have drifted beyond the boundaries of this space
        '''
        for point in self.points:
            if point[0] < self.domains['x']['min']: point[0] = self.domains['x']['min']
            if point[0] > self.domains['x']['max']: point[0] = self.domains['x']['max']
            if point[1] < self.domains['y']['min']: point[1] = self.domains['y']['min']
            if point[1] > self.domains['y']['max']: point[1] = self.domains['y']['max']


    def get_domains(self, arr):
        '''
        Return an object with the x, y domains of `arr`
        '''
        x = arr[:, 0]
        y = arr[:, 1]
        return {
          'x': {
            'min': min(x),
            'max': max(x),
          },
          'y': {
            'min': min(y),
            'max': max(y),
          }
        }

    def get_domains1024(self,buffer):
        '''
        Harcoding domains based on 1024 x 1024 image
        added by Angueyra
        '''
        return {
          'x': {
            'min': 0-buffer,
            'max': 1024+buffer,
          },
          'y': {
            'min': 0-buffer,
            'max': 1024+buffer,
          }
        }

    def get_bb_points(self, arr):
        '''
        Given an array of 2D points, return the four vertex bounding box
        '''
        return np.array([
          [self.domains['x']['min'], self.domains['y']['min']],
          [self.domains['x']['max'], self.domains['y']['min']],
          [self.domains['x']['min'], self.domains['y']['max']],
          [self.domains['x']['max'], self.domains['y']['max']],
        ])

    def get_bb_pointsClockwise(self):
        '''
        Return the four vertex bounding box
        added by Angueyra
        '''
        return np.array([
          [self.domains['x']['min'], self.domains['y']['min']],
          [self.domains['x']['min'], self.domains['y']['max']],
          [self.domains['x']['max'], self.domains['y']['max']],
          [self.domains['x']['max'], self.domains['y']['min']],
        ])


    def build_voronoi(self):
        '''
        Build a voronoi map from self.points. For background on
        self.voronoi attributes, see: https://docs.scipy.org/doc/scipy/
          reference/generated/scipy.spatial.Voronoi.html
        '''
        # build the voronoi tessellation map
        try:
            self.voronoi = scipy.spatial.Voronoi(self.points, qhull_options='Qbb Qc Qx')
        except:
            self.voronoi = scipy.spatial.Voronoi(self.points, qhull_options='QJ')

        # constrain voronoi vertices within bounding box
        if self.constrain:
            for idx, vertex in enumerate(self.voronoi.vertices):
                x, y = vertex
                if x < self.domains['x']['min']:
                    self.voronoi.vertices[idx][0] = self.domains['x']['min']
                if x > self.domains['x']['max']:
                    self.voronoi.vertices[idx][0] = self.domains['x']['max']
                if y < self.domains['y']['min']:
                    self.voronoi.vertices[idx][1] = self.domains['y']['min']
                if y > self.domains['y']['max']:
                    self.voronoi.vertices[idx][1] = self.domains['y']['max']


    def points_contain_duplicates(self):
        '''
        Return a boolean indicating whether self.points contains duplicates
        '''
        vals, count = np.unique(self.points, return_counts=True)
        return np.any(vals[count > 1])


    def find_centroid(self, vertices):
        '''
        Find the centroid of a Voroni region described by `vertices`,
        and return a np array with the x and y coords of that centroid.
        The equation for the method used here to find the centroid of a
        2D polygon is given here: https://en.wikipedia.org/wiki/
          Centroid#Of_a_polygon
        @params: np.array `vertices` a numpy array with shape n,2
        @returns np.array a numpy array that defines the x, y coords
          of the centroid described by `vertices`
        '''
        area = 0
        centroid_x = 0
        centroid_y = 0
        for i in range(len(vertices)-1):
            step = (vertices[i  , 0] * vertices[i+1, 1]) - \
                 (vertices[i+1, 0] * vertices[i  , 1])
            area += step
            centroid_x += (vertices[i, 0] + vertices[i+1, 0]) * step
            centroid_y += (vertices[i, 1] + vertices[i+1, 1]) * step
        area /= 2
        # prevent division by zero - equation linked above
        if area == 0: area += 0.0000001
        centroid_x = (1.0/(6.0*area)) * centroid_x
        centroid_y = (1.0/(6.0*area)) * centroid_y
        # prevent centroids from escaping bounding box
        if self.constrain:
            if centroid_x < self.domains['x']['min']: centroid_x = self.domains['x']['min']
            if centroid_x > self.domains['x']['max']: centroid_x = self.domains['x']['max']
            if centroid_y < self.domains['y']['min']: centroid_y = self.domains['y']['min']
            if centroid_y > self.domains['y']['max']: centroid_y = self.domains['y']['max']
        return np.array([centroid_x, centroid_y])


    def relax(self):
        '''
        Moves each point to the centroid of its cell in the voronoi
        map to "relax" the points (i.e. jitter the points so as
        to spread them out within the space).
        Ignores vertices pushed to infinity by Voronoi algorithm
        Centroids that fall beyond bounding box are pushed to limits of bounding box
        Jan_2023 (Angueyra):
            added clipping using border as polygon.
            points at infinity are still a problem with clipping
        '''
        centroids = []
        clip = PolygonClipper() # added by Angueyra
        for idx in self.voronoi.point_region:
            # the region is a series of indices into self.voronoi.vertices
            # remove point at infinity, designated by index -1
            region = [i for i in self.voronoi.regions[idx] if i != -1]
            # enclose the polygon
            region = region + [region[0]]
            # get the vertices for this region
            verts = self.voronoi.vertices[region]
            verts = clip(verts,self.get_bb_pointsClockwise()) # added by Angueyra
            # find the centroid of those vertices
            centroids.append(self.find_centroid(verts))
        self.points = np.array(centroids)
        self.constrain_points()
        self.jitter_points()
        self.build_voronoi()

    def relaxAlt(self):
        '''
        Moves each point to the centroid of its cell in the voronoi
        map to "relax" the points (i.e. jitter the points so as
        to spread them out within the space).
        Jan_2023 (Angueyra):
            Using voronoi_finite_polygons_2d to deal with infinite regions
            Subsequently using border for clipping before relaxation
        '''
        centroids = []

        border = self.get_bb_pointsClockwise()
        regions, vertices = voronoi_finite_polygons_2d(self.voronoi)
        clip = PolygonClipper()
        # colorize
        for region in regions:
            verts = vertices[region]
            verts = clip(verts,border)
            # find the centroid of those vertices
            centroids.append(self.find_centroid(verts))
        self.points = np.array(centroids)
        self.constrain_points()
        self.jitter_points()
        self.build_voronoi()

    def get_points(self):
        '''
        Return the input points in the new projected positions
        @returns np.array a numpy array that contains the same number
          of observations in the input points, in identical order
        '''
        return self.points



class PolygonClipper:
# POINTS NEED TO BE PRESENTED CLOCKWISE OR ELSE THIS WONT WORK
    def __init__(self,warn_if_empty=True):
        self.warn_if_empty = warn_if_empty

    def is_inside(self,p1,p2,q):
        R = (p2[0] - p1[0]) * (q[1] - p1[1]) - (p2[1] - p1[1]) * (q[0] - p1[0])
        if R <= 0:
            return True
        else:
            return False

    def compute_intersection(self,p1,p2,p3,p4):

        """
        given points p1 and p2 on line L1, compute the equation of L1 in the
        format of y = m1 * x + b1. Also, given points p3 and p4 on line L2,
        compute the equation of L2 in the format of y = m2 * x + b2.

        To compute the point of intersection of the two lines, equate
        the two line equations together

        m1 * x + b1 = m2 * x + b2

        and solve for x. Once x is obtained, substitute it into one of the
        equations to obtain the value of y.

        if one of the lines is vertical, then the x-coordinate of the point of
        intersection will be the x-coordinate of the vertical line. Note that
        there is no need to check if both lines are vertical (parallel), since
        this function is only called if we know that the lines intersect.
        """

        # if first line is vertical
        if p2[0] - p1[0] == 0:
            x = p1[0]

            # slope and intercept of second line
            m2 = (p4[1] - p3[1]) / (p4[0] - p3[0])
            b2 = p3[1] - m2 * p3[0]

            # y-coordinate of intersection
            y = m2 * x + b2

        # if second line is vertical
        elif p4[0] - p3[0] == 0:
            x = p3[0]

            # slope and intercept of first line
            m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
            b1 = p1[1] - m1 * p1[0]

            # y-coordinate of intersection
            y = m1 * x + b1

        # if neither line is vertical
        else:
            m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
            b1 = p1[1] - m1 * p1[0]

            # slope and intercept of second line
            m2 = (p4[1] - p3[1]) / (p4[0] - p3[0])
            b2 = p3[1] - m2 * p3[0]

            # x-coordinate of intersection
            x = (b2 - b1) / (m1 - m2)

            # y-coordinate of intersection
            y = m1 * x + b1

        intersection = (x,y)

        return intersection

    def clip(self,subject_polygon,clipping_polygon):

        final_polygon = subject_polygon.copy()

        for i in range(len(clipping_polygon)):

            # stores the vertices of the next iteration of the clipping procedure
            next_polygon = final_polygon.copy()

            # stores the vertices of the final clipped polygon
            final_polygon = []

            # these two vertices define a line segment (edge) in the clipping
            # polygon. It is assumed that indices wrap around, such that if
            # i = 1, then i - 1 = K.
            c_edge_start = clipping_polygon[i - 1]
            c_edge_end = clipping_polygon[i]

            for j in range(len(next_polygon)):

                # these two vertices define a line segment (edge) in the subject
                # polygon
                s_edge_start = next_polygon[j - 1]
                s_edge_end = next_polygon[j]

                if self.is_inside(c_edge_start,c_edge_end,s_edge_end):
                    if not self.is_inside(c_edge_start,c_edge_end,s_edge_start):
                        intersection = self.compute_intersection(s_edge_start,s_edge_end,c_edge_start,c_edge_end)
                        final_polygon.append(intersection)
                    final_polygon.append(tuple(s_edge_end))
                elif self.is_inside(c_edge_start,c_edge_end,s_edge_start):
                    intersection = self.compute_intersection(s_edge_start,s_edge_end,c_edge_start,c_edge_end)
                    final_polygon.append(intersection)

        return np.asarray(final_polygon)

    def __call__(self,A,B):
        clipped_polygon = self.clip(A,B)
        if len(clipped_polygon) == 0 and self.warn_if_empty:
            warnings.warn("No intersections found. Are you sure your \
                          polygon coordinates are in clockwise order?")

        return clipped_polygon

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


class NNDcalculator():
    '''
    Calculates NND (nearest-neighbour distances) and idealized NND (after running Lloyd's algorithm) using napari points or labels layer
    Created by Angueyra (Jan 2023)
    '''

    def __init__(self, dPath = '', celltype = 'U', datatype = 'points', bg = 'N_mip', nIter = 20, saveLloyd = False, overwriteFlag = False, **kwargs):
        '''
        descriptionHere
        '''
        self.dPath = dPath
        self.celltype = celltype
        self.datatype = datatype
        if (datatype=='labels'):
            self.fPath = self.celltype + '_seg_curated.tiff'
        elif (datatype=='points'):
            if (self.celltype=='U'):
                self.fPath = 'U_missing.csv'
            else:
                self.fPath = self.celltype + '_points.csv'
        self.ifPath = self.celltype + '_lloyd.npy'
        self.bg = bg + '.tiff'
        self.p2um = 63.64/1024
        self.kN = 1
        self.nIter = nIter
        self.Centroids = self.getCentroids();
        if (np.size(self.Centroids)>0): # skip any calculations if filepath is empty
            self.vor = scipy.spatial.Voronoi(self.Centroids)
            self.kdt = scipy.spatial.KDTree(self.Centroids)
            self.nnd = self.calcNND()
            self.iCentroids = self.calcLloyd(nIter = self.nIter, overwriteFlag = overwriteFlag)
            self.ivor = scipy.spatial.Voronoi(self.iCentroids)
            self.ikdt = scipy.spatial.KDTree(self.iCentroids)
            self.innd = self.calcNNDi()
            if saveLloyd:
                self.saveLloyd()

    def getCentroids(self):
        # return empty array if file does not exist
        if (os.path.isfile(self.dPath + self.fPath)==False):
            warnings.warn('{0} does not exist'.format(self.dPath + self.fPath))
            Centroids = [];
        else:
            if (self.datatype=='points'): # assumes input is structured like a napari points layer
                df =  pd.read_csv(self.dPath + self.fPath)
                df = df.drop(columns = 'index')
                Centroids = df.to_numpy()
                Centroids = Centroids[:,[1,0]]
            elif(self.datatype=='labels'): # input is an image file akin to napari labels layer
                labelData = imread (self.dPath + self.fPath)
                Centroids = calculateCentroidsFromLabels(labelData)
        return Centroids

    def calcNND(self):
        nnd = self.findNND(self.Centroids, self.kdt, self.kN)
        return nnd

    def calcNNDi(self):
        nndi = self.findNND(self.iCentroids, self.ikdt, self.kN)
        return nndi

    @staticmethod
    def findNND(centroids, kdt, kN):
        # initialize emtpy arrays
        nnd = np.empty((len(centroids),kN))
        nnd[:]=np.NaN
        # find k nearest neighbor distances
        for n in range(0,len(centroids)):
            temp = kdt.query(centroids[n],kN+1)
            for k in range (0,kN):
                nnd[n,k] = temp[0][k+1]
        return nnd


    def calcLloyd(self,nIter=20, overwriteFlag = False):
        if ((os.path.isfile(self.dPath + self.ifPath)) & (overwriteFlag == False)):
            print('Loaded saved Lloyd relaxation')
            iCentroids = np.load(self.dPath + self.ifPath, allow_pickle=True)
        else:
            #start from original locations
            field = Field(self.Centroids, buffer = 20)
            # # #start from center of the field
            # rng = np.random.default_rng()
            # rndCentroids = rng.integers(500,525, size=np.shape(self.Centroids))
            # field = Field(rndCentroids, buffer = 20)
            with tqdm(total=nIter, file=sys.stdout) as progBar:
            # run Lloyd relaxation on the field of points
                for i in np.arange(0,nIter):
                    field.relax()
                    progBar.update(1)
            # get the resulting point positions
            iCentroids = field.get_points()
        return iCentroids

    def saveLloyd(self):
        np.save(self.dPath + self.ifPath, self.iCentroids)
        np.save(self.dPath + self.celltype + '_nnd.npy', self.nnd)
        np.save(self.dPath + self.celltype + '_innd.npy', self.innd)
        print('Saved results for {0}'.format(self.dPath))

    def plotAll(self):
        fH, axH = plt.subplots(2, 2, figsize= [24,24], gridspec_kw={'height_ratios': [5, 1]})
        axH[0,1].sharey(axH[0,0])
        axH[1,1].sharey(axH[1,0])
        self.plotVor(fH, axH[0,0]);
        self.plotLloyd(fH, axH[0,1]);
        self.plotNDD(fH,axH[1,0]);
        self.plotNDDLloyd(fH,axH[1,1]);
        fH.tight_layout()
        return fH, axH

    def plotVor(self, fH, axH, plotRegions = False):
        '''
        plotter for N_mip, U locations and voronoi partition
        '''

        if (os.path.isfile(self.dPath + self.bg)):
            pH = axH.imshow(imread(self.dPath + self.bg), cmap='bone')

        pH = axH.scatter(self.Centroids[:,0],self.Centroids[:,1], color=zfC[self.celltype], zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = 0.99);
        pH = scipy.spatial.voronoi_plot_2d(self.vor, ax = axH, show_points = False, show_vertices=False, line_colors=zfC[self.celltype], line_width=2, line_alpha=0.99, point_size=2)

        if plotRegions:
            regions, vertices = voronoi_finite_polygons_2d(self.vor)
            # colorize
            for region in regions:
                polygon = vertices[region]
                axH.fill(*zip(*polygon), alpha=0.2)

        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_xlim([-100,1124])
        axH.set_ylim([-100,1124])

        def p2umfx(x, pos): # formatter function to convert pixels to microns
            s = '{:.2f}'.format(np.round(x * self.p2um,decimals=2))
            return s

        axH.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        axH.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        return fH, axH

    def plotLloyd(self, fH, axH):
        '''
        plotter for N_mip, U locations and voronoi partition after Lloyd's relaxation
        '''

        if (os.path.isfile(self.dPath + self.bg)):
            pH = axH.imshow(imread(self.dPath + self.bg), cmap='bone')

        pH = axH.scatter(self.Centroids[:,0],self.Centroids[:,1], color=zfC[self.celltype], zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = .2);
        pH = scipy.spatial.voronoi_plot_2d(self.vor, ax = axH, show_points = False, show_vertices=False, line_colors=zfC[self.celltype], line_width=2, line_alpha=0.2, point_size=2)

        pH = axH.scatter(self.iCentroids[:,0],self.iCentroids[:,1], color='#FFFFFF', zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = .99);
        pH = scipy.spatial.voronoi_plot_2d(self.ivor, ax = axH, show_points = False, show_vertices=False, line_colors='#FFFFFF', line_width=2, line_alpha=0.99, point_size=2)

        # regions, vertices = voronoi_finite_polygons_2d(self.vor)
        # # colorize
        # for region in regions:
        #     polygon = vertices[region]
        #     plt.fill(*zip(*polygon), alpha=0.2)

        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_xlim([-100,1124])
        axH.set_ylim([-100,1124])

        def p2umfx(x, pos): # formatter function to convert pixels to microns
            s = '{:.2f}'.format(np.round(x * self.p2um,decimals=2))
            return s

        axH.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        axH.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        return fH, axH

    def plotNDD(self, fH, axH):
        '''
        plotter for NDD distribution
        '''
        for k in range (0,self.kN):
            pH = axH.hist(self.nnd[:,k]*self.p2um, bins=20, color=zfC[self.celltype], edgecolor=baseColor, linewidth=0.5, alpha = .8, histtype='stepfilled');

        axH.set_xlim([0,np.max((self.nnd*self.p2um)*1.2)])
        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('nnd ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('counts', fontproperties=fontLabels)

        nndmedian = np.median(self.nnd[:,k])
        pH = axH.axvline(nndmedian*self.p2um, color = lighten_color(zfC[self.celltype],0.5), linestyle='--')

        # summarize using median
        print ('\n number of UV cones = {0}'.format(len(self.Centroids)))
        for k in range (0,self.kN):
            print ('Median for {0}-NDD (Raw):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,nndmedian,nndmedian*self.p2um))
            # print ('Median for {0}-NDD UV-cones (Lloyd):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,inndmedian,inndmedian*self.p2um))

    def plotNDDLloyd(self, fH, axH):
        '''
        plotter for NDD distribution after Lloyd's relaxation
        '''
        for k in range (0,self.kN):
            pH = axH.hist(self.nnd[:,k]*self.p2um, bins=20, color=zfC[self.celltype], edgecolor=baseColor, linewidth=0.5, alpha = .8, histtype='stepfilled');

            pH = axH.hist(self.innd[:,k]*self.p2um, bins=20, color='#FFFFFF', edgecolor=baseColor, linewidth=0.5, alpha = .8, histtype='stepfilled');

        axH.set_xlim([0,np.max((self.nnd*self.p2um)*1.2)])
        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('nnd ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('counts', fontproperties=fontLabels)

        nndmedian = np.median(self.nnd[:,k])
        pH = axH.axvline(nndmedian*self.p2um, color = lighten_color(zfC[self.celltype],0.5), linestyle='--')


        inndmedian = np.median(self.innd[:,k])
        pH = axH.axvline(inndmedian*self.p2um, color = lighten_color('#ffffff',0.5), linestyle=':')

        # summarize using median
        for k in range (0,self.kN):
            print ('Median for {0}-NDD (Lloyd):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,inndmedian,inndmedian*self.p2um))


def calculateCentroidsFromLabels(labelsData):
    props = regionprops_table(labelsData,properties=('label','centroid'))
    centroidData = np.column_stack((props['centroid-1'],props['centroid-0']))
    return centroidData

"""
class NNDUcalculator():
    '''
    initial class intended to calculate NND and idealized NDD for UV cones starting with points data
    later decided to generalize class to take more flexible inputs capable of dealing with other photoreceptor files
    Created by Angueyra (Jan 2023)
    '''

    def __init__(self, dPath = '', nIter = 20, saveLloyd = False, overwriteFlag = False, **kwargs):
        '''
        descriptionHere
        '''
        self.dPath = dPath
        self.uPath = 'U_missing.csv'
        self.uiPath = 'U_lloyd.npy'
        self.N_mip = 'N_mip.tiff'
        self.p2um = 63.64/1024
        self.kN = 1
        self.nIter = nIter
        self.U = self.getU();
        self.vorU = scipy.spatial.Voronoi(self.U)
        self.kdtU = scipy.spatial.KDTree(self.U)
        self.nndU = self.calcNDD()
        self.Ui = self.calcLloyd(nIter = self.nIter, overwriteFlag = overwriteFlag)
        self.vorUi = scipy.spatial.Voronoi(self.Ui)
        self.kdtUi = scipy.spatial.KDTree(self.Ui)
        self.nndUi = self.calcNDDi()
        if saveLloyd:
            self.saveLloyd()

    def getU(self):
        if (os.path.isfile(self.dPath + self.uPath)==False):
            raise TypeError('{0} does not exist'.format(self.dPath + self.uPath))
        df =  pd.read_csv(self.dPath + self.uPath)
        df = df.drop(columns = 'index')
        U=df.to_numpy()
        U = U[:,[1,0]]
        return U

    def calcNDD(self):
        # initialize emtpy arrays
        nndU = np.empty((len(self.U),self.kN))
        nndU[:]=np.NaN

        # find k nearest neighbor distances
        for n in range(0,len(self.U)):
            temp = self.kdtU.query(self.U[n],self.kN+1)
            for k in range (0,self.kN):
                nndU[n,k] = temp[0][k+1]
        return nndU

    def calcNDDi(self):
        # initialize emtpy arrays
        nndUi = np.empty((len(self.Ui),self.kN))
        nndUi[:]=np.NaN

        # find k nearest neighbor distances
        for n in range(0,len(self.Ui)):
            temp = self.kdtUi.query(self.Ui[n],self.kN+1)
            for k in range (0,self.kN):
                nndUi[n,k] = temp[0][k+1]
        return nndUi

    def calcLloyd(self,nIter=20, overwriteFlag = False):
        if ((os.path.isfile(self.dPath + self.uiPath)) & (overwriteFlag == False)):
            print('Loaded saved Lloyd relaxation')
            Ui = np.load(self.dPath + self.uiPath, allow_pickle=True)
        else:
            field = Field(self.U, buffer = 20) #start from original locations of UV cones
            with tqdm(total=nIter, file=sys.stdout) as progBar:
            # run Lloyd relaxation on the field of points
                for i in np.arange(0,nIter):
                    field.relax()
                    progBar.update(1)
            # get the resulting point positions
            Ui = field.get_points()
        return Ui

    def saveLloyd(self):
        np.save(self.dPath + self.uiPath, self.Ui)
        np.save(self.dPath + 'nndU.npy', self.nndU)
        np.save(self.dPath + 'nndUi.npy', self.nndUi)
        print('Saved results for {0}'.format(self.dPath))

    def plotAll(self):
        fH, axH = plt.subplots(2, 2, figsize= [24,24], gridspec_kw={'height_ratios': [5, 1]})
        axH[0,1].sharey(axH[0,0])
        axH[1,1].sharey(axH[1,0])
        self.plotVor(fH, axH[0,0]);
        self.plotLloyd(fH, axH[0,1]);
        self.plotNDD(fH,axH[1,0]);
        self.plotNDDLloyd(fH,axH[1,1]);
        fH.tight_layout()
        return fH, axH

    def plotVor(self, fH, axH, plotRegions = False):
        '''
        plotter for N_mip, U locations and voronoi partition
        '''

        if (os.path.isfile(self.dPath + self.N_mip)):
            pH = axH.imshow(imread(self.dPath + self.N_mip), cmap='bone')

        pH = axH.scatter(self.U[:,0],self.U[:,1], color=zfC['U'], zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = 0.99);
        pH = scipy.spatial.voronoi_plot_2d(self.vorU, ax = axH, show_points = False, show_vertices=False, line_colors=zfC['U'], line_width=2, line_alpha=0.99, point_size=2)

        if plotRegions:
            regions, vertices = voronoi_finite_polygons_2d(self.vorU)
            # colorize
            for region in regions:
                polygon = vertices[region]
                axH.fill(*zip(*polygon), alpha=0.2)

        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_xlim([-100,1124])
        axH.set_ylim([-100,1124])

        def p2umfx(x, pos): # formatter function to convert pixels to microns
            s = '{:.2f}'.format(np.round(x * self.p2um,decimals=2))
            return s

        axH.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        axH.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        return fH, axH

    def plotLloyd(self, fH, axH):
        '''
        plotter for N_mip, U locations and voronoi partition after Lloyd's relaxation
        '''

        if (os.path.isfile(self.dPath + self.N_mip)):
            pH = axH.imshow(imread(self.dPath + self.N_mip), cmap='bone')

        pH = axH.scatter(self.U[:,0],self.U[:,1], color=zfC['U'], zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = .2);
        pH = scipy.spatial.voronoi_plot_2d(self.vorU, ax = axH, show_points = False, show_vertices=False, line_colors=zfC['U'], line_width=2, line_alpha=0.2, point_size=2)

        pH = axH.scatter(self.Ui[:,0],self.Ui[:,1], color=zfC['S'], zorder=8, marker = 'o', s=50, edgecolor=baseColor, linewidth=0.5, alpha = .99);
        pH = scipy.spatial.voronoi_plot_2d(self.vorUi, ax = axH, show_points = False, show_vertices=False, line_colors=zfC['S'], line_width=2, line_alpha=0.99, point_size=2)

        # regions, vertices = voronoi_finite_polygons_2d(self.vorU)
        # # colorize
        # for region in regions:
        #     polygon = vertices[region]
        #     plt.fill(*zip(*polygon), alpha=0.2)

        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('d ($\mu$m)', fontproperties=fontLabels)
        axH.set_xlim([-100,1124])
        axH.set_ylim([-100,1124])

        def p2umfx(x, pos): # formatter function to convert pixels to microns
            s = '{:.2f}'.format(np.round(x * self.p2um,decimals=2))
            return s

        axH.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        axH.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(p2umfx))
        return fH, axH

    def plotNDD(self, fH, axH):
        '''
        plotter for NDD distribution
        '''
        for k in range (0,self.kN):
            pH = axH.hist(self.nndU[:,k]*self.p2um, bins=20, color=zfC['U'], edgecolor=baseColor, linewidth=0.5, alpha = .8, histtype='stepfilled');

        axH.set_xlim([0,np.max((self.nndU*self.p2um)*1.2)])
        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('nnd ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('counts', fontproperties=fontLabels)

        nndUmedian = np.median(self.nndU[:,k])
        pH = axH.axvline(nndUmedian*self.p2um, color = '#ffffff80', linestyle='--')

        # summarize using median
        print ('\n number of UV cones = {0}'.format(len(self.U)))
        for k in range (0,self.kN):
            print ('Median for {0}-NDD UV-cones (Raw):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,nndUmedian,nndUmedian*self.p2um))
            # print ('Median for {0}-NDD UV-cones (Lloyd):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,nndUimedian,nndUimedian*self.p2um))

    def plotNDDLloyd(self, fH, axH):
        '''
        plotter for NDD distribution after Lloyd's relaxation
        '''
        for k in range (0,self.kN):
            pH = axH.hist(self.nndU[:,k]*self.p2um, bins=20, color=zfC['U'], edgecolor=baseColor, linewidth=0.5, alpha = .6, histtype='stepfilled');

            pH = axH.hist(self.nndUi[:,k]*self.p2um, bins=20, color=zfC['S'], edgecolor=baseColor, linewidth=0.5, alpha = .8, histtype='stepfilled');

        axH.set_xlim([0,np.max((self.nndU*self.p2um)*1.2)])
        fontLabels = formatFigureMain(fH, axH, pH);
        axH.set_xlabel('nnd ($\mu$m)', fontproperties=fontLabels)
        axH.set_ylabel('counts', fontproperties=fontLabels)

        nndUmedian = np.median(self.nndU[:,k])
        pH = axH.axvline(nndUmedian*self.p2um, color = '#ffffff80', linestyle='--')

        nndUimedian = np.median(self.nndUi[:,k])
        pH = axH.axvline(nndUimedian*self.p2um, color = '#00ffff80', linestyle=':')

        # summarize using median
        for k in range (0,self.kN):
            print ('Median for {0}-NDD UV-cones (Lloyd):\t {1:.2f} pixels\t {2:.2f} um'.format(k+1,nndUimedian,nndUimedian*self.p2um))
"""
