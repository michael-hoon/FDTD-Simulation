import numpy as np
import time
import matplotlib.pyplot as plt

#Only changes to be made to yapic.py is to copy 'our stuff' in, initialise epsilon and mirror object finally add * epsilon behind the *dus
#Initialise tgt with E and B dicts

#our stuff
class Epsil2D():
    def __init__(self, sh, bgep = 1):
        #Takes in shape(dimens) of grid and background epsil value (air)
        self._objects = {}
        self.bgep = 1
        self.epsil = np.full(sh, self.bgep, dtype = float) 

    @property
    def objects(self):
        return self._objects.keys()

    def add_object(self, object):
        if isinstance(object, Object):
            self._objects[object.name] = object
            self._adjust_epsil(object)
        else:
            raise TypeError("Only Objects can be added")

    def remove_object(self, object):
        if object.name in self._objects:
            self._objects.pop(object.name)
            self._adjust_epsil(object, revert = True)
        else:
            print("Object does not exist")

    def get_object(self, name):
        if name in self._objects:
            return self._objects[name]
        else:
            print("Object does not exist")

    def _adjust_epsil(self, object, revert = False):
        #To add feature to detect and handle overlap, objects may need attribute to allow/disallow overlap?
        if revert:
            self.epsil[object.pos] = self.bgep
        else:
            self.epsil[object.pos] = object.ep
        
    def show(self):
        #Visualise objects
        #Extra: option to display edgelengths/points using the object's methods???
        plt.imshow(self.epsil, interpolation = None)
        plt.colorbar(location = 'right')
        plt.show()


class Object():
    def __init__(self, name, ep, sh, exceed = False):
        #initialise with epsil value for object and a grid with same shape as background
        #Child classes will override, calls super().__init__, store additional arguments and call _set_pos with those args
        #Additional args will be dimensions of object with values btwn 0 and 1
        #Object created with dims/size relative to grid
        self.ep = ep
        self.pos = np.zeros(sh, dtype = bool)
        self.exceed = exceed
        self.name = name

    def _set_pos(self):
        #To be called in __init__ in child classes
        #If self.exceed == True, allow for making object that exceeds grid but cuts off
            #Use is_within_grid to adjust algorithm for cutting off
        #Else raise error if object exceeds grid, use is_within_grid
        pass

    @property
    def edges(self):
        #To determine edges
        #Return list of tuples of (min, max) edges in each dim to ensure is within range
        #values of integars
        pass

    @property
    def edgepoints(self):
        #might not need
        #returns list of min/max points in each dim to ensure is within range
        #Values are array indices (coords)
        pass

    @property
    def is_within_grid(self):
        #uses dims of objects to calc edge points to check if within grid
        #Returns true if within grid
        pass

    def show(self):
        vis = np.zeros(self.pos.shape)
        vis[self.pos] = 1
        plt.imshow(vis, interpolation = None)
        plt.colorbar(orientation = 'vertical')
        plt.show()

class Diag2D(Object):
    #Just create a diagonal object spanning entire grid, only 45 degrees, purely created for the mirror 
    #Is there a better way to implement using in-built numpy diag stuff instead of just iterating? I cant find a suitable func (Darren)
    #required dims: center (tuple/list of coords normalised 0~1)
    #               width (normalised width (float) ranging from 0~1)
    #               topleft if true span from top left to bot right, if false span from top right to bot left, default True
    def __init__(self, name, ep, sh, width, center, exceed = False, topleft = True):
        assert len(sh) == 2, f"Diag2D only valid for 2D"
        assert len(sh) == len(center), f"wrong dimensions, grid has shape {sh} while center is {center}"
        assert all([(i >= 0 and i <= 1 for i in center)]), "Center not within grid"
        super().__init__(name, ep, sh, exceed)
        self.center = list(map(lambda size, scale: int(size * scale), sh, center))
        self.width = width * np.sqrt(sum([x ** 2 for x in sh]))
        self.topleft = topleft
        self._set_pos()

    def _set_pos(self):
        if self.exceed:
            pass
        else:
            assert self.is_within_grid, "Object exceeds grid"
            it = np.nditer(self.pos, flags=['multi_index'])
            for i in it:
                x,y = it.multi_index
                if self.is_within_edges(x, y):
                    self.pos[x, y] = True
                    
    @property
    def intercepts(self):
        #replaces edges or edgepoints
        #Returns (minc, maxc) in terms of y = mx + c for the upper and lower bounds of the diagonal
        #m = 1/-1 automatically because only 45 deg depending on if topleft
        margin = self.width // np.sqrt(2)
        x, y = self.center[0], self.center[1]
        minx, miny = x - margin, y - margin
        maxx, maxy = x + margin, y + margin
        if self.topleft:
            return (miny + minx, maxy + maxx)
        else:
            return (miny - maxx, maxy - minx)

    def is_within_edges(self, x, y):
        minc, maxc = self.intercepts
        if self.topleft:
            return y >= minc - x and y <= maxc - x
        else:
            return y >= x + minc and y <= x + maxc

    @property
    def is_within_grid(self):
        minc, maxc = self.intercepts
        maxx, maxy = self.pos.shape
        if self.topleft:
            return maxy > maxc - maxx and 0 < minc
        else:
            return maxy > maxc and 0 < maxx + minc

class CenterRectangle(Object):
    #required dims: dims (tuple/list of normalised lengths, e.g. for 200x50 rect in 500x500 grid, len = (0.4, 0.1))
    #               center (tuple/list of coords also normalised 0~1)
    #center and dims will round down to closest grid index, object may end up smaller than expected if input is not precise
    def __init__(self, name, ep, sh, dims, center, exceed = False):
        assert len(sh) == len(center) and len(sh) == len(dims), f"wrong dimensions, grid has shape {sh} while dims is {dims} and center is {center}"
        assert all([(i >= 0 and i <= 1 for i in center)]), "Center not within grid"
        super().__init__(name, ep, sh, exceed)
        self.center = list(map(lambda size, scale: int(size * scale), sh, center))
        self.dims = list(map(lambda size, scale: size * scale, sh, dims)) #note may not be int, convert when calculating edges/edgepoints
        self._set_pos() 

    def _set_pos(self):
        if self.exceed:
            pass
        else:
            assert self.is_within_grid, "Object exceeds grid"
            if len(self.edges) == 2:
                xmin, xmax = self.edges[0][0], self.edges[0][1]
                ymin, ymax = self.edges[1][0], self.edges[1][1]
                self.pos[xmin:xmax+1, ymin:ymax+1] = True
            elif len(self.edges) == 3:
                xmin, xmax = self.edges[0][0], self.edges[0][1]
                ymin, ymax = self.edges[1][0], self.edges[1][1]
                zmin, zmax = self.edges[2][0], self.edges[2][1]
                self.pos[xmin:xmax+1, ymin:ymax+1, zmin:zmax+1] = True

    @property
    def edges(self):
        #Because rectangle, returns list of tuples in [(min, max), ...]
        return list(map(lambda center, dim: (int(center - dim//2), int(center + dim//2)), self.center, self.dims))

    @property    
    def edgepoints(self):
        assert len(self.edges) == 2, "Only implemented for 2D so far"
        if len(self.edges) == 2:
            xmin, xmax = self.edges[0][0], self.edges[0][1]
            ymin, ymax = self.edges[1][0], self.edges[1][1]
            return {"topleft": (xmin, ymax),
                    "topright": (xmax, ymax),
                    "botleft": (xmin, ymin),
                    "botright": (xmax, ymin)}
        elif len(self.edges) == 3:
            pass
            
    @property
    def is_within_grid(self):
        shape = self.pos.shape
        for i in range(len(shape)):
            if self.edges[i][0] < 0 or self.edges[i][1] > shape[i]:
                return False
        return True

class FancyRectangle2D(Object):
    #required dims: dims (tuple/list of normalised lengths, e.g. for 200x50 rect in 500x500 grid, len = (0.4, 0.1))
    #               center (tuple/list of coords also normalised 0~1)
    #               angle of rotation in radians: anticlockwise rotation from x-axis, where dims is (x,y)
    #center and dims will round down to closest grid index, object may end up smaller than expected if input is not precise
    def __init__(self, name, ep, sh, dims, center, angle, exceed = False):
        assert len(sh) == len(center) and len(sh) == len(dims), f"wrong dimensions, grid has shape {sh} while dims is {dims} and center is {center}"
        assert all([(i >= 0 and i <= 1 for i in center)]), "Center not within grid"
        super().__init__(name, ep, sh, exceed)
        self.center = list(map(lambda size, scale: int(size * scale), sh, center))
        self.dims = list(map(lambda size, scale: size * scale, sh, dims)) #note may not be int, convert when calculating edges/edgepoints
        self.angle = angle
        self._set_pos()

    def _set_pos(self):
        pass

    @property
    def edges(self):
        pass

    @property
    def edgepoints(self):
        pass

    @property
    def is_within_grid(self):
        pass 

#stuff copied from yapicc to setup
dim = 2
xlim = [ -25.0e-4,  25.0e-4,
         -25.0e-4,  25.0e-4]
xspacing = [9,9]
dxs = np.array([
    (xlim[2*i+1] - xlim[2*i])/xspacing[i] for i in range(dim)
    ])
posx = np.mgrid[
    xlim[0]:xlim[1]:(xspacing[0]+1)*1j,
    xlim[2]:xlim[3]:(xspacing[1]+1)*1j]
posxh = np.array([
    ix + idx/2.0 for ix,idx in zip(posx,dxs) ])
sh = posx.shape[1:]

#Tests
cenrec1 = CenterRectangle('cenrec1', 1, (100,100), (0.5, 0.5), (0.5, 0.5))
assert cenrec1.edges == [(25, 75),(25, 75)], f'{cenrec1.edges}'
assert cenrec1.edgepoints == {"topleft": (25, 75),
                            "topright": (75, 75),
                            "botleft": (25, 25),
                            "botright": (75, 25)}, f'{cenrec1.edgepoints}'
assert cenrec1.is_within_grid
# cenrec1.show()

cenrec2 = CenterRectangle('cenrec2', 1,  (100,100), (1, 0.5), (0.255, 0.755), exceed = True)
assert cenrec2.edges == [(-25, 75),(50, 100)], f'{cenrec2.edges}'
assert not cenrec2.is_within_grid

# cenrec3 = CenterRectangle('cenrec3', 1, (11,11), (0.5, 0.5), (0.5, 0.5))
# print(cenrec3.pos)
# print(cenrec3.edges)
# cenrec3.show()

diag1 = Diag2D('diag1', 0.5, (10,10), 0.05, (0.5, 0.5), topleft = False)
# print(diag1.intercepts)
# print(diag1.center)
assert diag1.is_within_grid
# diag1.show()

epsil = Epsil2D((10,10), 1)
epsil.add_object(diag1)
# epsil.show()

# FIX diag2d math- the width
# Make wall objects - init epsil with wall objects
# Need to handle overlapping objects, overlap = True/False, true-true, false-false, true-false