import numpy as np
import time
import matplotlib.pyplot as plt

#Only changes to be made to yapic.py is to copy 'our stuff' in, initialise epsilon and mirror object finally add * epsilon behind the *dus
#Initialise tgt with E and B dicts

#TO BE IMPLEMENTED
# feature to handle overlapping objects, overlap = True/False, true-true, false-false, true-false
# feature to handle objects that exceed the grid, an alt set_pos that cuts off
# Fancy rectangle

#our stuff
class Epsil2D():
    """
    A class to represent the epsilon value.

    '''

    Attributes 
    ----------
        sh : array
            Shape (dimensions) of the grid 
        bgep : int
            Background epsilon value (relative permittivity) of the surrounding material (air = 1)
        boundary : bool
            Determines if boundaries of the simulation are set, default is True

    Methods
    -------
        objects()
            Returns the key values of the object ##

        add_object(object)
            Initiates an object 

        remove_object(object)
            Removes an existing object 

        get_object(name)
            Returns the named object 

        _adjust_epsil(object, revert = False)
            Adjusts the epsilon value to detect and handle overlaps ## Not completed yet

    show()
        Visualise objects using the matplotlib plot function 

    """
    def __init__(self, sh, bgep = 1, boundary = True):
        """
        Constructs all the necessary attributes required for the Epsil2D object.
        If argument 'bgep' isn't passed in, the default (air = 1) is used.

        Parameters
        ----------
            sh : array
                Shape (dimensions) of the grid 
            bgep : int
                Background epsilon value (relative permittivity) of the surrounding material.
            boundary : bool
                Default value is True, will fix all edge values at 1 regardless of the object added.
                If set to false, will initially set at bgep and be affected by objects added.

        """
        
        self._objects = {}
        self.bgep = 1
        self.boundary = boundary
        self.epsil = np.full(sh, self.bgep, dtype = float)
        if self.boundary:
            self.set_bounds() 

    @property
    def objects(self):
        """
        nil 
        
        """
        return self._objects.keys()

    def add_object(self, object):
        """
        Instantiates an object 
        
        Parameters
        ----------
            object : Object
        
        Raises
        ------
            TypeError
                If the object raised is not of the Object type.

        """
        if isinstance(object, Object):
            self._objects[object.name] = object
            self._adjust_epsil(object)
        else:
            raise TypeError("Only Objects can be added")

    def remove_object(self, object):
        """
        Removes an existing object.
        
        """
        if object.name in self._objects:
            self._objects.pop(object.name)
            self._adjust_epsil(object, revert = True)
        else:
            print("Object does not exist")

    def get_object(self, name):
        """
        Returns the name of the object. 

        """
        if name in self._objects:
            return self._objects[name]
        else:
            print("Object does not exist")

    def _adjust_epsil(self, object, revert = False):

        """
        Adjusts the epsilon value to detect and handle overlaps, ## Not complete.

        Parameters
        ----------
            object: Object
                Initiated object.
            revert : bool
                Changes the current epsilon value to the initial epsilon value (default is False).

        """

        #To add feature to detect and handle overlap, objects may need attribute to allow/disallow overlap?
        if revert:
            self.epsil[object.pos] = self.bgep
        else:
            self.epsil[object.pos] = object.ep
        if self.boundary:
            self.set_bounds()

    def set_bounds(self, ep = 1):
        """
        Sets the boundaries of the simulation. This can also be called after all objects are added, to set at different ep value if required. 

        Parameters
        ----------
            ep : int
                Value for the boundaries of the simulation

        """

        self.epsil[0, :] = ep
        self.epsil[-1, :] = ep
        self.epsil[:, 0] = ep
        self.epsil[:, -1] = ep
        
    def show(self):
        """
        Visualise the objects via the matplotlib imshow(), colorbar() functions, and outputs the graph using the show() function. 

        """

        #Extra: option to display edgelengths/points using the object's methods???
        plt.imshow(self.epsil, interpolation = None)
        plt.colorbar(location = 'right')
        plt.show()

class Object():
    """
    A Class to represent the different objects initialised.
    
    '''

    Attributes
    ----------
        name: str
            Name of the object to be initialised

        ep: int
            Epsilon value for the object 

        sh : array
            Shape (dimensions) of the grid, with the same shape as background

        exceed : bool


    Methods
    -------
        set_pos():
            To be called in .__init__ of the child classes.  
        
        edges():
            Determines the edges of the object.
            Returns a list of tuples of (min, max) edges in each dimension to ensure that it is within the range of values of integers.

        edgepoints():
            Returns a list of min/max points in each dimension to ensure it is within the range. 
            Values are array indices (coordinates).

        is_within_grid():
            Checks if the edge points of the object are within the grid dimensions.

        show():
            Uses the matplotlib imshow(), colorbar() functions to visualise the object.

    """

    def __init__(self, name, ep, sh, exceed = False):
        """
        Constructs all the necessary attributes for the Object object. Child classes will override the attributes, calls super().__init__, and store the additional arguments. 
        _set_pos is then called with those arguments. The additional arguments will be the dimensions of the object, with values between 0 and 1. 
        Objects are created with the dimensions relative to the grid.
        
        Parameters
        ----------
            name: str
                Name of the object to be initialised

            ep: int
                Epsilon value for the object 

            sh : array
                Shape (dimensions) of the grid, with the same shape as background
        
            exceed : bool


        """

        self.ep = ep
        self.pos = np.zeros(sh, dtype = bool)
        self.exceed = exceed
        self.name = name

    def _set_pos(self):
        """
        To be called in __init__ in the child classes. If self.exceed == True, allows for making an object that exceeds the grid, but gets cut off. 
        Use is_within_grid to adjust algorithm for cutting off.

        Raises
        ------


        """
        #To be called in __init__ in child classes
        #If self.exceed == True, allow for making object that exceeds grid but cuts off
            #Use is_within_grid to adjust algorithm for cutting off
        #Else raise error if object exceeds grid, use is_within_grid
        pass

    @property
    def edges(self):
        """
        Determines the edges of the object.
        Returns a list of tuples of (min, max) edges in each dimension to ensure that it is within the range of values of integers.
        
        """
        pass

    @property
    def edgepoints(self):
        """
        Returns a list of min/max points in each dimension to ensure it is within the range. 
        Values are array indices (coordinates) 
        
        """
        #might not need

        pass

    @property
    def is_within_grid(self):
        """
        Checks if the edge points of the object are within the grid dimensions

        Returns
        -------
            True: if edge points are within the grid

        """
        #uses dims of objects to calc edge points to check if within grid

        pass

    def show(self):
        """
        Uses the matplotlib imshow(), colorbar() functions to visualise the object.

        """
        vis = np.zeros(self.pos.shape)
        vis[self.pos] = 1
        plt.imshow(vis, interpolation = None)
        plt.colorbar(orientation = 'vertical')
        plt.show()

class Diag2D(Object):
    """
    A class that represents the diagonal object initialised in the simulation. 
    Created purely to test out mirror object in grid, initialised with 45 degree angle across, spanning the grid.
    Inherits the name, ep, sh, and exceed parameters from the Object class.

    '''

    Attributes
    ----------
        name : str
            Name of the diagonal object to be initialised.

        ep: int
            Epsilon value for the object 

        sh : int
            Shape (dimensions) of the grid, with the same shape as background

        width : float
            Normalised width of the diagonal ranging from 0 - 1
        
        center : tuple
            Normalised tuple of coordinates ranging from 0 - 1

        exceed : bool

        topleft : bool
            If True, spans the diagonal from the top left of the boundary to the bottom right. Vice versa. 
            Defaults to true if argument is not specified. 

    
    Methods
    -------
        _set_pos():
            Checks if the positions of the diagonal exceeds the grid boundaries. 
            If exceed is True, 
        
        intercepts():
            Replaces edges or edgepoints method, and returns the intercepts for the upper and lower bounds of the diagonal. 



    """

    #Is there a better way to implement using in-built numpy diag stuff instead of just iterating? I cant find a suitable func (Darren)

    def __init__(self, name, ep, sh, width, center, exceed = False, topleft = True):
        """
        Constructs the necessary attributes for the Diag2D class. 

        Parameters
        ----------
            name : str
                Name of the diagonal object to be initialised.

            ep: int
                Epsilon value for the object 

            sh : int
                Shape (dimensions) of the grid, with the same shape as background. Only valid for 2 Dimensions.

            width : float
                Normalised width of the diagonal ranging from 0 - 1
            
            center : tuple
                Normalised tuple of coordinates ranging from 0 - 1

            exceed : bool

            topleft : bool
                If True, spans the diagonal from the top left of the boundary to the bottom right, vice versa.
                (Default is True) 
        """
        
        assert len(sh) == 2, f"Diag2D only valid for 2D"
        assert len(sh) == len(center), f"wrong dimensions, grid has shape {sh} while center is {center}"
        assert all([(i >= 0 and i <= 1 for i in center)]), "Center not within grid"
        super().__init__(name, ep, sh, exceed)
        self.center = list(map(lambda size, scale: int(size * scale), sh, center))
        self.width = width * np.sqrt(sum([x ** 2 for x in sh])) / 2
        self.topleft = topleft
        self._set_pos()

    def _set_pos(self):
        """
        
        """
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
        """
        Replace edges or edgepoints methods.
        Returns (minc, maxc) in terms of y = mx + c for the upper and lower bounds of the diagonal.
        m = 1 or -1, depending on the value of topleft, as the diagonal is set to 45 degrees.

        Returns
        -------
        (minc, maxc)
            tuple: 
                (miny + minx, maxy + maxx) if topleft 

            tuple:
                (miny - maxx, maxy - minx) if not topleft

        """

        margin = self.width // np.sqrt(2)
        x, y = self.center[0], self.center[1]
        minx, miny = x - margin, y - margin
        maxx, maxy = x + margin, y + margin
        if self.topleft:
            return (miny + minx, maxy + maxx)
        else:
            return (miny - maxx, maxy - minx)

    def is_within_edges(self, x, y):
        """
        Checks if the diagonal is within the edges of the . 

        Returns 
        -------
        bool:
            True if within edges

        """
        minc, maxc = self.intercepts
        if self.topleft:
            return y >= minc - x and y <= maxc - x
        else:
            return y >= x + minc and y <= x + maxc

    @property
    def is_within_grid(self):
        """
        Checks if the object is initiated within the grid.

        Returns
        -------
        bool:
            True if within grid

        """
        minc, maxc = self.intercepts
        maxx, maxy = self.pos.shape
        if self.topleft:
            return maxy > maxc - maxx and 0 < minc
        else:
            return maxy > maxc and 0 < maxx + minc

class CenterRectangle(Object):
    """
    A class that represents the CenterRectangle object.

    '''

    Attributes
    ----------
        name : str

        ep : int

        sh : int

        dims : tuple
            List of normalised lengths of the rectangle, e.g. for 200x50 rect in 500x500 grid, len = (0.4, 0.1).

        center : tuple
            List of coordinates of the center of the rectangle, normalised between 0 - 1.

        exceed : bool
            Default is False, 

    Methods
    -------
        _set_pos():

        edges():

        edgepoints():

        is_within_grid():


    """
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
        """
        Returns a list of tuples in [(min, max), ...] as the edges of a rectangle.

        """

        return list(map(lambda center, dim: (int(center - dim//2), int(center + dim//2)), self.center, self.dims))

    @property    
    def edgepoints(self):
        """
        Returns a tuple of the coordinates of the edges of the Centre Rectangle. 

        """
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
        """
        
        """
        shape = self.pos.shape
        for i in range(len(shape)):
            if self.edges[i][0] < 0 or self.edges[i][1] > shape[i]:
                return False
        return True

class FancyRectangle2D(Object):
    """
    A class that represents the FancyRectangle object.
    Makes the rectangle easily adaptable, and rotate to any angle desired for the simulation.
    
    '''

    Attributes
    ----------
        name : str

        ep : int

        sh : int

        dims : tuple
            tuple of normalised lengths, e.g. for 200x50 rect in 500x500 grid, len = (0.4, 0.1))

        center : tuple

        angle : float
            angle of rotation in radians, anticlockwise from the x-axis, where dims is (x, y)

        
    """
    #required dims: dims (tuple/list of normalised lengths, e.g. for 200x50 rect in 500x500 grid, len = (0.4, 0.1))
    #               center (tuple/list of coords also normalised 0~1)
    #               angle of rotation in radians: anticlockwise rotation from x-axis, where dims is (x,y)
    #center and dims will round down to closest grid index, object may end up smaller than expected if input is not precise
    def __init__(self, name, ep, sh, dims, center, angle, exceed = False):
        """
        
        """
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
xspacing = [500,500]
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

diag1 = Diag2D('diag1', 0.5, sh, 0.05, (0.5, 0.5), topleft = True)
# print(diag1.intercepts)
# print(diag1.center)
assert diag1.is_within_grid
# diag1.show()

epsil = Epsil2D(sh, 1)
epsil.add_object(diag1)
# epsil.show()
