import numpy as np
import time

#Only changes to be made to yapic.py is to copy 'our stuff' in, initialise epsilon and mirror object finally add * epsilon behind the *dus
#Initialise tgt with E and B dicts

#our stuff
class Epsil2D():
    def __init__(self, sh, bgep = 1):
        #Takes in shape(dimens) of grid and background epsil value (air)
        self._objects = {}
        self.bgep = 1
        self._epsil = np.full(sh, self.bgep) 

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
        if revert:
            pass
            #code to revert values in self._epsil according to object._pos to bgep
        else:
            pass
            #code to change values in self._epsil according to object._pos to object.ep


class Object():
    def __init__(self, ep, sh, exceed = True):
        #initialise with epsil value for object and a grid with same shape as background
        #Child classes will override, calls super().__init__, store additional arguments and call _set_pos with those args
        #Additional args will be dimensions of object with values btwn 0 and 1
        #Object created with dims/size relative to grid
        self.ep = ep
        self._pos = np.full(sh, False)
        self.exceed = exceed

    def _set_pos(self):
        #To be called in __init__ in child classes
        #If self.exceed == True, allow for making object that exceeds grid but cuts off
            #Use is_within_grid to adjust algorithm for cutting off
        #Else raise error if object exceeds grid, use is_within_grid
        pass

    @property
    def edgepoints(self):
        #might not need
        #returns tuple of min/max points in each dim to ensure is within range
        #Values are array indices
        pass

    @property
    def is_within_grid(self):
        #uses dims of objects to calc edge points to check if within grid
        #Returns true if within grid
        pass

class Ractangle(Object):
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