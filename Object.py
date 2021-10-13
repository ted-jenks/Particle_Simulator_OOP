
"""
Module for object class
Edward Jenks 2nd year computing

Initially the container was going to inherit from the ball class. In this case, 
it would have had to have a very large mass. To be more accurate for very many 
collisions with very many balls I chose to have an object class they both 
inherited from so the container could have infinite mass.
"""
#%%
import scipy as sp
import pylab as pl
#%%

class object:
    
    """
    Object class
    (mass, radius, position, velocity)
    This class contains the information the ball and the container can inherit from
    """
    
    def __init__(self, mass, radius, p1, p2, v1, v2):
        '''
        Initialises an object.
        
        mass, radius, position and velocity can be selected. For calculations, 
        radius must be a float. If it is not, an error is raised.
        '''
        self._m = mass
        self._R = radius
        self._v = sp.array([v1,v2])
        self._r = sp.array([p1,p2])
        if not isinstance(self._R,float):
            raise TypeError('Radius must be a float.')
    
    def pos(self):
        '''
        Returns the object's position.
        '''
        return self._r
    
    
    def vel(self):
        '''
        Returns the object's velocity.
        '''
        return self._v
    
    def get_patch(self):
        '''
        Returns the object's patch.
        '''
        if self._R >= 10:
            return pl.Circle(self._r, self._R, ec = 'Black', 
                             fill = False, ls='solid')
        else:
            return pl.Circle(self._r, self._R, fc='darkred')
            
    def __repr__(self):
        '''
        Returns information about the object.
        '''
        if self._R >= 10:
            return ("%s(Mass=%s, Radius=%s, Position=%s)" % 
                    ("Container", self._m, self._R, self._r))
        else:
            return ("%s(Mass=%s, Radius=%s, Position=%s, Velocity=%s)" % 
                    ("Ball", self._m, self._R, self._r, self._v))
    
    def __str__(self):
        '''
        Returns information about the object.
        '''
        if self._R >= 10:
            return ("(M=%g, R=%g, P=%gi+%sj)" % 
                    (self._m, self._R, self._r[0],self._r[1], self._v[0], self._v[1]))
        else:
            return ("(M=%g, R=%g, P=%gi+%sj, V=%si,%sj)" % 
                    (self._m, self._R, self._r[0],self._r[1], self._v[0], self._v[1]))
   
