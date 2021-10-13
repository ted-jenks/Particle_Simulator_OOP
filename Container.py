
"""
Module for container class
Edward Jenks 2nd year computing
"""

#%%
import Object as ob
#%%

class container(ob.object): #Inherits from object class
    
    """
    Container class
    (radius)
    This class contains the information for the container.
    It inherits from the object class.
    """
    
    def __init__(self, radius):
        '''
        Initialises the container at 0,0 with no velocity or mass. 
        
        A value for radius can be selected.
        '''
        mass = 0
        p1= 0
        p2 = 0
        v1 = 0
        v2 = 0
        ob.object.__init__(self, mass, radius, p1, p2, v1, v2)
    