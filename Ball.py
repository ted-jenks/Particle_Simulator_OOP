
"""
Module for ball class
Edward Jenks 2nd year computing
"""

#%%
import scipy as sp
import Container as ct
import Object as ob
#%%

class ball(ob.object): #Inherits from object class
    
    """
    Ball class
    (self, mass, radius, x position, y position, x velocity, y velocity)
    This object contains information for the balls.
    It inherits from the object class.
    """
    
    Collisions = 0 #Counts the total collisions that have taken place
    Time = 0 #Counts the time the simulation has been running for
    MomChangeTot = 0 #Stores the total change in momentum
    TempTime = 0 #A reseting time
    TempMomChangeTot = 0 #A reseting momentum change
    PressureOverTime = [] #Pressure stored every 10s
    
    def __init__(self, mass, radius, p1, p2, v1, v2):
        '''
        Initialises a ball.
        
        Values for mass, radius, position and velocity can be selected.
        '''
        ob.object.__init__(self, mass, radius, p1, p2, v1, v2)
    
    
    def move(self, dt):
        '''
        moves the ball with constant velocity for dt.
        '''
        self._r += self._v*dt
    
        
    def time_to_collision(self,other):
        '''
        Returns the time to collision with another object.
        '''
        if isinstance(other, ct.container):
            R = other._R - self._R
        else:
            R = self._R + other._R
        
        r = self._r - other._r
        v = self._v - other._v
            
        a = sp.dot(v,v)
        b = 2 * sp.dot(r,v)
        c = (sp.dot(r,r)-R**2)
        dt1 = (-b + sp.sqrt(b**2 -4*a*c))/(2*a)
        dt2 = (-b - sp.sqrt(b**2 -4*a*c))/(2*a)
        dtT = []  #Stores both possible dts
        dtT.append(dt1)
        dtT.append(dt2)
        for n in range(2):
            if dtT[n] <= 0 or isinstance(dtT[n],complex):
                dtT[n] = 10e40
        dtT.sort()  #Finds minimum dt
        dt = dtT[0]
        return (dt)
    
    
    def collide(self, other):
        '''
        Collides the ball with another object.
        
        Also calculates momentum change after each collision to ensure KE and
        momentum are conserved. If they are not an erro will be raised.
        '''
        if isinstance(other, ct.container):   #Collisions with container
            k1 = self._m * (self._v[0]**2+self._v[1]**2)
            r21 = self._r
            Dir = ((r21)/sp.sqrt( r21[0]**2 + r21[1]**2 ))
            V1_init = sp.dot(self._v,Dir)
            V1_finmag = -V1_init
            V1_perp = self._v - (V1_init*Dir)
            V1_fin = V1_finmag*Dir
            self._v = V1_perp + V1_fin
            k2 = self._m * (self._v[0]**2+self._v[1]**2)
            MomChange = 2 * self._m * abs(V1_init)
            
            if ball.Time >= 200:
                ball.MomChangeTot += MomChange
            ball.TempMomChangeTot += MomChange
        else:   #Collisions with balls
            k1 = ((self._m * (self._v[0]**2+self._v[1]**2)) + 
                  (other._m * (other._v[0]**2+other._v[1]**2)))
            r21 = other._r-self._r
            Dir = ((r21)/sp.sqrt( r21[0]**2 + r21[1]**2 ))
            V1_init = sp.dot(self._v,Dir)
            V2_init = sp.dot(other._v,Dir)
            V2_finmag = (((2*self._m)/(self._m+other._m))*V1_init - 
                        ((self._m-other._m)/(self._m+other._m))*V2_init)
            V1_finmag = (((self._m-other._m)/(self._m+other._m))*V1_init + 
                        ((2*self._m)/(self._m+other._m))*V2_init)
            V1_perp = self._v - (V1_init*Dir)
            V2_perp = other._v - (V2_init*Dir)
            V1_fin = V1_finmag*Dir
            V2_fin = V2_finmag*Dir
            self._v = V1_perp + V1_fin
            other._v = V2_perp + V2_fin
            k2 = ((self._m * (self._v[0]**2+self._v[1]**2)) + 
                  (other._m * (other._v[0]**2+other._v[1]**2)))   
        ball.Collisions += 1
        print(ball.Collisions)
        if round(k1,6) != round(k2,6):
            TypeError('Kinetic energy and momentum not conserved')
    
                
