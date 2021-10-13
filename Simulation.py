'''
Module to hold simulation class
Edward Jenks 2nd year computing

This is a simulation of atoms in a circluar container. It can be used to demonstrate
thermodynamic phenominon such as the ideal gas law and the Maxwell-Boltzmann 
distribution. It also confirms the conservation of momentum and kinetic energy.
'''
#%%
import scipy as sp
import Ball as bl
import Container as ct
import pylab as pl
import matplotlib.pyplot as plt
import random
#%%

class simulation:
    
    """
    Simulation class
    (container radius, ball mass, ball radiu, number of balls)
    This class runs and analyzes the simulation.
    """
    
    init = 0 #Changes to 1 after the original dt matrix is formed
    Radius = 0 #The radius of the container
    Extensions = [] #List of ball distances from center
    Separations = [] #List of bal to ball distances
    KineticEnergies = [] #List of ball kinetic energies
    Speeds = [] #List of ball speeds
    TotalKE = 0 #Total KE of the system
    Temperature = 0 #Temperature of the system
    Pressure = 0 #Pressure of the system
    
    def __init__(self,radiusC,mass, radiusB, N, Vmag):
        '''
        Initialises the simulation
        
        Values for the radius of the container, mass of the balls, radius of the balls 
        and number of balls must be input. 
        
        The balls are given random velocities of average ~ 0 and are placed in a 
        grid that scales depending on ball radius andthe number of balls.
        '''
        self._ball = []
        self._N = N
        self._m = mass
        p = sp.linspace(-5,5, int(round((10/(3*radiusB)))))
        if self._N > int(round((10/(3*radiusB))))**2:
            raise TypeError('Too many balls. Reduce N or ball radius.')
        n = 0
        x = -1
        for i in range(self._N):
            if x <= (int(round((10/(3*radiusB))))) - 2:
                x += 1
                ball = (bl.ball(mass, radiusB, p[x], p[n], 
                        Vmag*random.randrange(-100,100)/100, 
                        Vmag*random.randrange(-100,100)/100))
                self._ball.append(ball)
            else:
                x = 0
                n += 1
                ball = (bl.ball(mass, radiusB, p[x], p[n], 
                        Vmag*random.randrange(-100,100)/100, 
                        Vmag*random.randrange(-100,100)/100))
                self._ball.append(ball)                                               
        simulation.Radius = radiusC
        self._cont = ct.container(radiusC)
    
        
    def next_collision(self):
        '''
        Calculates the time until the next collisions.
    
        The first time the method is run all possible collision times are calculated 
        and stored in a matrix. After that, only the two balls that collided have their
        collision times recalculated for efficiency.
    
        dt is the minimum of the matrix and the position of it gives the objects involved.
        '''   
        for i in range(self._N):
            self._ball[i].move(0.00001)
            
        dtall = sp.zeros(shape = [self._N, self._N + 1]) #Creates matrix to store collision times
        
        if simulation.init == 0:
            for i in range(self._N):
                dtall[i][self._N] = (self._ball[i].time_to_collision(self._cont))       
                for n in range(i+1, self._N):
                    dtall[i][n] = (self._ball[i].time_to_collision(self._ball[n]))      
            for i in range(self._N):
                for n in range(self._N):
                    if dtall[i][n] == 0:
                        dtall[i][n] = 10e40 #Get rid of non calculated collisions
            simulation.dtall = dtall
            simulation.init +=1
            
        else:   #Only recalculates dts for the two that just collided
            simulation.dtall = simulation.dtall -simulation.dtp - 0.00001
            simulation.dtall[simulation.lasti][self._N] = self._ball[simulation.lasti].time_to_collision(self._cont)
            for n in range(self._N):
                if n != simulation.lasti:
                    simulation.dtall[simulation.lasti][n] = self._ball[simulation.lasti].time_to_collision(self._ball[n])
            
            if simulation.lastj != self._N:
                simulation.dtall[simulation.lastj][self._N] = self._ball[simulation.lastj].time_to_collision(self._cont)
                for n in range(self._N):
                    if n != simulation.lastj:
                        simulation.dtall[simulation.lastj][n] = self._ball[simulation.lastj].time_to_collision(self._ball[n])
            for y in range(self._N):
                for x in range(self._N+1):
                    if simulation.dtall[y][x] == 0 or simulation.dtall[y][x]==0.0:
                        simulation.dtall[y][x] = 10e40
        dt = simulation.dtall.min()
        simulation.dtp = dt                          #Store last dt
        bl.ball.Time += (dt)                         #Add time to total time   
        bl.ball.TempTime += (dt)
        for i in range(self._N):
            self._ball[i].move(dt)                   #Move balls
        ij_min = sp.where(simulation.dtall==dt)      #Location of minimum dt
        i, j=ij_min[0][0], ij_min[1][0]
        simulation.lasti, simulation.lastj = i, j
        self.extension()
        self.separation()
        self.pressure_in_time()
        self.ke()
        if j == self._N:
            self._ball[i].collide(self._cont)        #Carry out container collision
        else:
            self._ball[i].collide(self._ball[j])     #Carry out ball collision
                                                       
              
    def extension(self):
        '''
        Stores the positions of the balls relative to the container at each collision.
        '''
        if bl.ball.Time > 200:    
            for i in range(self._N):
                Distance = sp.sqrt(self._ball[i]._r[0]**2 + self._ball[i]._r[1] **2)    
                simulation.Extensions.append(Distance)
    
    
    def separation(self):
        '''
        Stores the distances between balls at each collision.
        '''
        if bl.ball.Time > 200:
            for i in range(self._N):
                for n in range(i+1, self._N):
                    r = (self._ball[i]._r - self._ball[n]._r)                           
                    Distance = sp.sqrt(r[0]**2 + r[1] **2)      
                    simulation.Separations.append(Distance)
   
    
    def pressure_in_time(self): 
        '''
        Stores the pressure at a given time.
        '''                                                        
        if bl.ball.TempTime >= 0.5:
            pressure = (bl.ball.TempMomChangeTot / (bl.ball.Time * 
                        2.0 * 3.141592653589793 * simulation.Radius))
            bl.ball.PressureOverTime.append(pressure)
            bl.ball.TempTime = 0


    def ke(self): #Method to find the kinetic energy of the balls
        '''
        Stores the kinetic energy of each ball at each collision
        '''
        if bl.ball.Time > 200:
            Total = 0
            for n in range(self._N):
                s = sp.sqrt(self._ball[n]._v[0] **2 + self._ball[n]._v[1] **2)
                k = 0.5 * self._ball[n]._m * s**2
                Total += k
                simulation.Speeds.append(s)
                simulation.KineticEnergies.append(k) 
            simulation.TotalKE = Total
   
        
    def pressure_time_graph(self):
        '''
        Plots a graph of pressure development with time.
        '''
        plt.figure()                                                                    
        x = sp.linspace(0,bl.ball.Time,len(bl.ball.PressureOverTime))
        plt.title('Pressure Development with Time')
        plt.xlabel('Time(s)')
        plt.ylabel('Pressure(Pa)')
        plt.plot(x, bl.ball.PressureOverTime)
    
        
    def info(self):
        '''
        Prints the Pressure, Temperature, total kinetic energy and Maxwell-Boltzmann
        varience of the system after the simulation has run.
        '''
        simulation.Pressure = (bl.ball.MomChangeTot / ((bl.ball.Time-200) * 
                                2.0 * 3.141592653589793 * simulation.Radius))
        print('Average pressure:',simulation.Pressure)
        simulation.Temperature = simulation.TotalKE / (self._N * 1.380648e-23)
        print('Temperaure:',simulation.Temperature)
        print('Total Kinetic Energy:', simulation.TotalKE)
        print('Maxwell-Boltzmann Varience:', sp.var(simulation.Speeds))
        print(bl.ball.Time)
        Er = ((self._N * 1.380648e-23 * simulation.Temperature)/
              (simulation.Pressure * 3.141592653589793 * simulation.Radius**2))
        print('NkT/PV = ', Er)
    
        
    def extensions_graph(self):
        '''
        Plots the positions of the balls relative to the container at each collision.
        '''
        plt.figure()                                                                
        numbins = 19
        plt.title('Ball Extensions from Centre')
        plt.xlabel('Distance (m)')
        plt.ylabel('Frequancy Density')
        n, bins, patches = (plt.hist(simulation.Extensions, numbins, 
                            facecolor = 'b', alpha = 0.5, histtype = 'bar', ec = 'k'))
        plt.show()
    
    
    def separations_graph(self):
        '''
        Stores the positions of the balls relative to the container at each collision.
        '''
        plt.figure()                                                                   
        numbins = 18
        plt.title('Ball Separations')
        plt.xlabel('Distance (m)')
        plt.ylabel('Frequancy Density')
        n, bins, patches = plt.hist(simulation.Separations, numbins, 
                                    facecolor = 'b', alpha = 0.5, histtype = 'bar', ec = 'k')
        plt.show()
    
        
    def ke_graph(self):
        '''
        Plots the frequancy of kinetic energies in the balls.
        '''
        plt.figure()                                                                  
        numbins = 50
        plt.title('Ball Kinetic Energies')
        plt.xlabel('Energy (J)')
        plt.ylabel('Frequancy Density')
        n, bins, patches = plt.hist(simulation.KineticEnergies, numbins, 
                                 facecolor = 'b', alpha = 0.5, histtype = 'bar', ec = 'k')
        plt.xlim(0,3)
        plt.show()
    
       
    def MB_graph(self):
        '''
        Plots the probability of finding balls at certain velocities. Pronted over is
        the predicted result, the Maxwell-Boltzmann distribution.
        '''
        plt.figure()                                               
        PDF_x = sp.linspace(0,30,20000)
        PDF_y = PDF_x * sp.exp(-(0.5*self._m*PDF_x**2)/(1.380648e-23*simulation.Temperature)) 
        simulation.expvar = sp.var(PDF_y)
        numbins = 40
        weights = 5*sp.ones_like(simulation.Speeds)/len(simulation.Speeds)
        plt.title('Maxwell Boltzman Distribution')
        n, bins, patches = (plt.hist(simulation.Speeds, numbins,facecolor = 'b', 
                           alpha = 0.5, histtype = 'bar', ec = 'k', 
                           weights =weights, label = 'Data'))
        plt.plot(PDF_x,PDF_y,"r--", label = 'Maxwell-Boltzmann Distribution', 
                 linewidth = 3)            #Maxwell Bolztaman Trendline
        plt.legend()
        plt.xlabel('Speed (m/s)')
        plt.ylabel('Probability')
        plt.xlim(0,10)
        plt.show()
    
    
    def show_graphs(self):   #Method to plot the data gathered above
        '''
        Runs the methods to generate graphs. can be edited depending on what graphs 
        you want shown.
        '''
        self.info()
        self.pressure_time_graph()
        self.extensions_graph()
        self.separations_graph()
        self.ke_graph() #Only for large numbers of collisions
        self.MB_graph() #Only for large number of collisions
        print(simulation.expvar)
    
    
    def run(self, num_frames, animate=False, show_graphs=False): #Choose to show animation
        '''
        Runs the simulation. number of frames must be chosen, this will be equal to 
        the number of collisions you want to run plus one. Additionally, you can choose 
        not to run the animation for efficiencey when running many balls for many collisions.
        '''
        if animate:
            plt.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_artist(self._cont.get_patch())
            for i in range(self._N):
                ax.add_patch(self._ball[i].get_patch())
        for frame in range(num_frames):
            self.next_collision()
            if animate:
                pl.pause(.01)
        if animate:
            pl.show()
        if show_graphs:
            self.show_graphs()
    


            
            