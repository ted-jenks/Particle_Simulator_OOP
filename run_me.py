
"""
Edward Jenks
2nd year computing

Simulation of a thermodynamic system.

This module is for you you experiment with the code. Various parameters of the 
particles and the container can be altered from here and run.
"""
import Simulation as sim


sim = sim.simulation(10.0,1,0.3,50,1)
'''
Simulation class
Takes values (container radius, ball mass, ball radius, number of balls, Max velocity).
'''
sim.run(10000, animate=False, show_graphs=True)
'''
If any parameters are entered incorrectly or you try to load too many balls for the 
radius you have selected an error should be raised.

As the simulation waits for thermal equilibrium to be reached before taking data, 
only use show_graphs for more than
'''