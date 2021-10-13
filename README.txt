To run the simulation open run_me.py and run it.
Starting parameters are already entered but the docstring explains how to change 
it if you wish.

There should be 4 modules for the simulation:
Object.py
Container.py
Ball.py
Simulation.py

And also run_me.py to run it.

Initially the container was going to inherit from the ball class. In this case, 
it would have had to have a very large mass. To be more accurate for very many 
collisions with very many balls I chose to have an object class they both 
inherited from so the container could have infinite mass.
