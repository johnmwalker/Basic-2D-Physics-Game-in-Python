# Basic-2D-Physics-Game-in-Python

This is the final draft of the game I made for PHYS 398 in Spring 2016
(sophomore year). I have not bothered to update it much since then. See other 
projects for better coding.

To play, run the poorly-named TestStation.py file.

Requires:

PyGame

Features: 

N-body gravity
Basic collisions
RK4 Solver

Goal of game: 

Suck the energy out of the planets to make them fall into the black hole.
You'll lose if you fall in yourself.

Controls: 

Numpad + and - speed up and slow down the simulation.
> and < (. and ,) are "cheats" to give your thruster more power.
Spacebar activates the "energy sucking" function. Planets in the crosshairs will slow down.
Arrow keys control motion of ship.
    
Issues:

I couldn't figure out how to fix the rotation animation of the ship, so it
doesn't really rotate about the center of the image. I think this may have
been due to a bug in PyGame.

Collisions were poorly implemented. For some reason I decided to use trig
to calculate forces. Later projects built on this code use vectors for
greatly improved accuracy.

This was my first big project, so many coding conventions were likely
ignored.
