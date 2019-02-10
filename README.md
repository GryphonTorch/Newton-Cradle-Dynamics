# Newton-Cradle-Dynamics
This python 2.7 code simulates the motion of balls in a Newton's Cradle, a common physics toy used to demonstrate momentum conservation. 
In real-life, the balls lose energy in collisions and eventually come to a stop. 
The code uses Euler-Richardson numerical integration to calculate the balls' paths, subject to a Hertz-Kurobawa viscoelastic damping force.

Air drag dissipation is included. The collision model is similar to a damped spring, but the contact force is proportioanl to the displacement^1.5. This exponent, rather than 1, takes into account the spherical nature of the colliding surfaces where we expect an increasingly larger repulsion force as the balls impinge on each other. 

Object-oriented programming used so as to allow for easy extension of balls, to 3,4,5 .etc balls. 
Prepared for the 2019 International Young Physicists' Tournament Problem 15, Newton's Cradle (http://old.iypt.org/images/9/9d/problems2019_signed.pdf). 
