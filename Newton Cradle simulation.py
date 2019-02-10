"""
Newton's Cradle Simulation

EULER-RICHARDSON METHOD FOR VELOCITY DEPENDENT ACCELERATION with ADAPTIVE TIME
This simulation collides three balls using the Hertzian spring model with 
air drag and viscoelastic dissipation terms included. Note over time, the balls
start to move in phase (clumped together) due to the viscoelastic damping. 
"""
import numpy
import matplotlib.pyplot as plt
import csv          # to export data

#Physics constants depending on model
g = 9.81            # gravitational acceleration (m/s^2)
L = 0.35             # length of the pendulums (m)
alpha = 1.5         # Hertzian collision power parameter
beta = 1.5          # Hertz-Kuwabara model viscoelastic dissipation power
E = 200*10**9       # Young's modulus of steel ball (200 Gigapascal)
R = 0.025           # radius of steel ball, say 2.5 cm
v = 0.30            # Poisson's ratio, say 0.30
m = 0.05
k1 = 1.01E10  #(2*R)**0.5*E/(3*(1-v**2))    # spring constant; see paper
k2 = m*g/L                        # gravitational restoring force constant mg/L
eta = 0.00107261                      # air drag constant; now from paper
gamma = 2600                       # viscoelastic dissipation constant

timescale = (m**2/(k1**2*0.5))**0.2/10
print "Recommended timescale, from collisions, is", timescale

#Computational constants
deltaT = 0.000001     # duration of each timestep (seconds)
count = 2000000       # for time, this is an integer
times = []           # to plot for time elapsed

"1. =====Making Balls class====="
class Ball(object):
    """class for building different balls and their data required.
    Each Ball object collects floats of its center position, initial (equilibrium) 
    position, linear velocity, mass and radius. 
    For defensive programming, only access data attributes through the various 
    set and get functions defined below to prevent inadvertent changes."""    
    def __init__(self, pos0, posTop, vel0, mass, radius):
        self.position = pos0        # initial pendulum ball x-position, raised up, in m
        self.positionTop = posTop   # initial string top position, in m     
        self.velocity = vel0        # initial velocity 
        self.mass = mass            # mass, in kg
        self.radius = radius        # radius, in m
        self.results = []           # to export data into csv (a list of lists)
        self.epsilon = 0            # this ball's overlap with the adjacent left
    
    def setPos(self, newPos):
        """sets position to the given input"""
        self.position = newPos
    
    def setVel(self, newVel):
        """sets velocity to the given input"""
        self.velocity = newVel
    
    def setEpsilonRight(self, newEp):
        """sets epsilon (overlap) to the given input, assuming this ball
        is on the right??"""
        self.epsilon = newEp        
    
    def setPath(self, time, newPos):
        """allows a position and time list to be appended to the results list
        a list of list, for convenience when exported to csv file"""
        self.results.append([time,newPos]) 
    
    def getPos(self):
        """returns ball center's current horizontal displacement"""
        return self.position
    
    def getPosTop(self):
        """returns ball center's initial horizontal displacement"""
        return self.positionTop
        
    def getVel(self):
        """returns ball's current velocity"""
        return self.velocity
    
    def getEpsilonRight(self):
        """returns epsilon (overlap) with the left, assuming this ball
        is on the right"""
        return self.epsilon           
    
    def getRadius(self):
        """returns ball's radius"""
        return self.radius    
        
    def getMass(self):
        """returns ball's mass"""
        return self.mass    
       
    def getPath(self):
        """returns trajectory as list"""
        return self.results
    
"2. =====Setting up the scene====="
# Create the pendulum bobs
#Ball_One = Ball(-0.0527 , -0.017, 0.70, 0.067, 0.0125)      # displaced to the right
#Ball_Two = Ball(0.0073 , .0079, 0, 0.067, 0.0125)    # center ball
#Ball_Three = Ball(0.0325 , 0.0325, 0, 0.067, 0.0125)  # leftmost ball ???

Ball_One = Ball(0.025 , 0.025, 0, 0.067, 0.0125)      # displaced to the right
Ball_Two = Ball(0.0 , 0.0, 0, 0.067, 0.0125)    # center ball
Ball_Three = Ball(-0.067 , -0.025, 0, 0.067, 0.0125)  # leftmost ball ???

#initialise some empty lists, for plots later
pos_Ball_One = []    # only position data, no time
pos_Ball_Two = []
pos_Ball_Three = []
times = []
epsilons12 = []        # for debugging
epsilons23 = []        # for debugging
epsilons012 = []       # checking
visc012 = [] #check

energies = []        # to check
collisionSteps = 0   # initialise number of colisions, an int

#x0 notation denotes equilibrium position of Balls, which are constant.
x0_Ball_One = Ball_One.getPosTop()  
x0_Ball_Two = Ball_Two.getPosTop()  
x0_Ball_Three = Ball_Three.getPosTop()  

timeNow = 0

"3. =====Action=====" 
#refer to Verlet integration algorithm wikipedia page for more info. 
for n in range(count):
    #x0, x, x1, x2 refer to equilibrium, current, mid and next timestep values
    """Part I: Calculate a(t)"""
    
    #set initial boundary values
    x_Ball_One = Ball_One.getPos()    # position of Ball One 
    v_Ball_One = Ball_One.getVel()    # linear velocity of Ball One
    ep0_Ball_One = Ball_One.getEpsilonRight()  #earlier overlap of One with Two
    x_Ball_Two = Ball_Two.getPos()  
    v_Ball_Two = Ball_Two.getVel()  
    ep0_Ball_Two = Ball_Two.getEpsilonRight()  #earlier overlap of Two with Three
    x_Ball_Three = Ball_Three.getPos()  
    v_Ball_Three = Ball_Three.getVel()  
    
    if abs(x_Ball_One - x_Ball_Two) > 0.026 and  abs(x_Ball_Two - x_Ball_Three) > 0.026:
        deltaT = 0.001   #let's speed things up cos balls are far away
    else:
        deltaT = 0.000001   # balls are in contact
         
    #calculate current overlaps, if any, for force, and viscoelastic dissipation
    epsilon12 = Ball_One.getRadius() + Ball_Two.getRadius() - abs(x_Ball_One - x_Ball_Two)
    if epsilon12 >= 0:
        collisionSteps += 1    # allowed value of epsilon between Balls One and Two
    else:
        epsilon12 = 0   # balls not touching, set value to 0
    epsilons12.append(epsilon12) #debugging
    epsilons012.append(ep0_Ball_One) #debugging. is it the same? 

    visc12 = gamma*(epsilon12**beta - ep0_Ball_One**beta)/deltaT  #this line has problems
    visc012.append(visc12)  # debugging
    # linearised approximation for the viscoelastic dissipation term   
    
    epsilon23 = Ball_Two.getRadius() + Ball_Three.getRadius() - abs(x_Ball_Two - x_Ball_Three)
    if epsilon23 >= 0:
        pass            # allowed value of epsilon between Balls Two and Three
    else:
        epsilon23 = 0   # balls not touching, set value to 0
    epsilons23.append(epsilon23) #debugging
    visc23 = gamma*(epsilon23**beta - ep0_Ball_Two**beta)/deltaT  # viscoelastic dissipation
    
    #next, calculate relevant accelerations at this time, with drag
    acc_Ball_One = (k1*(epsilon12**alpha) + k2*(x0_Ball_One - x_Ball_One) - \
    eta*v_Ball_One + visc12)/Ball_One.getMass()
    
    acc_Ball_Two = (-k1*(epsilon12**alpha) + k1*(epsilon23**alpha) + \
    k2*(x0_Ball_Two - x_Ball_Two) - eta*v_Ball_Two - visc12 + visc23) / Ball_Two.getMass()
    
    acc_Ball_Three = (- k1*(epsilon23**alpha) + k2*(x0_Ball_Three - x_Ball_Three) - \
    eta*v_Ball_Three - visc23)/Ball_Three.getMass()
    

    """Part II: Calculate v1, x1 midstep"""
    #now, calculate mid-step velocity and position by Euler method directly
    v1_Ball_One = v_Ball_One + 0.5*acc_Ball_One*deltaT   
    v1_Ball_Two = v_Ball_Two + 0.5*acc_Ball_Two*deltaT 
    v1_Ball_Three = v_Ball_Three + 0.5*acc_Ball_Three*deltaT 
    
    x1_Ball_One = x_Ball_One + 0.5*v1_Ball_One*deltaT
    x1_Ball_Two = x_Ball_Two + 0.5*v1_Ball_Two*deltaT
    x1_Ball_Three = x_Ball_Three + 0.5*v1_Ball_Three*deltaT


    """Part III: Calculate acc1 midstep; check the viscoelastic dissipation term"""
    #redo epsilon overlap, and viscoelastic dissipation terms
    epsilon12_new = Ball_One.getRadius() + Ball_Two.getRadius() - abs(x1_Ball_One - x1_Ball_Two)
    if epsilon12_new >= 0:
        pass            # allowed value of epsilon between Balls One and Two
    else:
        epsilon12_new = 0   # balls not touching, set value to 0  
    visc12_new = gamma*(epsilon12_new**beta - ep0_Ball_One**beta)/deltaT*2  # viscoelastic dissipation

    epsilon23_new = Ball_Two.getRadius() + Ball_Three.getRadius() - abs(x1_Ball_Two - x1_Ball_Three)
    if epsilon23_new >= 0:
        pass            # allowed value of epsilon between Balls Two and Three
    else:
        epsilon23_new = 0   # balls not touching, set value to 0
    visc23_new = gamma*(epsilon23_new**beta - ep0_Ball_Two**beta)/deltaT*2  # viscoelastic dissipation
        
    #redo accelerations
    acc1_Ball_One = (k1*(epsilon12_new**alpha) + k2*(x0_Ball_One - x_Ball_One) - \
    eta*v_Ball_One + visc12_new)/Ball_One.getMass()
    
    acc1_Ball_Two = (-k1*(epsilon12_new**alpha) + k1*(epsilon23_new**alpha) + \
    k2*(x0_Ball_Two - x_Ball_Two) - eta*v_Ball_Two - visc12_new + visc23_new) / Ball_Two.getMass()
    
    acc1_Ball_Three = (- k1*(epsilon23_new**alpha) + k2*(x0_Ball_Three - x_Ball_Three) - \
    eta*v_Ball_Three - visc23_new)/Ball_Three.getMass()    

    
    """Part IV: Calculate next timestep v2. x2"""
    v2_Ball_One = v_Ball_One + acc1_Ball_One*deltaT
    v2_Ball_Two = v_Ball_Two + acc1_Ball_Two*deltaT
    v2_Ball_Three = v_Ball_Three + acc1_Ball_Three*deltaT

    energies.append(v2_Ball_One**2 + v2_Ball_Two**2 + v2_Ball_Three**2) # check
    
    x2_Ball_One = x_Ball_One + v1_Ball_One*deltaT
    x2_Ball_Two = x_Ball_Two + v1_Ball_Two*deltaT
    x2_Ball_Three = x_Ball_Three + v1_Ball_Three*deltaT    
    
    """Part V: Update position and velocity back to Ball object"""
    Ball_One.setPos(x2_Ball_One)
    Ball_Two.setPos(x2_Ball_Two)
    Ball_Three.setPos(x2_Ball_Three)

    Ball_One.setVel(v2_Ball_One)
    Ball_Two.setVel(v2_Ball_Two)
    Ball_Three.setVel(v2_Ball_Three)

    Ball_One.setEpsilonRight(epsilon12_new)
    Ball_Two.setEpsilonRight(epsilon23_new)

    Ball_One.setPath(n*deltaT, x2_Ball_One)
    Ball_Two.setPath(n*deltaT, x2_Ball_Two)
    Ball_Three.setPath(n*deltaT, x2_Ball_Three)

    pos_Ball_One.append(x2_Ball_One)
    pos_Ball_Two.append(x2_Ball_Two)
    pos_Ball_Three.append(x2_Ball_Three)

    times.append(timeNow + deltaT)
    timeNow = timeNow + deltaT
    
    #call Next iteration
    n = n + 1
    
    p = count/10
    if n in (p,2*p,3*p,4*p,5*p,6*p,7*p,8*p,9*p):     # for reference
        print float(n)/count*100, "% calculated"

"4. =====Plot and export results====="    

print collisionSteps, "timesteps between ball one and two"

def plotPath(positions, name):
    """plots a graph of position against time"""
    plt.figure(figsize=(7,7))        
    plt.plot(times, positions)
    plt.ylabel(name)
    plt.xlabel('Time (sec)')
    plt.show()    

def writePath(Ball, name):
    "exports position data to a csv file in the folder you saved this code in"
    myFile = open(name + '.csv', 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(Ball.getPath())    
    print("File writing complete")

writePath(Ball_One, "Ball One path")
writePath(Ball_Two, "Ball Two path")
writePath(Ball_Three, "Ball Three path")

print "=====Simulation complete====="

plotPath(pos_Ball_One, 'Ball One position (m)')
plotPath(pos_Ball_Two, 'Ball Two position (m)')
plotPath(epsilons12, 'Epsilon value (m)')  # Debugging
plotPath(epsilons23, 'Epsilon value (m)')  # Debugging
plotPath(energies, 'Energies (J)')       # For now, no dissipation terms

plt.figure(figsize=(7,7))        
plt.plot(times, pos_Ball_One, 'r')
plt.plot(times, pos_Ball_Two, 'g')
plt.plot(times, pos_Ball_Three, 'b')
plt.ylabel('All Balls positions (m)')
plt.xlabel('Time (sec)')
plt.show()   