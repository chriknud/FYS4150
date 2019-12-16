import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def readfile(inputfile, number_of_objects, integration_points):
    file = open(inputfile, 'r')
    planetpositions = np.empty([number_of_objects, 3, integration_points])
    for j in range (integration_points):
        for i in range (number_of_objects):
            values = file.readline().split()

            planetpositions[i,0, j] = float(values[0])
            planetpositions[i,1, j] = float(values[1])
            planetpositions[i,2, j] = float(values[2])

    return planetpositions


planets = ['Mercury','Venus','Earth','Mars','Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']

def plotter_with_sun(inputfile, number_of_objects, integration_points, Label):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    i = 0
    plt.plot(planetpositions[0,0,:], planetpositions[0,1,:], 'o',label = 'Sun')
    for j in range (len(planetpositions)-1):
        i+=1
        if number_of_objects == 10:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = planets[i-1])
        else:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = Label)
    plt.legend(loc='best', fontsize = 23)
    #plt.show()

def plotter(inputfile, number_of_objects, integration_points, Label):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    i = 0
    #plt.plot(planetpositions[0,0,:], planetpositions[0,1,:], 'o',label = 'Sun')
    for j in range (len(planetpositions)-1):
        i+=1
        if number_of_objects == 10:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = planets[i-1])
        else:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = Label)
    plt.legend(loc='best', fontsize = 23)
    #plt.show()

def plot_angular_energy(inputfile1,inputfile2, year,dt):
    angular_momentum_verlet = np.genfromtxt(inputfile1, usecols = 0)
    angular_momentum_euler = np.genfromtxt(inputfile2, usecols = 0)

    energy_verlet = np.genfromtxt(inputfile1, usecols = 1)
    energy_euler = np.genfromtxt(inputfile2, usecols = 1)

    timestep = np.linspace(0,year, len(energy_verlet))
    plt.plot(timestep,angular_momentum_verlet, label = 'Verlet')
    plt.plot(timestep,angular_momentum_euler, label = 'Euler')

    plt.legend(loc='best', fontsize = 23)
    plt.title('Angular momentum over 2 years for Verlet and Euler', fontsize = 23)
    plt.xlabel('Time',fontsize = 23)
    plt.ylabel('Angular momentum',fontsize = 23)

    plt.show()

    plt.plot(timestep,energy_verlet, label = 'Verlet')
    plt.plot(timestep,energy_euler, label = 'Euler')
    plt.title('Total energy over 2 years', fontsize = 23)
    plt.legend(loc='best', fontsize = 23)
    plt.xlabel('Time',fontsize = 23)
    plt.ylabel('Energy',fontsize = 23)

    plt.show()


years = 2
dt0 = 1e0
dt1 = 1e-1
dt2 = 1e-2
dt3 = 1e-3
dt4 = 1e-4
dt5 = 1e-5

number_of_planets = 2

'''plotter_with_sun('positions_Earth_Sun_e1.xyz', number_of_planets, int(years/dt1)-1, Label = 'dt = 1e-1')
plotter('positions_Earth_Sun_e2.xyz', number_of_planets, int(years/dt2)-1, Label = 'dt = 1e-2')
plotter('positions_Earth_Sun_e3.xyz', number_of_planets, int(years/dt3)-1, Label = 'dt = 1e-3')
plotter('positions_Earth_Sun_e4.xyz', number_of_planets, int(years/dt4)-1, Label = 'dt = 1e-4')
plotter('positions_Earth_Sun_e5.xyz', number_of_planets, int(years/dt5)-1, Label = 'dt = 1e-5')
plt.title ('Sun Earth system with different dt - Verlet', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)

plt.show()'''

#plot_angular_energy('energy_and_momentum_verlet.xyz', 'energy_and_momentum_euler.xyz', years,dt4)

years2 = 2
'''plotter_with_sun('positions_vel04.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi * 0.4$')
plotter('positions_vel08.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi * 0.8$')
plotter('positions_velpi.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi$')
plotter('positions_vel12.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi * 1.2$')
plotter('positions_vel141.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi * 1.41$')
plotter('positions_vel16.xyz', number_of_planets, int(years2/dt4)-1, Label = r'Velocity = $2\pi * 1.6$')
plt.title ('Sun Earth system with different initial velocities', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)

plt.show()'''

years3 = 1
'''plotter_with_sun('beta30.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 2.0$')
plotter('beta32.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 2.2$')
plotter('beta34.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 2.4$')
plotter('beta36.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 2.6$')
plotter('beta38.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 2.8$')
plotter('beta40.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 3.0$')
#plotter('beta100.xyz', number_of_planets, int(years3/dt4)-1, Label = r'$\beta = 100.0$')

plt.title (r'Sun Earth system with different $\beta$', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)

plt.show()
'''
def plotter_with_sun_and_planets(inputfile, number_of_objects, integration_points, Labels):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    i = 0
    plt.plot(planetpositions[0,0,:], planetpositions[0,1,:], 'o',label = 'Sun')
    for j in range (len(planetpositions)-1):
        i+=1
        if number_of_objects == 10:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = planets[i-1])
        else:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = Labels[j])
    plt.legend(loc='best', fontsize = 23)
"""
years4 = 30
Labels = ['Earth', 'Jupiter']

plotter_with_sun_and_planets('Earth_Jupiter_Sun_fixed.xyz', 3, int(years4/dt5)-1, Labels)
plt.title (r'Sun Earth Jupiter system, non-fixed sun', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)
plt.show()

plotter_with_sun_and_planets('Earth_Jupiter_Sun10.xyz', 3, int(years4/dt5)-1, Labels)
plt.title (r'Sun Earth Jupiter system with 10 times the mass of Jupiter', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)
plt.show()

plotter_with_sun_and_planets('Earth_Jupiter_Sun1000.xyz', 3, int(years4/dt5)-1, Labels)
plt.title (r'Sun Earth Jupiter system with 1000 times the mass of Jupiter', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)
plt.show()"""

'''plotter_with_sun('Solar_system.xyz', 10, int(17/dt4)-1, Label = 'dt = 1e-1')
plt.title ('Solar System', fontsize = 23)
plt.xlabel('x-position [AU]', fontsize = 23)
plt.ylabel('y-position [AU]', fontsize = 23)

plt.show()'''

planets2 = ['Sun','Mercury','Venus','Earth','Mars','Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']
dt5 = 1e-5
def plotter3d(inputfile, number_of_objects, integration_points):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    for i in range (len(planetpositions)):
        ax.plot(planetpositions[i, 0, :],planetpositions[i, 1, :], planetpositions[i, 2, :], label = planets2[i] )
    ax.legend(loc = 'best', fontsize = 23)
    ax.set_xlabel('x-position [AU]', fontsize = 23)
    ax.set_ylabel('y-position [AU]', fontsize = 23)
    ax.set_zlabel('z-position [AU]', fontsize = 23)
    plt.title('Solar system' , fontsize = 23)
    plt.show()
plotter3d('SolarSystem.xyz', 10, int(200/dt5)-1)
