import matplotlib.pyplot as plt
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

def plotter(inputfile, number_of_objects, integration_points):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    i = 0
    plt.plot(planetpositions[0,0,:], planetpositions[0,1,:], 'o',label = 'Sun')
    for j in range (len(planetpositions)-1):
        i+=1
        if number_of_objects == 10:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = planets[i-1])
        else:
            plt.plot(planetpositions[i,0,:], planetpositions[i,1,:], label = 'object:{}'.format(i))
    plt.legend(loc='best', fontsize = 17)
    plt.show()
years = 80
dt = 1e-4
plotter('positions.xyz', 10, int(years/dt)-1)
