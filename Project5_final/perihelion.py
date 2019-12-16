import numpy as np
import matplotlib.pyplot as plt

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

def find_perihelion(inputfile, number_of_objects, integration_points):
    planetpositions = readfile(inputfile, number_of_objects, integration_points)
    plt.plot(planetpositions[1, 0, :],planetpositions[1, 1, :], label = 'Mercury')
    plt.plot(planetpositions[0, 0, :],planetpositions[0, 1, :], 'o', label = 'Sun')
    r = np.sqrt(planetpositions[1, 0, :]**2+planetpositions[1, 1, :]**2+planetpositions[1, 2, :]**2)
    rs = np.sqrt(planetpositions[0, 0, :]**2+planetpositions[0, 1, :]**2+planetpositions[0, 2, :]**2)
    DR = np.abs(r-rs)
    print ('DR0',DR[0])
    round = 0
    p = 0
    dr = np.abs(planetpositions[1,1,0])
    eps = 1e-4
    x_p = []
    y_p = []
    p = 0.4
    for i in range (len(planetpositions[0,0,:])-1):
        if planetpositions[1,1,i-1]>0 and planetpositions[1,1,i]<0 and planetpositions[1,0,i]<0:
            round +=1

        if planetpositions[1,1,i]<0 and  planetpositions[1,0,i]<0:
            if np.abs(DR[i])<=p:
                #print ('hei', i)
                p = DR[i]
                new_x = planetpositions[1,0,i]
                new_y = planetpositions[1,1,i]


            if planetpositions[1,0,i]<0 and planetpositions[1,0,i+1]>0:
                x_p.append(new_x)
                y_p.append(new_y)

                p = 0.4
            '''if np.abs(DR[i])<=dr+eps:
                 x_p.append(planetpositions[1, 0, i])
                 y_p.append(planetpositions[1, 1, i])'''
            '''if DR[i] < DR[i-1] and DR[i] < DR[i+1]:
                x_temp.append(planetpositions[1, 0, i])
                y_temp.append(planetpositions[1, 1, i])'''



    x_p = np.array(x_p)
    y_p = np.array(y_p)
    print ('Number of rounds = ',round, " and length of x is ", len(x_p) )
    angle = np.arctan(y_p[-2]/x_p[-2])*3600
    degrees = np.arctan(y_p[-2]/x_p[-2])
    plt.plot(x_p,y_p, 'o', label = 'Perihelion points')
    plt.legend(loc='best',fontsize = 23)
    plt.xlabel('x-position [AU]',fontsize = 23)
    plt.ylabel('y-position [AU]',fontsize = 23)
    plt.title('Trajectory of Mercury with perihelion points', fontsize = 23)

    print ('Angles at perihelion',angle, 'degrees: ', degrees)
    plt.show()

dt4 = 1e-4
dt5 = 1e-6

find_perihelion("Mercury_perihelion_relativistic.xyz", 2, int(100/dt5)-1)
plt.show()
