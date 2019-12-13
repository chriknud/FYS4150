#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H


#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>


class SolarSystem
{
public:
    SolarSystem();
    CelestialBody &createCelestialBody(vec3 position, vec3 velocity, double mass);
    void calculateForcesAndEnergy();
    int numberOfBodies() const;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);
    vec3 angularMomentum() const;
    std::vector<CelestialBody> &bodies();
    double G = (4 * M_PI * M_PI);
    void calculateForcesAndEnergyRelativistic();
private:
    std::vector<CelestialBody> m_bodies;
    vec3 m_angularMomentum;
    std::ofstream m_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
};

#endif // SOLARSYSTEM_H
