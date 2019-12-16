#include "solarsystem.h"
#include <iostream>
#include <math.h>
using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            // Calculate the force and potential energy here
            vec3 F = body1.mass*body2.mass*G*deltaRVector/(pow(dr,3));

            body1.force -= F;

            body2.force += F;


            m_potentialEnergy -= 2*body1.mass*body2.mass*G/dr;

        }
        m_angularMomentum += (body1.position).cross(body1.mass*body1.velocity);
        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();

    }
}

void SolarSystem::calculateForcesAndEnergyRelativistic()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];

        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 m_angularMomentum2 = (body2.position).cross(body2.velocity);

            double l = m_angularMomentum2.length();

            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();
            // Calculate the force and potential energy here


            vec3 F = body1.mass*body2.mass*G*deltaRVector/pow(dr, 3) * (1+3*pow(l, 2)/pow(dr*c,2));


            body1.force -= F;

            body2.force += F;


            m_potentialEnergy -= 2*body1.mass*body2.mass*G/dr;
            m_angularMomentum += (body1.position).cross(body1.mass*body1.velocity);

        }

        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}

int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    //m_file << numberOfBodies() << endl;
    //m_file << "Comment line that needs to be here." << endl;
    for(CelestialBody &body : m_bodies) {
        m_file << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n";
    }
}

void SolarSystem::writeToFile_energy_and_momentum(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    //m_file << numberOfBodies() << endl;
    //m_file << "Comment line that needs to be here." << endl;
    double momentum = m_angularMomentum.length();

    m_file << momentum << " " << m_potentialEnergy +m_kineticEnergy << "\n";
}


vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}
