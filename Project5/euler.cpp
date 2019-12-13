#include "euler.h"
#include "solarsystem.h"

Euler::Euler(double dt) :
    m_dt(dt)
{

}

void Euler::integrateOneStep(SolarSystem &system)
{
    system.calculateForcesAndEnergy();

    for(CelestialBody &body : system.bodies()) {
        body.position += body.velocity*m_dt;
        body.velocity += body.force / body.mass * m_dt;
    }
}



void Euler::velocityVerletOneStep(SolarSystem &system){


    double m_dt2 = m_dt*m_dt/2;
    vec3 old_a[system.numberOfBodies()];

    for (int i = 0; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies()[i];
        old_a[i] = body.force/body.mass;
}

    for(CelestialBody &body : system.bodies()) {
        body.position += body.velocity*m_dt +m_dt2*body.force/body.mass;
    }
    system.calculateForcesAndEnergy();

    for(int i = 0; i < system.numberOfBodies(); i++) {
        CelestialBody &body = system.bodies()[i];
        body.velocity += m_dt/2 * (body.force/body.mass+old_a[i]);
         }

}

void Euler::velocityVerletOneStepSunFixed(SolarSystem &system){


    double m_dt2 = m_dt*m_dt/2;
    vec3 old_a[system.numberOfBodies()-1];

    for (int i = 1; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies()[i];
        old_a[i] = body.force/body.mass;
}

    for (int i = 1; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies()[i];
        body.position += body.velocity*m_dt +m_dt2*body.force/body.mass;
    }
    system.calculateForcesAndEnergy();

    for(int i = 1; i < system.numberOfBodies(); i++) {
        CelestialBody &body = system.bodies()[i];
        body.velocity += m_dt/2 * (body.force/body.mass+old_a[i]);
         }

}













