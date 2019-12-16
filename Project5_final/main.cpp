#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.h"
#include "euler.h"
using namespace std;

void print_output(class SolarSystem &solarSystem);
void personal_input();
void NASA_values(class SolarSystem &solarSystem);
void Earth_values(class SolarSystem &solarSystem);



int main(int numArguments, char **arguments)
{
    double years = 30;
    double dt = 1e-5;
    int numTimesteps = years/dt;


    SolarSystem solarSystem;

    CelestialBody &sun = solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );
    NASA_values(solarSystem);
    //Earth_values(solarSystem);
    double Msun = 2e30;
    double MEarth = 6e24/Msun;
    double MJupiter = 1.9e27/Msun;
    double MMercury = 3.3e23/Msun;

    double xEarth=3.948527228009325E-01, yEarth=9.100160380472437E-01, zEarth=-2.709495540997714E-05;
    double vxEarth=-1.603566066362973E-02, vyEarth=6.880628606437826E-03, vzEarth=6.059839199455330E-07;
    double xJupiter=2.771209156933313E-01, yJupiter=-5.224508231691265, zJupiter=1.546777941340911E-02;
    double vxJupiter=7.443826199722129E-03, vyJupiter=7.587138148929383E-04, vzJupiter=-1.696810240208194E-04;

    //Mercury
    //solarSystem.createCelestialBody(vec3(-0.3075, 0, 0),vec3(0,-12.44, 0),MMercury );
    //Earth
    //solarSystem.createCelestialBody( vec3(1,0,0), vec3(0, 2*M_PI, 0), MEarth);
    //solarSystem.createCelestialBody( vec3(xEarth,yEarth,zEarth), vec3(vxEarth, vyEarth, vzEarth)*365.25, MEarth);

    //Jupiter
    //solarSystem.createCelestialBody( vec3(5.2,0,0), vec3(0,0.548*M_PI,0), MJupiter);
    //solarSystem.createCelestialBody( vec3(xJupiter,yJupiter,zJupiter), vec3(0.9*vxJupiter, vyJupiter, vzJupiter)*365.25, MJupiter);


    vector<CelestialBody> &bodies = solarSystem.bodies();

    double total_mass = 0;
    vec3 R_cm, V_cm;
    for(int i=0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i];
        total_mass += body.mass;
        R_cm += body.position * body.mass;
        V_cm += body.velocity * body.mass;
    }

    R_cm /= total_mass;
    V_cm /= total_mass;
    for(int i = 0; i<bodies.size(); i++) {
      CelestialBody &body = bodies[i]; // Reference to this body
      body.position += R_cm;
        if (i==0){
            body.velocity += V_cm;
        }
        cout << "R_cm = " << R_cm << ", V_cm = " << V_cm << endl;
       cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;

    }

    Euler integrator(dt);
    solarSystem.calculateForcesAndEnergy();
    //solarSystem.calculateForcesAndEnergyRelativistic();
    print_output(solarSystem);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        //integrator.integrateOneStep(solarSystem);
        integrator.velocityVerletOneStep(solarSystem);
        //integrator.velocityVerletOneStepRelativistic(solarSystem);

        solarSystem.writeToFile("SolarSystem.xyz");
        //solarSystem.writeToFile("Mercury_perihelion_relativistic.xyz");
        //solarSystem.writeToFile_energy_and_momentum("energy_and_mome//ntum_verlet.xyz");
    }
    print_output(solarSystem);
    cout << "I just created my first solar system that has " << solarSystem.bodies().size() << " objects." << endl;

    return 0;
}

void Earth_values(SolarSystem &solarSystem){

    double Msun = 2e30;
    double MEarth = 6e24/Msun;
    double xEarth=3.948527228009325E-01, yEarth=9.100160380472437E-01, zEarth=-2.709495540997714E-05;
    double vxEarth=-1.603566066362973E-02, vyEarth=6.880628606437826E-03, vzEarth=6.059839199455330E-07;

    //Earth
    solarSystem.createCelestialBody( vec3(xEarth,yEarth,zEarth), vec3(vxEarth, vyEarth, vzEarth)*365.25, MEarth);


}

void NASA_values(SolarSystem &solarSystem){

    double Msun = 2e30;

    double MEarth = 6e24/Msun, MJupiter = 1.9e27/Msun, MMars = 6.6e23/Msun, MVenus = 4.9e24/Msun, MSaturn = 5.5e26/Msun;
    double MMercury = 3.3e23/Msun, MUranus = 8.8e25/Msun, MNeptun = 1.03e26/Msun, MPluto = 1.31e22/Msun;

    double xMercury = -3.089137495084154E-01, yMercury = 1.744886373010318E-01, zMercury = 4.167600354497743E-02;
    double vxMercury = -1.928258980107407E-02, vyMercury = -2.350312105925493E-02, vzMercury = -1.520556440066312E-04;

    double xVenus = 4.814605067455450E-01, yVenus = -5.345470370402726E-01, zVenus = -3.540916553726071E-02;
    double vxVenus=1.493272115673404E-02, vyVenus=1.341199462523215E-02, vzVenus=-6.779258843704102E-04;

    double xEarth=3.948527228009325E-01, yEarth=9.100160380472437E-01, zEarth=-2.709495540997714E-05;
    double vxEarth=-1.603566066362973E-02, vyEarth=6.880628606437826E-03, vzEarth=6.059839199455330E-07;

    double xMars=-1.543208932754952, yMars=-5.042083307102040E-01, zMars=2.707009173557715E-02;
    double vxMars=4.926983218616618E-03, vyMars=-1.208451455394788E-02, vzMars=-3.740559366586062E-04;

    double xJupiter=2.771209156933313E-01, yJupiter=-5.224508231691265, zJupiter=1.546777941340911E-02;
    double vxJupiter=7.443826199722129E-03, vyJupiter=7.587138148929383E-04, vzJupiter=-1.696810240208194E-04;

    double xSaturn = 3.632628879585697E+00, ySaturn = -9.348288811274543E+00, zSaturn = 1.793542343960454E-02;
    double vxSaturn = 4.891570166847385E-03, vySaturn = 2.004764720827292E-03, vzSaturn = -2.299017365589604E-04;

    double xUranus = 1.629688404988837E+01, yUranus = 1.128605542338266E+01, zUranus = -1.692115522793217E-01;
    double vxUranus = -2.268171563746997E-03, vyUranus = 3.050128900966884E-03, vzUranus = 4.061969844748497E-05;

    double xNeptun = 2.921750763559268E+01, yNeptun = -6.461552366128481E+00, zNeptun = -5.402832642395710E-01;
    double vxNeptun =6.568727514842341E-04, vyNeptun = 3.084041276878159E-03, vzNeptun = -7.824234958842697E-05;

    double xPluto = 1.287478886188476E+01, yPluto = -3.138014806326925E+01, zPluto = -3.662996467798016E-01;
    double vxPluto = 2.965077404403603E-03, vyPluto = 5.169914189921419E-04, vzPluto = -9.099615783803918E-04;
    //Mercury
    solarSystem.createCelestialBody( vec3(xMercury,yMercury,zMercury), vec3(vxMercury, vyMercury, vzMercury)*365.25, MMercury);

    //Venus
    solarSystem.createCelestialBody( vec3(xVenus,yVenus,zVenus), vec3(vxVenus, vyVenus, vzVenus)*365.25, MVenus);

    //Earth
    solarSystem.createCelestialBody( vec3(xEarth,yEarth,zEarth), vec3(vxEarth, vyEarth, vzEarth)*365.25, MEarth);

    //Mars
    solarSystem.createCelestialBody( vec3(xMars,yMars,zMars), vec3(vxMars, vyMars, vzMars)*365.25, MMars);

    //Jupiter
    solarSystem.createCelestialBody( vec3(xJupiter,yJupiter,zJupiter), vec3(vxJupiter, vyJupiter, vzJupiter)*365.25, MJupiter);

    //Saturn
    solarSystem.createCelestialBody( vec3(xSaturn,ySaturn,zSaturn), vec3(vxSaturn, vySaturn, vzSaturn)*365.25, MSaturn);

    //Uranus
    solarSystem.createCelestialBody( vec3(xUranus,yUranus,zUranus), vec3(vxUranus, vyUranus, vzUranus)*365.25, MUranus);

    //Neptun
    solarSystem.createCelestialBody( vec3(xNeptun,yNeptun,zNeptun), vec3(vxNeptun, vyNeptun, vzNeptun)*365.25, MNeptun);

    //Pluto
    solarSystem.createCelestialBody( vec3(xPluto,yPluto,zPluto), vec3(vxPluto, vyPluto, vzPluto)*365.25, MPluto);

}

void personal_input(){
    SolarSystem solarSystem;
    double Msun = 2e30;
    double MEarth = 6e24/Msun, MJupiter = 1.9e27/Msun, MMars = 6.6e23/Msun, MVenus = 4.9e24/Msun, MSaturn = 5.5e26/Msun;
    double MMercury = 3.3e23/Msun, MUranus = 8.8e25/Msun, MNeptun = 1.03e26/Msun, MPluto = 1.31e22/Msun;

    //Mercury
    solarSystem.createCelestialBody( vec3(-0.39,0, 0), vec3(0, 2*M_PI, 0), MMercury);

    //Venus
    solarSystem.createCelestialBody( vec3(0, -0.72, 0), vec3(2*M_PI, 0, 0), MVenus);

    //Earth
    solarSystem.createCelestialBody( vec3(1, 0, 0), vec3(0, 2*M_PI, 0), MEarth );
    //solarSystem.createCelestialBody( vec3(3.948527228009325e-1, 9.100160380472437e-1, -2.709495540997714e-5), vec3(1.603566066362973e-2,6.880628606437826e-3, 6.059839199455330e-7), 3e-6 );

    //Mars
    solarSystem.createCelestialBody( vec3(0, 1.52, 0), vec3(2*M_PI, 0, 0), MMars);

    //Jupiter
    solarSystem.createCelestialBody( vec3(5.2, 0, 0), vec3(0, 0.6*M_PI, 0), MJupiter*1000);
    //solarSystem.createCelestialBody( vec3(5.2, 0, 0), vec3(0, 0.88*M_PI, 0), MJupiter);
    //solarSystem.createCelestialBody( vec3(2.771209156933313, -5.224508231691265, 1.546777941340911e-2), vec3(7.443826199722129e-1, 7.587138148929383e-4, -1.696810240208194e-4), 1e-3 );

    //Saturn
    solarSystem.createCelestialBody( vec3(-9.54,0, 0), vec3(0, 0.5*M_PI, 0), MSaturn);

    //Uranus
    solarSystem.createCelestialBody( vec3(19.19,0, 0), vec3(0, 0.5*M_PI, 0), MUranus);

    //Neptun
    solarSystem.createCelestialBody( vec3(0,30.6, 0), vec3(0.4*M_PI,0, 0), MNeptun);

    //Pluto
    solarSystem.createCelestialBody( vec3(0,-39.53, 0), vec3(0.3*M_PI,0, 0), MPluto);
}




void print_output(SolarSystem &solarSystem){
    cout << "Total energy after x years " << solarSystem.totalEnergy() << endl;
    cout << "Kinetic energy after x years " << solarSystem.kineticEnergy() << endl;
    cout << "Potential energy after x years " << solarSystem.potentialEnergy() << endl;
    cout << "Angular momentum after x years " << solarSystem.angularMomentum() << endl;

}
