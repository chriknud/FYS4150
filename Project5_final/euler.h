#ifndef EULER_H
#define EULER_H


class Euler
{
public:
    double m_dt;
    Euler(double dt);
    void integrateOneStep(class SolarSystem &system);
    void velocityVerletOneStep(class SolarSystem &system);
    void velocityVerletOneStepSunFixed(class SolarSystem &system);
    void velocityVerletOneStepRelativistic(class SolarSystem &system);
};

#endif // EULER_H
