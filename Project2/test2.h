#ifndef TEST2_H
#define TEST2_H
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "armadillo"
using namespace std;
using namespace arma;

class test2
{
public:
    test2();
    int value;
    int func(int value);
};

#endif // TEST2_H
