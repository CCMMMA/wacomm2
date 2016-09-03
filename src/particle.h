#ifndef PARTICLE_H
#define PARTICLE_H

#include "wacomm.h"

#include <iostream>
#include <fstream>

using namespace std;

class Particle {
    float survprob;
    int id;
    float i,j,k;
    float health; // Status { 0: died , 1: alive }
    int time; // Age 
    float pstatus; // Survive probability

    float health0;
    float tau0;

    public:
        Particle(int,float,float,float,float,float);
        float getI();
        float getJ();
        float getK();
        int getTime();
        float getHealth(); 
        float getPStatus();
        void write(std::ofstream&);
        void move(float3d, float3d, float3d, float3d,double2d, double2d, double2d,float2d zeta,double1d depth,double2d lon_u,double2d lat_v,double dti, double deltat);
        bool count(double2d, float3d&);
};

#endif
