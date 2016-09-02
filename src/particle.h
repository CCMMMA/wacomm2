#ifndef PARTICLE_H
#define PARTICLE_H

#include "wacomm.h"

class Particle {
    float survprob;
    int id;
    float i,j,k;
    int health; // Status { 0: died , 1: alive }
    int time; // Age 
    float pstatus; // Survive probability

    public:
        Particle(int,float,float,float);
        float getI();
        float getJ();
        float getK();
        int getTime();
        int getHealth(); 
        float getPStatus();
        void move(float3d, float3d, float3d, float3d,double2d, double2d, double2d,float conc,float2d zeta,double2d s_w,double2d s_rho,double2d lon_u,double2d lat_v,double dti, double deltat);
};

#endif
