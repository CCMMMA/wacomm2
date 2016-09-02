#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

Particle::Particle(int id, float i, float j, float k) {
    this->id=id;
    this->i=i;
    this->j=j;
    this->k=k;
    health=1;
    time=0;
    pstatus=1.0f;
    survprob=1.0f;
};

float Particle::getI() { return i; }
float Particle::getJ() { return j; }
float Particle::getK() { return k; }
int Particle::getTime() { return time; }
int Particle::getHealth() { return health; }
float Particle::getPStatus() { return pstatus; }

void Particle::move(float3d u, float3d v, float3d w, float3d akt,double2d mask_u, double2d mask_v, double2d mask_rho, float conc,float2d zeta,double2d s_w,double2d s_rho,double2d lon_u,double2d lat_v,double dti, double deltat) {
    boost::minstd_rand intgen;
    boost::uniform_01<boost::minstd_rand> gen(intgen);

    unsigned int wdims=u.shape()[0];
    float idet=0,jdet=0,kdet=0;
    float crid=1;

    float vs=0;
    double iint=deltat/dti;
    double t=0;
    printf ("Steps:%f\n",iint);
    for (unsigned t=0;t<iint;t++) {
        if (health<survprob) {
            pstatus=0;
            continue;
        }
        int iI=(int)i; float iF=i-iI;
        int jI=(int)j; float jF=j-jI;
        int kI=(int)k; float kF=k-kI;

        if (mask_rho[iI][jI]<=0.0 || j<0 || i<0 || j>=mask_u.shape()[0]||i<mask_u.shape()[1]) {
            health=-1;
            pstatus=0;
            continue;
        }

        float u1=u[kF][jF][iF]*(1-iF)*(1-jF);
        float u2=u[kF][jF+1][iF]*(1-iF)*jF;
        float u3=u[kF][jF+1][iF+1]*iF*jF;
        float u4=u[kF][jF][iF+1]*iF*(1-jF);
        float uu=u1+u2+u3+u4;

        float v1=v[kF][jF][iF]*(1-iF)*(1-jF);
        float v2=v[kF][jF+1][iF]*(1-iF)*jF;
        float v3=v[kF][jF+1][iF+1]*iF*jF;
        float v4=v[kF][jF][iF+1]*iF*(1-jF);
        float vv=v1+v2+v3+v4;

        float w1=w[kF]  [jF]  [iF]*(1-iF)*(1-jF)*(1-kF);
        float w2=w[kF][jF+1]  [iF]*(1-iF)*    jF*(1-kF);
        float w3=w[kF][jF+1][iF+1]*    iF*    jF*(1-kF);
        float w4=w[kF]  [jF][iF+1]*    iF*(1-jF)*(1-kF);
        float w5=w[kF+1]  [jF]  [iF]*(1-iF)*(1-jF)*kF;
        float w6=w[kF+1][jF+1]  [iF]*(1-iF)*    jF*kF;
        float w7=w[kF+1][jF+1][iF+1]*    iF*    jF*kF;
        float w8=w[kF+1]  [jF][iF+1]*    iF*(1-jF)*kF;
        float ww=w1+w2+w3+w4+w5+w6+w7+w8;

        float ileap=uu*dti;
        float jleap=vv*dti;
        float kleap=(vs+ww)*dti;

        double sigmaprof=3.46*(1+kdet/wdims);
        double gi=0,gj=0,gk=0;
        for (int i=0;i<12;i++) {
            gi=gi+gen()-0.5;
            gj=gj+gen()-0.5;
            gk=gk+gen()-0.5;
        }
        float rileap=gi*sigmaprof;
        float rjleap=gj*sigmaprof;

        float a1=akt[kF]  [jF]  [iF]*(1-iF)*(1-jF)*(1-kF);
        float a2=akt[kF][jF+1]  [iF]*(1-iF)*    jF*(1-kF);
        float a3=akt[kF][jF+1][iF+1]*    iF*    jF*(1-kF);
        float a4=akt[kF]  [jF][iF+1]*    iF*(1-jF)*(1-kF);
        float a5=akt[kF+1]  [jF]  [iF]*(1-iF)*(1-jF)*kF;
        float a6=akt[kF+1][jF+1]  [iF]*(1-iF)*    jF*kF;
        float a7=akt[kF+1][jF+1][iF+1]*    iF*    jF*kF;
        float a8=akt[kF+1]  [jF][iF+1]*    iF*(1-jF)*kF;
        float aa=a1+a2+a3+a4+a5+a6+a7+a8;

        float rkleap=gk*aa*crid;

        ileap=ileap+rileap;
        jleap=jleap+rjleap;
        kleap=kleap+rkleap;

        double d1=(lat_v[jI+1][iI]-lat_v[jI][iI]);
        double d2=(lon_u[jI][iI+1]-lon_u[jI][iI]);
        d1=pow(sin(0.5*d1),2) + pow(sin(0.5*d2),2)* cos(lat_v[jI+1][iI])*cos(lat_v[jI][iI]);
        double dist=2.0*atan2(pow(d1,.5),pow(1.0-d1,.5))*6371.0;
        idet=i+0.001*ileap/dist


    }
}
