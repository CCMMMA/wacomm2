#include "particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>

using namespace std;

double mod(double a, double p) { return a-p*(int)(a/p); }
double sign(double a, double b) { return abs(a)*sgn(b); }

Particle::Particle(int id, float i, float j, float k, float tau0, float survprob) {
    health0=1.0f;
    this->tau0=tau0;

    this->id=id;
    this->i=i;
    this->j=j;
    this->k=k;
    health=1.0f;
    time=0;
    pstatus=1.0f;
    this->survprob=survprob;


};

float Particle::getI() { return i; }
float Particle::getJ() { return j; }
float Particle::getK() { return k; }
int Particle::getTime() { return time; }
float Particle::getHealth() { return health; }
float Particle::getPStatus() { return pstatus; }

void Particle::write(std::ofstream& myfile) {
    if (health>0) {
        myfile << id << ";" << k << ";" << j << ";" << i << ";" << health << endl;;
    }
}

bool Particle::count(double2d mask_rho, float3d& conc) {
    if(health > survprob) {
        unsigned int iI=(int unsigned)i; 
        unsigned int jI=(int unsigned)j; 
        unsigned int kI=(int unsigned)k; 

 	if ( mask_rho[jI][iI] > 0.0 ) {
 	    conc[kI][jI][iI]=conc[kI][jI][iI]+1.0;
            //printf("%u,%u,%u:%f\n",kI,jI,iI,conc[kI][jI][iI]);
 	} else {
 	    health=-2.0;
 	}
    } else {
        health=-2.0;
    }
    if (health<=0) { 
        return false;
    } else {
        return true;
    }
}

void Particle::move(float3d u, float3d v, float3d w, float3d akt,double2d mask_u, double2d mask_v, double2d mask_rho,float2d zeta,double1d depth,double2d lon_u,double2d lat_v,double dti, double deltat) {

    boost::minstd_rand intgen;
    boost::uniform_01<boost::minstd_rand> gen(intgen);

    float idet=0,jdet=0,kdet=0;
    float crid=1;
    double vs=0;

    double iint=deltat/dti;
    double t=0;
    
    for (unsigned t=0;t<iint;t++) {
        if (health<survprob) {
            pstatus=0;
            return;
        }
        unsigned int iI=(unsigned int)i; float iF=i-iI;
        unsigned int jI=(unsigned int)j; float jF=j-jI;
        unsigned int kI=(unsigned int)k; float kF=k-kI;

        if (mask_rho[iI][jI]<=0.0 || jI>=mask_rho.shape()[1]|| iI>=mask_rho.shape()[2]) {
            health=-1;
            pstatus=0;
            return;
        }

        if (k>=w.shape()[0]) {
            health=-13;
            pstatus=0;
            return;
        }

        //printf("Alive:%u kI:%u jI:%u iI:%u kF:%f jF:%f iF:%f\n",id,kI,jI,iI,kF,jF,iF);

        float u1=u[kF][jF][iF]*(1-iF)*(1-jF);
        float u2=u[kF][jF+1][iF]*(1-iF)*jF;
        float u3=u[kF][jF+1][iF+1]*iF*jF;
        float u4=u[kF][jF][iF+1]*iF*(1-jF);
        float uu=u1+u2+u3+u4;

        //printf("uu:%f\n",uu);

        float v1=v[kF][jF][iF]*(1-iF)*(1-jF);
        float v2=v[kF][jF+1][iF]*(1-iF)*jF;
        float v3=v[kF][jF+1][iF+1]*iF*jF;
        float v4=v[kF][jF][iF+1]*iF*(1-jF);
        float vv=v1+v2+v3+v4;

        //printf("vv:%f\n",vv);

        float w1=w[kF]  [jF]  [iF]*(1-iF)*(1-jF)*(1-kF);
        float w2=w[kF][jF+1]  [iF]*(1-iF)*    jF*(1-kF);
        float w3=w[kF][jF+1][iF+1]*    iF*    jF*(1-kF);
        float w4=w[kF]  [jF][iF+1]*    iF*(1-jF)*(1-kF);
        float w5=w[kF+1]  [jF]  [iF]*(1-iF)*(1-jF)*kF;
        float w6=w[kF+1][jF+1]  [iF]*(1-iF)*    jF*kF;
        float w7=w[kF+1][jF+1][iF+1]*    iF*    jF*kF;
        float w8=w[kF+1]  [jF][iF+1]*    iF*(1-jF)*kF;
        float ww=w1+w2+w3+w4+w5+w6+w7+w8;

        //printf("ww:%f\n",ww);

        float ileap=uu*dti;
        float jleap=vv*dti;
        float kleap=(vs+ww)*dti;

        //printf ("kleap:%f jleap:%f ileap:%f\n",kleap,jleap,ileap); 

        double sigmaprof=3.46*(1+kdet/w.shape()[0]);
        double gi=0,gj=0,gk=0;
        for (int a=0;a<12;a++) {
            gi=gi+gen()-0.5;
            gj=gj+gen()-0.5;
            gk=gk+gen()-0.5;
        }

        float a1=akt[kF]  [jF]  [iF]*(1-iF)*(1-jF)*(1-kF);
        float a2=akt[kF][jF+1]  [iF]*(1-iF)*    jF*(1-kF);
        float a3=akt[kF][jF+1][iF+1]*    iF*    jF*(1-kF);
        float a4=akt[kF]  [jF][iF+1]*    iF*(1-jF)*(1-kF);
        float a5=akt[kF+1]  [jF]  [iF]*(1-iF)*(1-jF)*kF;
        float a6=akt[kF+1][jF+1]  [iF]*(1-iF)*    jF*kF;
        float a7=akt[kF+1][jF+1][iF+1]*    iF*    jF*kF;
        float a8=akt[kF+1]  [jF][iF+1]*    iF*(1-jF)*kF;
        float aa=a1+a2+a3+a4+a5+a6+a7+a8;
        //printf("aa:%f\n",aa);

        float rileap=gi*sigmaprof;
        float rjleap=gj*sigmaprof;
        float rkleap=gk*aa*crid;

        //printf ("rkleap:%f rjleap:%f rileap:%f\n",rkleap,rjleap,rileap);

        ileap=ileap+rileap;
        jleap=jleap+rjleap;
        kleap=kleap+rkleap;

        double d1,d2,dist;

        d1=(lat_v[jI+1][iI]-lat_v[jI][iI]);
        d2=(lon_u[jI][iI+1]-lon_u[jI][iI]);
        d1=pow(sin(0.5*d1),2) + pow(sin(0.5*d2),2)* cos(lat_v[jI+1][iI])*cos(lat_v[jI][iI]);
        dist=2.0*atan2(pow(d1,.5),pow(1.0-d1,.5))*6371.0;
        idet=i+0.001*ileap/dist;

        d1=(lat_v[jI+1][iI]-lat_v[jI][iI]);
 	d2=(lon_u[jI][iI+1]-lon_u[jI][iI]);
 	d1=pow(sin(0.5*d1),2) + pow(sin(0.5*d2),2)* cos(lat_v[jI+1][iI])*cos(lat_v[jI][iI]);
 	dist=2.0*atan2(sqrt(d1),pow(1.0-d1,.5))*6371.0;
 	jdet=j+0.001*jleap/dist;

        //printf("depth:%f, zeta:%f\n",depth[kI],zeta[jI][iI]);
 	dist=depth[kI]*zeta[jI][iI];
        if (dist<=0) dist=1e-4;
        //printf("dist:%f, kleap:%f\n",dist,kleap);
 	if ( abs(kleap) > abs(dist) ) {
            kleap=sign(dist,kleap);
        }
 	kdet=k+kleap/dist;

        //printf("kdet:%f jdet:%f idet:%f\n",kdet,jdet,idet);

        if ( kdet >= w.shape()[0] ) {
            kdet=w.shape()[0]-1;
        }
 	if ( kdet < 0. ) {
            kdet=0;
        }

        //printf("kdet:%f jdet:%f idet:%f\n",kdet,jdet,idet);

        unsigned int jdetI=(unsigned int)jdet;
        unsigned int idetI=(unsigned int)idet;
        if ( mask_rho[jdetI][idetI] <= 0.0 ) {
 	    if ( idetI < iI ) {
 	        idet=(float)iI + abs(i-idet);
 	    } else if ( idetI > iI ) {
 		idet=(float)idetI- mod(idet,1.0);
 	    }
 	    if ( jdetI < jdet ) {
 	        jdet=(float)jdetI+ abs(j-jdet);
 	    } else if ( jdetI > jI ) {
 		jdet=(float)jdetI - mod(jdet,1.0);
 	    }
 	}

        i=idet;
        j=jdet;
        k=kdet;

        time=time+dti;
        health=health0*exp(-time/tau0);
    }
}

