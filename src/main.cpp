#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <netcdf>

#include <boost/multi_array.hpp>

#include "wacomm.h"
#include "json/json.h"
#include "source.h"
#include "particle.h"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


// Return this in event of a problem.
static const int NC_ERR = 2;

int main() {
    double dti=30;

    vector<Source*> sources;
    vector<Particle*> particles;

    std::ifstream file("conf.json");
    std::stringstream buffer;

    buffer << file.rdbuf();
    std::string config_doc = buffer.str();

    Json::Value root;   // will contains the root value after parsing.
    Json::Reader reader;
    bool parsingSuccessful = reader.parse( config_doc, root );
    if ( !parsingSuccessful ) {
        // report to the user the failure and their locations in the document.
        std::cout  << "Failed to parse configuration\n" << reader.getFormattedErrorMessages();
        return -1;
    }

    // Get the value of the member of root named 'encoding', return 'UTF-8' if there is no
    // such member.
    std::string nc_input = root.get("nc_input", "" ).asString();
    // Get the value of the member of root named 'encoding', return a 'null' value if
    // there is no such member.
    const Json::Value sources1 = root["sources"];
    for ( unsigned int index = 0; index < sources1.size(); ++index )  { // Iterates over the sequence elements.
        int id=sources1[index].get("id",NULL).asInt();
        float i=sources1[index].get("i",NULL).asFloat();
        float j=sources1[index].get("j",NULL).asFloat();
        float k=sources1[index].get("k",NULL).asFloat();
        unsigned int partsPerHour=sources1[index].get("partsPerHour",NULL).asFloat();
        Source *source=new Source(id,i,j,k,partsPerHour);
        sources.push_back(source);
        printf ("Source: %d\n",id);
    }
   


    // Open the file for read access
    NcFile dataFile(nc_input, NcFile::read);
   

    // Retrieve the variable named "data"
    unsigned int ocean_time = dataFile.getDim("ocean_time").getSize();
    int s_w        = dataFile.getDim("s_w").getSize();
    int s_rho      = dataFile.getDim("s_rho").getSize();
    int eta_u      = dataFile.getDim("eta_u").getSize();
    int xi_u       = dataFile.getDim("xi_u").getSize();
    int eta_v      = dataFile.getDim("eta_v").getSize();
    int xi_v       = dataFile.getDim("xi_v").getSize();
    int eta_rho    = dataFile.getDim("eta_rho").getSize();
    int xi_rho     = dataFile.getDim("xi_rho").getSize();

    double *_ocean_time=new double[ocean_time];
    double *_lat_u=new double[eta_u*xi_u];
    double *_lon_u=new double[eta_u*xi_u];
    double *_lat_v=new double[eta_v*xi_v];
    double *_lon_v=new double[eta_v*xi_v];
    double *_lat_rho=new double[eta_rho*xi_rho];
    double *_lon_rho=new double[eta_rho*xi_rho];

    NcVar dataOceanTime=dataFile.getVar("ocean_time");
    dataOceanTime.getVar(_ocean_time);
    double1d oceanTime(_ocean_time, boost::extents[ocean_time]);

    NcVar dataLatU=dataFile.getVar("lat_u");
    dataLatU.getVar(_lat_u);
    double2d lat_u(_lat_u, boost::extents[eta_u][xi_u]);

    NcVar dataLonU=dataFile.getVar("lon_u");
    dataLonU.getVar(_lon_u);
    double2d lon_u(_lon_u, boost::extents[eta_u][xi_u]);

    NcVar dataLatV=dataFile.getVar("lat_v");
    dataLatV.getVar(_lat_v);
    double2d lat_v(_lat_v, boost::extents[eta_v][xi_v]);

    NcVar dataLonV=dataFile.getVar("lon_v");
    dataLonV.getVar(_lon_v);
    double2d lon_v(_lon_v, boost::extents[eta_v][xi_v]);

    NcVar dataLatRho=dataFile.getVar("lat_rho");
    dataLatRho.getVar(_lat_rho);
    double2d lat_rho(_lat_rho, boost::extents[eta_rho][xi_rho]);

    NcVar dataLonRho=dataFile.getVar("lon_rho");
    dataLonRho.getVar(_lon_rho);
    double2d lon_rho(_lon_rho, boost::extents[eta_rho][xi_rho]);

    printf ("jj\n");

    double *_mask_u=new double[eta_u*xi_u];
    double *_mask_v=new double[eta_v*xi_v];
    double *_mask_rho=new double[eta_rho*xi_rho];

    NcVar dataMaskU=dataFile.getVar("mask_u");
    NcVar dataMaskV=dataFile.getVar("mask_v");
    NcVar dataMaskRho=dataFile.getVar("mask_rho");

    dataMaskU.getVar(_mask_u);
    dataMaskV.getVar(_mask_v);
    dataMaskRho.getVar(_mask_rho);

    double2d mask_u(_mask_u, boost::extents[eta_u][xi_u]);
    double2d mask_v(_mask_v, boost::extents[eta_v][xi_v]);
    double2d mask_rho(_mask_rho, boost::extents[eta_rho][xi_rho]);

    NcVar dataH=dataFile.getVar("h");
    double *_h=new double[eta_rho*xi_rho];
    dataH.getVar(_h);
    double2d h(_h,boost::extents[eta_rho][xi_rho]);

    NcVar dataZeta=dataFile.getVar("zeta");
    NcVar dataU=dataFile.getVar("u");
    NcVar dataV=dataFile.getVar("v");
    NcVar dataW=dataFile.getVar("w");
    NcVar dataAKt=dataFile.getVar("AKt");



    float *_zeta=new float[eta_rho*xi_rho];
    float *_u=new float[s_rho*eta_u*xi_u];
    float *_v=new float[s_rho*eta_v*xi_v];
    float *_w=new float[s_rho*eta_rho*xi_rho];
    float *_akt=new float[s_w*eta_rho*xi_rho];

    float *_conc=new float[s_rho*eta_rho*xi_rho];
    float2d conc(_conc, boost::extents[s_rho][eta_rho][xi_rho]);

    vector<size_t> countp_zeta,countp_u,countp_v,countp_w,countp_akt;
    countp_zeta.push_back(1);
    countp_zeta.push_back(eta_rho);
    countp_zeta.push_back(xi_rho);

    countp_u.push_back(1);
    countp_u.push_back(s_rho);
    countp_u.push_back(eta_u);
    countp_u.push_back(xi_u);
    countp_v.push_back(1);
    countp_v.push_back(s_rho);
    countp_v.push_back(eta_v);
    countp_v.push_back(xi_v);
    countp_w.push_back(1);
    countp_w.push_back(s_rho);
    countp_w.push_back(eta_rho);
    countp_w.push_back(xi_rho);
    countp_akt.push_back(1);
    countp_akt.push_back(s_w);
    countp_akt.push_back(eta_rho);
    countp_akt.push_back(xi_rho);

    unsigned int particleCount=1;

    double deltat=oceanTime[1]-oceanTime[0];

    for (unsigned int t=0;t<ocean_time;t++) {
        printf ("Time: %d/%d\n",(t+1),ocean_time);

        // Add the new particles
        for(unsigned int s=0;s<sources.size();s++) {
            Source *source=sources[s];
            for (unsigned int p=0;p<source->getPartsPerHour();p++) {
                Particle *particle=new Particle(particleCount,source->getI(),source->getJ(),source->getK());
                particles.push_back(particle);
                particleCount++;
            }
        }

        vector<size_t> startp_zeta,startp;
        startp_zeta.push_back(t);
        startp_zeta.push_back(0);
        startp_zeta.push_back(0);
        startp.push_back(t);
        startp.push_back(0);
        startp.push_back(0);
        startp.push_back(0);
        printf("Load\n");
        dataZeta.getVar(startp_zeta,countp_zeta,_zeta);
        dataU.getVar(startp,countp_u,_u);
        dataV.getVar(startp,countp_v,_v);
        dataW.getVar(startp,countp_w,_w);
        dataAKt.getVar(startp,countp_akt,_akt);

        float2d zeta(_zeta, boost::extents[eta_rho][xi_rho]);
        float3d u(_u, boost::extents[s_rho][eta_u][xi_u]);
        float3d v(_v, boost::extents[s_rho][eta_v][xi_v]);
        float3d w(_w, boost::extents[s_rho][eta_rho][xi_rho]);
        float3d akt(_akt, boost::extents[s_w][eta_rho][xi_rho]);

        printf ("Calc\n");

        for (unsigned int p=0;p<particles.size();p++) {
            Particle *particle=particles[p];
            particle->move(u,v,w,akt,mask_u,mask_v,mask_rho,conc,zeta,s_w,s_rho,lon_u,lat_v,dti,deltat);
        }

        // Check the values. 
        //for (int j= 0; j < eta_u; j++)
        //    for (int i = 0; i < xi_u; i++)
        //        printf("%2.4f\n",ut[k*eta_u*xi_u+j*eta_u+i]);
   
    }

    printf ("Latest particle id: %d\n",particleCount);

    delete _akt;
    delete _conc;
    delete _w;
    delete _v;
    delete _u;
    delete _h;
    delete _zeta;
    delete _mask_u;
    delete _mask_v;
    delete _mask_rho;
    delete _lat_u;
    delete _lon_u;
    delete _lat_v;
    delete _lon_v;
    delete _lat_rho;
    delete _lon_rho;
    delete _ocean_time;
    return 0;
}
