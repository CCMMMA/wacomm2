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



NcFile *createNetCDF(string file_name,double2d mask_rho, double2d lat_rho, double2d lon_rho, double1d s_rho) {
    // Names of things. 
    string  UNITS = "units";
    string  DEGREES_EAST = "degrees_east";
    string DEGREES_NORTH = "degrees_north";
    string XI_RHO_NAME="xi_rho";
    string ETA_RHO_NAME="eta_rho";
    string S_RHO_NAME="s_rho";
    string OCEAN_TIME_NAME="ocean_time";
    string LAT_RHO_NAME = "lat_rho";
    string LON_RHO_NAME ="lon_rho";
    string MASK_RHO_NAME ="mask_rho";
    string CONC_NAME ="conc";


    try {
        NcFile *sfc=new NcFile(file_name, NcFile::replace);

        NcDim oceanTimeDim = sfc->addDim(OCEAN_TIME_NAME);
        NcDim sRhoDim = sfc->addDim(S_RHO_NAME, s_rho.shape()[0]);
        NcDim etaRhoDim = sfc->addDim(ETA_RHO_NAME, lat_rho.shape()[0]);
        NcDim xiRhoDim = sfc->addDim(XI_RHO_NAME, lon_rho.shape()[1]);

        NcVar oceanTimeVar = sfc->addVar(OCEAN_TIME_NAME, ncDouble, oceanTimeDim);

        NcVar sRhoVar = sfc->addVar(S_RHO_NAME, ncDouble, sRhoDim);
        sRhoVar.putVar(s_rho.data());

        vector<NcDim> etaRhoXiRhoDims;
        etaRhoXiRhoDims.push_back(etaRhoDim);
        etaRhoXiRhoDims.push_back(xiRhoDim);
        NcVar latRhoVar = sfc->addVar(LAT_RHO_NAME, ncDouble, etaRhoXiRhoDims);
        NcVar lonRhoVar = sfc->addVar(LON_RHO_NAME, ncDouble, etaRhoXiRhoDims); 
        NcVar maskRhoVar = sfc->addVar(MASK_RHO_NAME, ncDouble, etaRhoXiRhoDims);
    
        latRhoVar.putVar(lat_rho.data());
        lonRhoVar.putVar(lon_rho.data());
        maskRhoVar.putVar(mask_rho.data());

        lonRhoVar.putAtt(UNITS,DEGREES_EAST);
        latRhoVar.putAtt(UNITS,DEGREES_NORTH);

        vector<NcDim> oceanTimeSRhoEtaRhoXiRhoDims;
        oceanTimeSRhoEtaRhoXiRhoDims.push_back(oceanTimeDim);
        oceanTimeSRhoEtaRhoXiRhoDims.push_back(sRhoDim);
        oceanTimeSRhoEtaRhoXiRhoDims.push_back(etaRhoDim);
        oceanTimeSRhoEtaRhoXiRhoDims.push_back(xiRhoDim);
        NcVar concVar = sfc->addVar(CONC_NAME, ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
        return sfc;
    }
    catch(NcException& e) {
        e.what(); 
        return NULL;
    }
}

int main() {
    boost::minstd_rand intgen;
    boost::uniform_01<boost::minstd_rand> gen(intgen);

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
    for ( unsigned int index = 0; index < sources1.size(); index++ )  { // Iterates over the sequence elements.
        int id=sources1[index].get("id",NULL).asInt();
        float i=sources1[index].get("i",NULL).asFloat();
        float j=sources1[index].get("j",NULL).asFloat();
        float k=sources1[index].get("k",NULL).asFloat();
        unsigned int partsPerHour=sources1[index].get("partsPerHour",0).asFloat();
        Source *source=new Source(id,i,j,k,partsPerHour);
        sources.push_back(source);
        printf ("Source: %d\n",id);
    }

    const Json::Value chm = root["chm"];
    float tau0=chm.get("tau0",86400.0).asFloat();
    float survprob=chm.get("survprob",1.0e-4).asFloat();
    printf ("Input:%s\n",nc_input.c_str());

    // Open the file for read access
    NcFile dataFile(nc_input, NcFile::read);
   

    // Retrieve the variable named "data"
    unsigned int ocean_time = dataFile.getDim("ocean_time").getSize();
    unsigned int s_w        = dataFile.getDim("s_w").getSize();
    unsigned int s_rho      = dataFile.getDim("s_rho").getSize();
    unsigned int eta_u      = dataFile.getDim("eta_u").getSize();
    unsigned int xi_u       = dataFile.getDim("xi_u").getSize();
    unsigned int eta_v      = dataFile.getDim("eta_v").getSize();
    unsigned int xi_v       = dataFile.getDim("xi_v").getSize();
    unsigned int eta_rho    = dataFile.getDim("eta_rho").getSize();
    unsigned int xi_rho     = dataFile.getDim("xi_rho").getSize();

    double *_ocean_time=new double[ocean_time];
    double *_s_w=new double[s_w];
    double *_s_rho=new double[s_rho];
    double *_lat_u=new double[eta_u*xi_u];
    double *_lon_u=new double[eta_u*xi_u];
    double *_lat_v=new double[eta_v*xi_v];
    double *_lon_v=new double[eta_v*xi_v];
    double *_lat_rho=new double[eta_rho*xi_rho];
    double *_lon_rho=new double[eta_rho*xi_rho];

    NcVar dataOceanTime=dataFile.getVar("ocean_time");
    dataOceanTime.getVar(_ocean_time);
    double1d oceanTime(_ocean_time, boost::extents[ocean_time]);

    NcVar dataSw=dataFile.getVar("s_w");
    dataSw.getVar(_s_w);
    double1d sW(_s_w, boost::extents[s_w]);

    NcVar dataSRho=dataFile.getVar("s_rho");
    dataSRho.getVar(_s_rho);
    double1d sRho(_s_rho, boost::extents[s_rho]);

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
    float3d conc(_conc, boost::extents[s_rho][eta_rho][xi_rho]);





    NcFile *outputFile=createNetCDF("output.nc4",mask_rho, lat_rho, lon_rho, sRho);
    NcVar outputOceanTimeVar=outputFile->getVar("ocean_time");
    NcVar outputConcVar=outputFile->getVar("conc");
    vector<size_t> startp_ocean_time,startp_conc,countp_ocean_time,countp_conc;
    startp_ocean_time.push_back(0);
    countp_ocean_time.push_back(1);
    startp_conc.push_back(0);
    startp_conc.push_back(0);
    startp_conc.push_back(0);
    startp_conc.push_back(0);
    countp_conc.push_back(1);
    countp_conc.push_back(s_rho);
    countp_conc.push_back(eta_rho);
    countp_conc.push_back(xi_rho);



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

    // lon_u, lat_v in radiants
    for (unsigned int j=0; j<lon_u.shape()[0];j++)
        for (unsigned int i=0;i<lon_u.shape()[1];i++)
            lon_u[j][i]=0.0174533*lon_u[j][i];
    for (unsigned int j=0; j<lat_v.shape()[0];j++)
        for (unsigned int i=0;i<lat_v.shape()[1];i++)
             lat_v[j][i]=0.0174533*lat_v[j][i];

    // Rho interpoted data
    float *_ucomp=new float[s_rho*eta_v*xi_u];
    float *_vcomp=new float[s_rho*eta_v*xi_u];
    float *_wcomp=new float[s_rho*eta_v*xi_u];
    float *_aktcomp=new float[s_rho*eta_v*xi_u];
    float3d ucomp(_ucomp,boost::extents[s_rho][eta_v][xi_u]);
    float3d vcomp(_vcomp,boost::extents[s_rho][eta_v][xi_u]);
    float3d wcomp(_wcomp,boost::extents[s_rho][eta_v][xi_u]);
    float3d aktcomp(_aktcomp,boost::extents[s_rho][eta_v][xi_u]);

    ocean_time=4;
    for (unsigned int t=0;t<ocean_time;t++) {
        printf ("Time: %d/%d -- %f -- Particles:",(t+1),ocean_time,oceanTime[t]);

        // Add the new particles
        for(unsigned int s=0;s<sources.size();s++) {
            Source *source=sources[s];
            for (unsigned int p=0;p<source->getPartsPerHour();p++) {
                double i=source->getI()+gen()*.5-.25;
                double j=source->getJ()+gen()*.5-.25;
                double k=source->getK()+gen()*.5-.25;
                if (k<0) k=0;
                if (j<0) j=0;
                if (i<0) i=0;
                Particle *particle=new Particle(particleCount,i,j,k,tau0,survprob);
                particles.push_back(particle);
                particleCount++;
            }
        }
        printf("%u\n",(unsigned int)particles.size());
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

        
        // Interpolate on internal rho points
        double uw1, uw2;
        printf("Interpolating: U\n");
        for (unsigned int k=0; k<s_rho;k++) {
            for (unsigned int j=0;j<eta_v;j++) {
                for (unsigned int i=0;i<xi_u;i++) {
                    if ( mask_rho[j][i] > 0.0 ) {
                        if ( mask_u[j][i] > 0.0 ) {
                            uw1=u[k][j][i];
                        } else {
                            uw1=0.0;
                        }
                        if ( j>0 && mask_u[j-1][i]> 0.0 ) {
                            uw2=u[k][j-1][i];
                        } else {
                            uw2=0.0;
                        }

                        ucomp[k][j][i]=0.5*(uw1+uw2);
                    } else {  
                        ucomp[k][j][i]=0.0;
                    }
                }
            }
        }

        double vw1,vw2;
        printf("Interpolating V\n");
        for (unsigned int k=0; k<s_rho;k++) {
            for (unsigned int j=0;j<eta_v;j++) {
                for (unsigned int i=0;i<xi_u;i++) {
                    if ( mask_rho[j][i] > 0.0 ) {
                        if ( mask_v[j][i] > 0.0 ) {
                            vw1=v[k][j][i];
                        } else {
                            vw1=0.0;
                        }

                        if ( i>0 && mask_v[j][i-1]> 0.0 ) {
                            vw2=u[k][j][i-1];
                        } else {
                            vw2=0.0;
                        }

                        vcomp[k][j][i]=0.5*(vw1+vw2);
                    } else {
                        vcomp[k][j][i]=0.0;
                    }
                }
            }
	}

        double u1,u2,u3;
        printf("Interpolating W\n");
        for (unsigned int k=0; k<s_rho;k++) {
            for (unsigned int j=0;j<eta_v;j++) {
                for (unsigned int i=0;i<xi_u;i++) {
                    if ( mask_rho[j][i] > 0.0 ) {
                        if ( i>0 && mask_rho[j][i-1] > 0.0 ) { 
                             u1=w[k][j][i-1];
                        } else {
                             u1=0;
                        }

                        if ( j>0 && mask_rho[j-1][i]> 0.0 ) {
                             u2=w[k][j-1][i];
                        } else {
                             u2=0.0;
                        }

                        if ( j>0 && i>0 && mask_rho[j-1][i-1]> 0.0 ) {
                            u3=w[k][j-1][i-1];
                        } else {
                            u3=0.0;
                        }

                        wcomp[k][j][i]=0.25*(u1+u2+u3+w[k][j][i]);
                    } else {
                        wcomp[k][j][i]=0.0;
                    }
                }
            }
	}

        printf("Interpolating AKT\n");
        for (unsigned int k=0; k<s_rho;k++) {
            for (unsigned int j=0;j<eta_v;j++) {
                for (unsigned int i=0;i<xi_u;i++) {
                    if ( mask_rho[j][i] > 0.0 ) {
                        if ( i>0 && mask_rho[j][i-1] > 0.0 ) {
                             u1=akt[k][j][i-1];
                        } else {
                             u1=0;
                        }

                        if ( j>0 && mask_rho[j-1][i]> 0.0 ) {
                             u2=akt[k][j-1][i];
                        } else {
                             u2=0.0;
                        }

                        if ( j>0 && i>0 && mask_rho[j-1][i-1]> 0.0 ) {
                            u3=akt[k][j-1][i-1];
                        } else {
                            u3=0.0;
                        }

                        aktcomp[k][j][i]=0.25*(u1+u2+u3+akt[k][j][i]);
                    } else {
                        aktcomp[k][j][i]=0.0;
                    }
                }
            }
	}

        double *_depth=new double[s_rho+1];
        for (unsigned int k=0;k<=s_rho;k++) {
            //printf("k:%u\n",k);
            if (k<s_rho)
            _depth[k]=sW[k+1]-sW[k];
            else _depth[k]=_depth[k-1];
            //printf("_depth[%u]=%f\n",k,_depth[k]);
        }
        double1d depth(_depth,boost::extents[s_rho+1]);
        
        printf ("Calc\n");
        // The particles could be distributed with mpi
        for (unsigned int p=0;p<particles.size();p++) {
            Particle *particle=particles[p];
            particle->move(ucomp,vcomp,wcomp,aktcomp,mask_u,mask_v,mask_rho,zeta,depth,lon_u,lat_v,dti,deltat);
        }

        printf ("Count:");
        vector<Particle*>::iterator it = particles.begin();

        while(it != particles.end()) {
            Particle *particle=*it;
            if (particle->count(mask_rho,conc)==false) {
                delete particle;
                it = particles.erase(it);
            }
            else ++it;
        }
        startp_ocean_time[0]=t;
        outputOceanTimeVar.putVar(startp_ocean_time,countp_ocean_time,&oceanTime[t]);
        startp_conc[0]=t;
        outputConcVar.putVar(startp_conc,countp_conc,conc.data());
        printf ("%u\n",(unsigned int)particles.size());
    }
    outputFile->close();
    printf ("Latest particle id: %d\n",particleCount);

    std::ofstream myfile;
    myfile.open ("restart.csv");
    myfile << "id;k;j;i;health" << endl;
    for (unsigned int p=0;p<particles.size();p++) {
        Particle *particle=particles[p];
        particle->write(myfile);
    }
    myfile.close();

    delete _wcomp;
    delete _vcomp;
    delete _ucomp;
    delete _aktcomp;
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
