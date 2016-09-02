#ifndef SOURCE_H
#define SOURCE_H

class Source {
    int id;
    float i,j,k;
    unsigned int partsPerHour;
    int mode;
    int start;
    int end;
    public:
        Source(unsigned int,float,float,float,unsigned int);
        float getI();
        float getJ();
        float getK();
        unsigned int getPartsPerHour();
};

#endif
