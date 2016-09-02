#include "source.h"

Source::Source(unsigned int id, float i, float j, float k, unsigned int partsPerHour) {
    this->id=id;
    this->i=i;
    this->j=j;
    this->k=k;
    this->partsPerHour=partsPerHour;
    mode=1;
    start=0;
    end=-1;
}

float Source::getI() { return i; }
float Source::getJ() { return j; }
float Source::getK() { return k; }
unsigned int Source::getPartsPerHour() { return partsPerHour; }
