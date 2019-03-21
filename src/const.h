
#ifndef _CONST_H_
#define _CONST_H_

typedef double rp;

const int ord = 5;        // Order of accuracy

const int hs = (ord-1)/2; // Halo size

const int tord = 3;       // Temporal order of accuracy

const int numVars = 3;    // Number of state variables

const int ID_H  = 0;      // Index for height
const int ID_HU = 1;      // Index for height * u-velocity
const int ID_HV = 2;      // Index for height * v-velocity:w


#endif
