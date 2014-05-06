//
//  kernels.h
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/26/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#ifndef SPHSim_C___kernels_h
#define SPHSim_C___kernels_h

#include <math.h>

class Kernel {
public:
    Kernel(){}
    Kernel(float h):h(h){
        invh = 1.0 / h;
        h2 = h*h;
        h3 = h2*h;
        h4 = h2*h2;
        h5 = h2*h3;
        h6 = h3*h3;
        h7 = h3*h4;
        h8 = h4*h4;
        h9 = h4*h5;
        
        Cpoly6 = 315.0/64.0*M_PI/h3;
        Cspiky = 45.0/M_PI/h5;
        Cvis = 45.0/M_PI/h5;
    }
    
    float poly6(float r) const {
        if( r > h ) return 0;
        float q = r / h;
        float u = 1 - q;
        
        return Cpoly6 * u * u * u;
    }
    
    float spiky_grad(float r) const {
        if( r > h ) return 0;
        float q = r * invh;
        float u = 1 - q * q;
        return Cspiky * u * u / q;
    }
    
    float viscosity_lap(float r) const {
        if( r > h ) return 0;
        float q = r * invh;
        return Cvis * (1 - q);
    }
    
private:
    double Cpoly6, Cspiky, Cvis;
    double invh;
    double h, h2, h3, h4, h5, h6, h7, h8, h9;
};


inline float poly6(float h, float r) {
    if( r > h ) return 0;
    
    const float h2 = h*h;
    const float h3 = h2*h;
    float q = r / h;
    float u = 1 - q;

    return 315.0 / 64.0 / M_PI / h3 * u * u * u;
}

inline float spiky_grad(float h, float r) {
    if( r > h ) return 0;
    
    const float h2 = h*h;
    const float h5 = h2*h2*h;
    float q = r / h;
    float u = 1 - q * q;
    
    return 45.0 / M_PI / h5 * u * u / q;
}

inline float viscosity_lap(float h, float r) {
    if( r > h ) return 0;
    
    const float h2 = h*h;
    const float h5 = h2*h2*h;
    float q = r / h;
    
    return 45.0 / M_PI / h5 * (1 - q);
}

#endif
