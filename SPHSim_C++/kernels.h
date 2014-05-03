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

#if 0
inline float poly6(float h, float r) {
    if( r > h ) return 0;
    
    const float h2 = h*h;
    float q = r / h;
    float u = 1 - q*q;
    return 4.0 / M_PI / h2 * u * u * u;
}

inline float spiky_grad(float h, float r) {
    if( r > h ) return 0;

    const float h2 = h*h;
    const float h4 = h2*h2;
    float q = r / h;
    float u = 1 - q;
    
    return 30.0 / M_PI / h4 * u * u / q;
}

inline float viscosity_lap(float h, float r) {
    if( r > h ) return 0;
    
    const float h2 = h*h;
    const float h4 = h2*h2;
    float q = r / h;
    
    return 40.0 / M_PI / h4 * (1-q);
}
#else
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

#endif
