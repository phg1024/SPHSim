//
//  SPHSystem.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include "SPHSystem.h"

SPHSystem::SPHSystem() {
    initParticles();
}

SPHSystem::SPHSystem(const string& filename) {
    params = Parameters(filename);
    initParticles();
}

void SPHSystem::writeHeader(ofstream& f) {
    float scale = 1.0;

    f.write(reinterpret_cast<const char*>(&n), sizeof(int));
    f.write(reinterpret_cast<const char*>(&params.nframes), sizeof(int));
    f.write(reinterpret_cast<const char*>(&scale), sizeof(float));
}

void SPHSystem::writeFrame(ofstream& f, float* c = 0) {
    for (int i = 0; i < n; ++i) {
        Vec x = p[i];
        float ci = c ? c[i] : 0;
        f.write(reinterpret_cast<const char*>(&x.x), sizeof(float));
        f.write(reinterpret_cast<const char*>(&x.y), sizeof(float));
        f.write(reinterpret_cast<const char*>(&ci), sizeof(float));
    }
}

void SPHSystem::run() {
    ofstream f(params.logfile, ios::binary);
    
    writeHeader(f);
    writeFrame(f);
    
    computeAcceleration();
    leapFrogStart();
    checkState();
    
    for (int frame = 1; frame < params.nframes; ++frame) {
        for (int i = 0; i < params.nsteps; ++i) {
            computeAcceleration();
            leapFrog();
            checkState();
        }
        writeFrame(f);
    }
}

void SPHSystem::resize(int count) {
    n = count;
    p.resize(n);
    v.resize(n);
    vh.resize(n);
    a.resize(n);
    density.resize(n);
    pressure.resize(n);
}

void SPHSystem::initParticles() {
    float h  = params.h;
    float hh = h/1.3;
    
    auto indicatef = [](float x, float y){ return x<0.5 && y<0.5; };
    
    // Count mesh points that fall in indicated region.
    int count = 0;
    for (float x = 0; x < 1; x += hh)
        for (float y = 0; y < 1; y += hh)
            count += indicatef(x,y);
    
    // Populate the particle data structure
    resize(count);
    
    int i = 0;
    for (float x = 0; x < 1; x += hh) {
        for (float y = 0; y < 1; y += hh) {
            if (indicatef(x,y)) {
                p[i] = Vec(x, y);
                v[i] = Vec(0, 0);
                ++i;
            }
        }
    }
    
    normalizeMass();
}

void SPHSystem::normalizeMass() {
    mass = 1.0;

    computeDensity();
    
    float rho0 = params.rho0;
    float rho2s = 0;
    float rhos  = 0;
    for (int i = 0; i < n; ++i) {
        rho2s += density[i]*density[i];
        rhos  += density[i];
    }

    mass *= ( rho0*rhos / rho2s );
}

void SPHSystem::computeDensity() {
    vector<float>& rho = density;
    
    float h  = params.h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * mass / M_PI / h8;
    
    memset(&(rho[0]), 0, sizeof(float)*rho.size());
    
    for (int i = 0; i < n; ++i) {
        rho[i] += 4 * mass / M_PI / h2;
        for (int j = i+1; j < n; ++j) {
            Vec dp = p[i] - p[j];
            float r2 = glm::dot(dp, dp);
            float z  = h2-r2;
            if (z > 0) {
                float rho_ij = C*z*z*z;
                rho[i] += rho_ij;
                rho[j] += rho_ij;
            }
        }
    }
}

void SPHSystem::computeAcceleration() {
    // Unpack basic parameters
    const float h    = params.h;
    const float rho0 = params.rho0;
    const float k    = params.k;
    const float mu   = params.mu;
    const float g    = params.g;
    const float h2   = h*h;
    
    // Unpack system state
    const vector<float>& rho = density;
    
    // Compute density and color
    computeDensity();
    
//    for ( auto x : rho ) {
//        cout << x << endl;
//    }
//    int dummy;
//    cin >> dummy;
    
    // Start with gravity and surface forces
    for (int i = 0; i < n; ++i) {
        a[i] = Vec(0, -g);
    }
    
    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;
    
    // Now compute interaction forces
    for (int i = 0; i < n; ++i) {
        const float rhoi = rho[i];
        for (int j = i+1; j < n; ++j) {
            Vec dp = p[i] - p[j];
            float r2 = glm::dot(dp, dp);
            if (r2 < h2) {
                const float rhoj = rho[j];
                float q = sqrt(r2)/h;
                float u = 1-q;
                float w0 = C0 * u/rhoi/rhoj;
                float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                float wv = w0 * Cv;
                Vec dv = v[i] - v[j];
                a[i].x += (wp*dp.x + wv*dv.x);
                a[i].y += (wp*dp.y + wv*dv.y);
                a[j].x -= (wp*dp.x + wv*dv.x);
                a[j].y -= (wp*dp.y + wv*dv.y);
            }
        }
    }
}

void SPHSystem::dampReflect(int side, float barrier, Vec& pos, Vec& velo, Vec& veloh) {
    // Coefficient of resitiution
    const float DAMP = 0.75;
    
    // Ignore degenerate cases
    if (velo[side] == 0)
        return;
    
    // Scale back the distance traveled based on time from collision
    float tbounce = (pos[side]-barrier)/velo[side];
    pos[0] -= velo[0]*(1-DAMP)*tbounce;
    pos[1] -= velo[1]*(1-DAMP)*tbounce;
    
    // Reflect the position and velocity
    pos[side]  = 2*barrier-pos[side];
    velo[side]  = -velo[side];
    veloh[side] = -veloh[side];
    
    // Damp the velocities
    velo[0] *= DAMP;  veloh[0] *= DAMP;
    velo[1] *= DAMP;  veloh[1] *= DAMP;
}

void SPHSystem::reflectOnBoundary() {
    // Boundaries of the computational domain
    const float XMIN = 0.0;
    const float XMAX = 1.0;
    const float YMIN = 0.0;
    const float YMAX = 1.0;
    
    for (int i = 0; i < n; ++i) {
        Vec& pos = p[i];
        Vec& velo = v[i];
        Vec& veloh = vh[i];
        if (pos.x < XMIN) dampReflect(0, XMIN, pos, velo, veloh);
        if (pos.x > XMAX) dampReflect(0, XMAX, pos, velo, veloh);
        if (pos.y < YMIN) dampReflect(1, YMIN, pos, velo, veloh);
        if (pos.y > YMAX) dampReflect(1, YMAX, pos, velo, veloh);
    }
}

void SPHSystem::leapFrogStart() {
    const float dt = params.dt;
    for (int i = 0; i < n; ++i) vh[i]  = v[i] + a[i] * dt / 2.0f;
    for (int i = 0; i < n; ++i) v[i]  += a[i]  * dt;
    for (int i = 0; i < n; ++i) p[i]  += vh[i] * dt;
    reflectOnBoundary();
}

void SPHSystem::leapFrog() {
    const float dt = params.dt;
    for (int i = 0; i < n; ++i) vh[i] += a[i]  * dt;
    for (int i = 0; i < n; ++i) v[i]   = vh[i] + a[i] * dt / 2.0f;
    for (int i = 0; i < n; ++i) p[i]  += vh[i] * dt;
    reflectOnBoundary();
}

void SPHSystem::checkState() {
    for (int i = 0; i < n; ++i) {
        float xi = p[i].x;
        float yi = p[i].y;
        assert( xi >= 0 || xi <= 1 );
        assert( yi >= 0 || yi <= 1 );
    }
}