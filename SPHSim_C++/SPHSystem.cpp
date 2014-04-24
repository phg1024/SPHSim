//
//  SPHSystem.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include "SPHSystem.h"
#define USE_BUCKETS 1

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
    writeFrame(f, &pressure[0]);
    
    computeAcceleration();
    leapFrogStart();
    checkState();
    
    for (int frame = 1; frame < params.nframes; ++frame) {
        for (int i = 0; i < params.nsteps; ++i) {
            computeAcceleration();
            leapFrog();
            checkState();
        }
        writeFrame(f, &pressure[0]);
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
    int gw = ceil(0.99/h);
    int gh = ceil(0.99/h);
    
    pgrid = Grid(gw, gh);
    
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

void SPHSystem::sortParticles() {
    pgrid.reset();
    
    for( int i=0;i<p.size();i++ ) {
        int gx = min(max((int)floor(p[i].x * pgrid.w), 0), pgrid.w-1);
        int gy = min(max((int)floor(p[i].y * pgrid.h), 0), pgrid.h-1);
        Grid::cell_t& cell = pgrid.getcell(gx, gy);
        cell.push_back(i);
    }
}

void SPHSystem::computeDensity() {
#if USE_BUCKETS
    sortParticles();
#endif

    vector<float>& rho = density;
    
    float h  = params.h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * mass / M_PI / h8;
    
    memset(&(rho[0]), 0, sizeof(float)*rho.size());
    
    const int neighbors[][2] = {
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1}, { 0, 0}, { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1}
    };

#if USE_BUCKETS
    for (int i = 0; i < n; ++i) {
        rho[i] += 4 * mass / M_PI / h2;
        int gx = max(floor(p[i].x * pgrid.w), 0.0f);
        int gy = max(floor(p[i].y * pgrid.h), 0.0f);
        
        const Grid::cell_t& cell = pgrid.getcell(gx, gy);
        for(int nidx=0;nidx<9;nidx++) {
            int nx = gx + neighbors[nidx][0];
            int ny = gy + neighbors[nidx][1];
            if( nx >= 0 && nx < pgrid.w && ny >= 0 && ny < pgrid.h ) {
                const Grid::cell_t& ncell = pgrid.getcell(nx, ny);
                
                for (int j : ncell) {
                    if (i >= j) {
                        continue;
                    }
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
    }
#else
    for (int i = 0; i < n; ++i) {
        rho[i] += 4 * mass / M_PI / h2;
        for (int j = i+1; j<n; ++j) {
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
#endif
    
    // update pressure
    for (int i=0; i<n; ++i) {
        pressure[i] = params.k * (rho[i] - params.rho0);
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
    
    // Start with gravity and surface forces
    for (int i = 0; i < n; ++i) {
        a[i] = Vec(0, -g);
    }
    
    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;
    
#if USE_BUCKETS
    const int neighbors[][2] = {
        {-1, -1}, {-1, 0}, {-1, 1},
        { 0, -1}, { 0, 0}, { 0, 1},
        { 1, -1}, { 1, 0}, { 1, 1}
    };
    
    // Now compute interaction forces
    for (int i = 0; i < n; ++i) {
        const float rhoi = rho[i];
        int gx = max(floor(p[i].x * pgrid.w), 0.0f);
        int gy = max(floor(p[i].y * pgrid.h), 0.0f);
        
        const Grid::cell_t& cell = pgrid.getcell(gx, gy);
        for(int nidx=0;nidx<9;nidx++) {
            int nx = gx + neighbors[nidx][0];
            int ny = gy + neighbors[nidx][1];
            if( nx >= 0 && nx < pgrid.w && ny >= 0 && ny < pgrid.h ) {
                const Grid::cell_t& ncell = pgrid.getcell(nx, ny);

                for (int j : ncell) {
                    if( i >= j ) continue;
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
    }
#else
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
#endif
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
    for (int i = 0; i < n; ++i) {
        p[i]  += vh[i] * dt;
    }
    reflectOnBoundary();
}

void SPHSystem::leapFrog() {
    const float dt = params.dt;
    
    for (int i = 0; i < n; ++i) vh[i] += a[i]  * dt;
    for (int i = 0; i < n; ++i) v[i]   = vh[i] + a[i] * dt / 2.0f;
    for (int i = 0; i < n; ++i) {
        p[i]  += vh[i] * dt;
    }
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