//
//  SPHSystem.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include "SPHSystem.h"
#include "kernels.h"

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
        //cout << p[i].x << ' ' << p[i].y << ' '<< p[i].z << ' ' << ci << endl;
        
        f.write(reinterpret_cast<const char*>(&x.x), sizeof(float));
        f.write(reinterpret_cast<const char*>(&x.y), sizeof(float));
        f.write(reinterpret_cast<const char*>(&x.z), sizeof(float));
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
    
    auto indicatef = [](float x, float y, float z){ return x>0.5 && y<0.5 && z>0.5;};
    
    // Count mesh points that fall in indicated region.
    int count = 0;
    for (float x = 0; x < 1; x += hh)
        for (float y = 0; y < 1; y += hh)
            for (float z = 0; z < 1; z += hh)
                count += indicatef(x, y, z);
    
    // Populate the particle data structure
    resize(count);
    
    int i = 0;
    for (float x = 0; x < 1; x += hh) {
        for (float y = 0; y < 1; y += hh) {
            for (float z = 0; z < 1; z += hh) {
                if (indicatef(x, y, z)) {
                    p[i] = Vec(x, y, z);
                    //cout << p[i].x << ' ' << p[i].y << ' '<< p[i].z << endl;
                    v[i] = Vec(0, 0, 0);
                    ++i;
                }
            }
        }
    }

    int gw = ceil(0.99/h);
    int gh = ceil(0.99/h);
    int gd = ceil(0.99/h);
    
    pgrid = Grid(gw, gh, gd);
    
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
    
    // mass = rho0 * sum( rho[i] ) / sum( rho[i] * rho[i] )
    mass *= ( rho0*rhos / rho2s );
}

void SPHSystem::sortParticles() {
    pgrid.reset();
    
    for( int i=0;i<p.size();i++ ) {
        int gx = min(max((int)floor(p[i].x * pgrid.w), 0), pgrid.w-1);
        int gy = min(max((int)floor(p[i].y * pgrid.h), 0), pgrid.h-1);
        int gz = min(max((int)floor(p[i].z * pgrid.d), 0), pgrid.d-1);
        Grid::cell_t& cell = pgrid.getcell(gx, gy, gz);
        cell.push_back(i);
    }
}

void SPHSystem::computeDensity() {
    sortParticles();
    
    vector<float>& rho = density;
    
    float h  = params.h;
    float h3 = h*h*h;
    
    memset(&(rho[0]), 0, sizeof(float)*rho.size());
    
    const int neighbors[][3] = {
        {-1, -1, -1}, {-1, 0, -1}, {-1, 1, -1},
        { 0, -1, -1}, { 0, 0, -1}, { 0, 1, -1},
        { 1, -1, -1}, { 1, 0, -1}, { 1, 1, -1},

        {-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0},
        { 0, -1, 0}, { 0, 0, 0}, { 0, 1, 0},
        { 1, -1, 0}, { 1, 0, 0}, { 1, 1, 0},

        {-1, -1, 1}, {-1, 0, 1}, {-1, 1, 1},
        { 0, -1, 1}, { 0, 0, 1}, { 0, 1, 1},
        { 1, -1, 1}, { 1, 0, 1}, { 1, 1, 1},
    };
    
    for (int i = 0; i < n; ++i) {
        rho[i] += mass * poly6(h, 0);
        
        int gx = max(floor(p[i].x * pgrid.w), 0.0f);
        int gy = max(floor(p[i].y * pgrid.h), 0.0f);
        int gz = max(floor(p[i].z * pgrid.d), 0.0f);
        
        //const Grid::cell_t& cell = pgrid.getcell(gx, gy);
        for(int nidx=0;nidx<27;nidx++) {
            int nx = gx + neighbors[nidx][0];
            int ny = gy + neighbors[nidx][1];
            int nz = gz + neighbors[nidx][2];
            
            if( nx >= 0 && nx < pgrid.w
             && ny >= 0 && ny < pgrid.h
             && nz >= 0 && nz < pgrid.d ) {
                const Grid::cell_t& ncell = pgrid.getcell(nx, ny, nz);
                
                for (int j : ncell) {
                    if (i >= j) {
                        continue;
                    }
                    Vec dp = p[i] - p[j];
                    
                    float r = glm::length(dp);
                    float rho_ij = mass * poly6(h, r);
                    rho[i] += rho_ij;
                    rho[j] += rho_ij;
                }
            }
        }
    }
    
    // update pressure
    const float B = params.k;
    const float gamma = 7.0;
    const float invRHO0 = 1.0 / params.rho0;
    for (int i=0; i<n; ++i) {
        //pressure[i] = params.k * (rho[i] - params.rho0);
        pressure[i] = B * (powf(rho[i] * invRHO0, gamma) - 1);
    }
}

void SPHSystem::computeAcceleration() {
    // Unpack basic parameters
    const float h    = params.h;
    //const float rho0 = params.rho0;
    //const float k    = params.k;
    const float mu   = params.mu;
    const float g    = params.g;
    
    // Unpack system state
    const vector<float>& rho = density;
    
    // Compute density and color
    computeDensity();
    
    // Start with gravity and surface forces
    for (int i = 0; i < n; ++i) {
        a[i] = Vec(0, -g, 0);
    }
    
    const int neighbors[][3] = {
        {-1, -1, -1}, {-1, 0, -1}, {-1, 1, -1},
        { 0, -1, -1}, { 0, 0, -1}, { 0, 1, -1},
        { 1, -1, -1}, { 1, 0, -1}, { 1, 1, -1},
        
        {-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0},
        { 0, -1, 0}, { 0, 0, 0}, { 0, 1, 0},
        { 1, -1, 0}, { 1, 0, 0}, { 1, 1, 0},
        
        {-1, -1, 1}, {-1, 0, 1}, {-1, 1, 1},
        { 0, -1, 1}, { 0, 0, 1}, { 0, 1, 1},
        { 1, -1, 1}, { 1, 0, 1}, { 1, 1, 1},
    };
    
    // Now compute interaction forces
    for (int i = 0; i < n; ++i) {
        const float rhoi = rho[i];
        int gx = max(floor(p[i].x * pgrid.w), 0.0f);
        int gy = max(floor(p[i].y * pgrid.h), 0.0f);
        int gz = max(floor(p[i].z * pgrid.d), 0.0f);
        
        //const Grid::cell_t& cell = pgrid.getcell(gx, gy);
        for(int nidx=0;nidx<27;nidx++) {
            int nx = gx + neighbors[nidx][0];
            int ny = gy + neighbors[nidx][1];
            int nz = gz + neighbors[nidx][2];
            
            if( nx >= 0 && nx < pgrid.w
             && ny >= 0 && ny < pgrid.h
             && nz >= 0 && nz < pgrid.d ) {
                const Grid::cell_t& ncell = pgrid.getcell(nx, ny, nz);
                
                for (int j : ncell) {
                    if( i >= j ) continue;
                    Vec dp = p[i] - p[j];
                    float r = glm::length(dp);
                    const float rhoj = rho[j];
                    
                    Vec fp = mass / rhoj / rhoi * 0.5f * (pressure[i] + pressure[j]) * spiky_grad(h, r) * dp;
                    
                    Vec dv = v[j] - v[i];

                    Vec fv = mu * mass / rhoj / rhoi * dv * viscosity_lap(h, r);
                    
                    Vec ftotal = fp + fv;
                    a[i] += ftotal;
                    a[j] -= ftotal;
                }
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
    
    pos -= velo * (1-DAMP) * tbounce;
    
    // Reflect the position and velocity
    pos[side]  = 2*barrier-pos[side];
    velo[side]  = -velo[side];
    veloh[side] = -veloh[side];
    
    // Damp the velocities
    velo *= DAMP;   veloh *= DAMP;
}

void SPHSystem::reflectOnBoundary() {
    // Boundaries of the computational domain
    const float XMIN = 0.0;
    const float XMAX = 1.0;
    const float YMIN = 0.0;
    const float YMAX = 1.0;
    const float ZMIN = 0.0;
    const float ZMAX = 1.0;

    for (int i = 0; i < n; ++i) {
        Vec& pos = p[i];
        Vec& velo = v[i];
        Vec& veloh = vh[i];
        if (pos.x < XMIN) dampReflect(0, XMIN, pos, velo, veloh);
        if (pos.x > XMAX) dampReflect(0, XMAX, pos, velo, veloh);
        if (pos.y < YMIN) dampReflect(1, YMIN, pos, velo, veloh);
        if (pos.y > YMAX) dampReflect(1, YMAX, pos, velo, veloh);
        if (pos.z < ZMIN) dampReflect(2, ZMIN, pos, velo, veloh);
        if (pos.z > ZMAX) dampReflect(2, ZMAX, pos, velo, veloh);
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