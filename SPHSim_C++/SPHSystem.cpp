//
//  SPHSystem.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include "SPHSystem.h"
#include "marchingcubes/CIsoSurface.h"

SPHSystem::SPHSystem() {
    kern = Kernel(params.h);
    initParticles();
}

SPHSystem::SPHSystem(const string& filename) {
    params = Parameters(filename);
    kern = Kernel(params.h);
    
    initParticles();
}

void SPHSystem::writeHeader() {
    float scale = 1.0;
    
    logstream.write(reinterpret_cast<const char*>(&n), sizeof(int));
    logstream.write(reinterpret_cast<const char*>(&params.nframes), sizeof(int));
    logstream.write(reinterpret_cast<const char*>(&scale), sizeof(float));
}

void SPHSystem::writeFrame(float* c = 0) {
    for (int i = 0; i < n; ++i) {
        Vec x = p[i];
        float ci = c ? c[i] : 0;
        //cout << p[i].x << ' ' << p[i].y << ' '<< p[i].z << ' ' << ci << endl;
        
        logstream.write(reinterpret_cast<const char*>(&x.x), sizeof(float));
        logstream.write(reinterpret_cast<const char*>(&x.y), sizeof(float));
        logstream.write(reinterpret_cast<const char*>(&x.z), sizeof(float));
        logstream.write(reinterpret_cast<const char*>(&ci), sizeof(float));
    }
    
    vector<float> voxels = voxelize(64, 64, 64);
    // write the voxels to a file
    stringstream ss;
    ss << "frame" << curframe << ".bin";
    cout << "writing file " << ss.str() << " with " << voxels.size() << " voxels." << endl;
    ofstream fout(ss.str(), ios::binary);
    fout.write(reinterpret_cast<const char*>(&voxels[0]), sizeof(float)*voxels.size());
    fout.close();
    
    CIsoSurface<float> mc;
    cout << "extracing mesh with marching cubes ..." << endl;
    mc.GenerateSurface(&voxels[0], 1.0f, 63, 63, 63, 1.0, 1.0, 1.0);
    cout << "done." << endl;
    stringstream ss2;
    ss2 << "frame" << curframe << ".obj";
    cout << "writing file " << ss2.str() << endl;
    mc.writeToFile(ss2.str());
}

void SPHSystem::init() {
    curframe = 0;

    if( params.writeLog ) {
        logstream = ofstream(params.logfile, ios::binary);
        // write initial state
        writeHeader();
        writeFrame(&pressure[0]);
    }
    
    // step once
    computeAcceleration();
    leapFrogStart();
    checkState();
}

void SPHSystem::step() {
    curframe++;
    cout << "frame " << curframe << " out of " << params.nframes << endl;
    if (curframe < params.nframes) {
        for (int i = 0; i < params.nsteps; ++i) {
            computeAcceleration();
            leapFrog();
            checkState();
        }
        
        if( params.writeLog )
            writeFrame(&pressure[0]);
    }
}

void SPHSystem::run() {
    init();
    for (int frame = 1; frame < params.nframes; ++frame) {
        step();
    }
}

void SPHSystem::resize(int count) {
    n = count;
    p.resize(n);
    v.resize(n);
    vh.resize(n);
    a.resize(n);
    cid.resize(n);
    nid.resize(n);
    density.resize(n);
    invrho.resize(n);
    pressure.resize(n);
}

void SPHSystem::initParticles() {
    float h  = params.h;
    float hh = h/1.3;
    
    auto indicatef = [](float x, float y, float z){ return x>0.875 && y<0.75;};
    
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
    
    mass *= ( rho0*rhos / rho2s );
}

void SPHSystem::sortParticles() {
    pgrid.reset();
    
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
    
    
    for( int i=0;i<p.size();i++ ) {
        int gx = min(max((int)floor(p[i].x * pgrid.w), 0), pgrid.w-1);
        int gy = min(max((int)floor(p[i].y * pgrid.h), 0), pgrid.h-1);
        int gz = min(max((int)floor(p[i].z * pgrid.d), 0), pgrid.d-1);
        Grid::cell_t& cell = pgrid.getcell(gx, gy, gz);
        cell.push_back(i);
        cid[i] = iVec(gx, gy, gz);
        
        // determine relevant neighbors
        nid[i].clear();
        nid[i].reserve(27);
        
        for(int nidx=0;nidx<27;nidx++) {
            int nx = gx + neighbors[nidx][0];
            int ny = gy + neighbors[nidx][1];
            int nz = gz + neighbors[nidx][2];
            
            if( nx >= 0 && nx < pgrid.w
               && ny >= 0 && ny < pgrid.h
               && nz >= 0 && nz < pgrid.d ) {
                nid[i].push_back(nz * pgrid.w * pgrid.h + ny * pgrid.w + nx);
            }
        }
    }
}

vector<float> SPHSystem::voxelize(int resX, int resY, int resZ) const {
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
    
    vector<float> voxels;
    voxels.reserve(resX*resY*resZ);
    
    float xStep = 1.0 / resX;
    float yStep = 1.0 / resY;
    float zStep = 1.0 / resZ;
    
    for (float zCoord = 0.5 * zStep; zCoord < 1.0; zCoord += zStep) {
        for (float yCoord = 0.5 * yStep; yCoord < 1.0; yCoord += yStep) {
            for (float xCoord = 0.5 * zStep; xCoord < 1.0; xCoord += xStep) {
                // evaluate the density at current point
                int gx = min(max((int)floor(xCoord * pgrid.w), 0), pgrid.w-1);
                int gy = min(max((int)floor(yCoord * pgrid.h), 0), pgrid.h-1);
                int gz = min(max((int)floor(zCoord * pgrid.d), 0), pgrid.d-1);
                
                Vec gp = Vec(xCoord, yCoord, zCoord);
                float rho = 0;
                
                for(int nidx=0;nidx<27;nidx++) {
                    int nx = gx + neighbors[nidx][0];
                    int ny = gy + neighbors[nidx][1];
                    int nz = gz + neighbors[nidx][2];
                    
                    if( nx >= 0 && nx < pgrid.w
                       && ny >= 0 && ny < pgrid.h
                       && nz >= 0 && nz < pgrid.d ) {
                        const Grid::cell_t& ncell = pgrid.getcell(nx, ny, nz);
                        
                        for (int j : ncell) {
                            Vec dp = gp - p[j];
                            
                            float r = glm::length(dp);
                            float rho_j = mass * kern.poly6(r);
                            rho += rho_j;
                        }
                    }
                }
                voxels.push_back(rho);
            }
        }
    }
    
    return voxels;
}

void SPHSystem::computeDensity() {
    sortParticles();
    
    vector<float>& rho = density;
    
    memset(&(rho[0]), 0, sizeof(float)*rho.size());
    
    for (int i = 0; i < n; ++i) {
        rho[i] += mass * kern.poly6(0);
        
        for(int nidx=0;nidx<nid[i].size();nidx++) {
            const Grid::cell_t& ncell = pgrid.getcell(nid[i][nidx]);
            
            for (int j : ncell) {
                if (i >= j) {
                    continue;
                }
                Vec dp = p[i] - p[j];
                
                float r = glm::length(dp);
                float rho_ij = mass * kern.poly6(r);
                rho[i] += rho_ij;
                rho[j] += rho_ij;
            }
        }
    }
    
    // update pressure
//    const float B = params.k;
//    const float gamma = 7.0;
//    const float invRHO0 = 1.0 / params.rho0;
    for (int i=0; i<n; ++i) {
        pressure[i] = params.k * (rho[i] - params.rho0);
//        pressure[i] = B * (powf(rho[i] * invRHO0, gamma) - 1);
        invrho[i] = 1.0 / rho[i];
    }
}

void SPHSystem::computeAcceleration() {
    // Unpack basic parameters
    const float mu   = params.mu;
    const float g    = params.g;
    
    // Compute density and color
    computeDensity();
    
    // Start with gravity and surface forces
    for (int i = 0; i < n; ++i) {
        a[i] = Vec(0, -g, 0);
    }
    
    // Now compute interaction forces
    for (int i = 0; i < n; ++i) {
        float w0i = mass * invrho[i];
        
        for(int nidx=0;nidx<nid[i].size();nidx++) {
            const Grid::cell_t& ncell = pgrid.getcell(nid[i][nidx]);
            
            for (int j : ncell) {
                if( i >= j ) continue;
                Vec dp = p[i] - p[j];
                float r = glm::length(dp);
                
                //Vec fp = mass / rhoj / rhoi * 0.5f * (pressure[i] + pressure[j]) * spiky_grad(h, r) * dp;
                float w0ij = w0i * invrho[j];
                Vec fp = w0ij * 0.5f * (pressure[i] + pressure[j]) * kern.spiky_grad(r) * dp;
                
                Vec dv = v[j] - v[i];
                
                Vec fv = mu * w0ij * dv * kern.viscosity_lap(r);
                
                Vec ftotal = fp + fv;
                a[i] += ftotal;
                a[j] -= ftotal;
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