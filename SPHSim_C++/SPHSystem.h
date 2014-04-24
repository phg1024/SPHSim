//
//  SPHSystem.h
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#ifndef __SPHSim_C____SPHSystem__
#define __SPHSim_C____SPHSystem__

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include <cstdio>

#include "glm.hpp"

class SPHSystem {
public:
    struct Parameters {
        Parameters() {
            setDefaultParams();
        }
        Parameters(const string& filename) {
            loadParamsFromFile(filename);
        }
        void setDefaultParams() {
            logfile = "result.bin";
            nframes = 400;
            nsteps = 100;
            dt = 1e-4;
            h = 5e-2;
            rho0 = 1000;
            k = 1e3;
            mu = 5.0;
            g = 9.8;
        }
        
        void loadParamsFromFile(const string& filename) {
            ifstream fin(filename);
            if( !fin ) {
                cerr << "Failed to load parameters from file " << filename << endl;
                setDefaultParams();
                return;
            }
            
            fin >> logfile
                >> nframes
                >> nsteps
                >> dt
                >> h
                >> rho0
                >> k
                >> mu
                >> g;
        }
        
        string logfile;
        int nframes;        // total number of frames
        int nsteps;         // number of steps per frame
        float dt;           // time step
        float h;            // particle size
        float rho0;         // rest density
        float k;            // bulk modulus
        float mu;           // viscosity
        float g;            // gravity magnitude
    };
    
    struct Grid {
        typedef vector<int> cell_t;
        Grid(){}
        Grid(int w, int h):w(w),h(h) { init(); }
        void init() {
            cells.resize(w * h);
        }
        void reset() {
            cells.clear();
            init();
        }
        cell_t& getcell(int x, int y) {
            int idx = y*w+x;
            if(idx >=0 && idx < w*h) return cells[idx];
            else throw "Try to access ghost cell.";
        }
        
        vector<cell_t> cells;
        int w, h;
    };
    
public:
    typedef glm::vec2 Vec;
    
    SPHSystem();
    SPHSystem(const string& filename);
    
    void run();
    
protected:
    void resize(int n);
    void initParticles();
    void normalizeMass();
    
    void sortParticles();
    
    void computeDensity();
    void computeAcceleration();
    
    void leapFrogStart();
    void leapFrog();
    
    void dampReflect(int, float, Vec&, Vec&, Vec&);
    void reflectOnBoundary();
    
    void checkState();
    
    void writeHeader(ofstream&);
    void writeFrame(ofstream&, float*);
    
private:
    Parameters params;
    
    vector<Vec> p, v, vh, a;
    
    vector<float> density;
    vector<float> pressure;
    
    float mass;
    unsigned int n;
    
    Grid pgrid;
};

#endif /* defined(__SPHSim_C____SPHSystem__) */
