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
            writeLog = true;
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
            
            string dooutput;
            fin >> dooutput
                >> logfile
                >> nframes
                >> nsteps
                >> dt
                >> h
                >> rho0
                >> k
                >> mu
                >> g;
            
            if( dooutput == "true" ) writeLog = true;
            else writeLog = false;
        }
        
        bool writeLog;
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
        Grid(int w, int h, int d):w(w),h(h),d(d) { init(); }
        void init() {
            cells.resize(w * h * d);
        }
        void reset() {
            cells.clear();
            init();
        }
        cell_t& getcell(int x, int y, int z) {
            int idx = w*h*z + y*w + x;
            if(idx >=0 && idx < w*h*d) return cells[idx];
            else throw "Try to access ghost cell.";
        }
        
        vector<cell_t> cells;
        int w, h, d;
    };
    
public:
    typedef glm::vec3 Vec;
    
    SPHSystem();
    SPHSystem(const string& filename);
    
    void init();
    void step();
    void run();
    
    const vector<Vec>& particles() const {
        return p;
    }
    
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
    
    void writeHeader();
    void writeFrame(float*);
    
private:
    Parameters params;
    
    vector<Vec> p, v, vh, a;
    
    vector<float> density;
    vector<float> pressure;
    
    float mass;
    unsigned int n;
    
    Grid pgrid;
    
    int curframe;
    ofstream logstream;
};

#endif /* defined(__SPHSim_C____SPHSystem__) */
