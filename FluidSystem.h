#pragma once

#include <vector>

class FluidSystem {
    public: 
    FluidSystem(int N, float dt, float diff, float visc, float force, float source
    , int dvel);

    int N;
    float dt;
    float diff;
    float visc;
    float force;
    float source;
    int dvel;
    std::vector<float> uForces;
    std::vector<float> vForces;
    std::vector<float> uForcesPrev;
    std::vector<float> vForcesPrev;
    std::vector<float> densities;
    std::vector<float> densitiesPrev;

    void free_data();
    void clear_data();
    int allocate_data();
    void draw_velocity();
    void draw_density();

    void add_source(std::vector<float> &x, std::vector<float> &s);
    void linear_solve(int b, std::vector<float> &x, std::vector<float> &x0, float a, float c);
    void diffuse(int b, std::vector<float> &x, std::vector<float> &x0);
    void advect(int b, std::vector<float> &d, std::vector<float> &d0, std::vector<float> &u, std::vector<float> &v);
    void project(std::vector<float> &u, std::vector<float> &v, std::vector<float> &p, std::vector<float> &div);
    void density_step();
    void velocity_step();
    void set_boundary(int b, std::vector<float> &x);
};