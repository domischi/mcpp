#ifndef MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H
#define MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H

#include <vector>
#include <random>
#include <cmath>
#include <iostream>

#include "xy_hamiltonian.h"
#include "../utilities.h"

// Describes a local shape anisotropy, such that a specific directions are more favorable in terms of an angle
// The parameter A describes how much this is actually the case
// The parameter p describes how many minima the actual shape anisotropy has.
class XY_Shape_Anisotropy : public XY_Hamiltonian<XY_Shape_Anisotropy>{
public:
    XY_Shape_Anisotropy (alps::Parameters const& params) : 
    XY_Hamiltonian<XY_Shape_Anisotropy>(params),
    A(params.value_or_default("Shape Anisotropy Strength",1.)),
    fixed_angle(params.defined("Shape Anisotropy Angle")),
    angle(params.value_or_default("Shape Anisotropy Angle",0.)),
    distribution(Init_Distribution_Type(params)),
    p_max(params["Shape Anisotropy p"])
    {
        Init_Lookup_Tables(params["DISORDER_SEED"], mcpp::init_N(params)); 
    }
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i) const {
        double e=0.;
        e -= A * std::cos(p[i]*(spins[i]-MinimumAtAngle[i]));
        return e;
    }
    
private:
    enum distribution_t {fixed_p, uniform_int};
    const distribution_t distribution;
    const int p_max;
    const double A;
    const bool fixed_angle;
    const double angle;
    std::vector<double> MinimumAtAngle;
    std::vector<int> p;

    distribution_t Init_Distribution_Type(alps::Parameters const& params) const {
        if(params.value_or_default("Shape Anisotropy Distribution Type", "Fixed") == "Fixed")
            return distribution_t::fixed_p;
        else if (params.value_or_default("Shape Anisotropy Distribution Type","Fixed") == "Uniform")
            return distribution_t::uniform_int;
        else{
            std::cerr << "Did not recognize Shape Anisotropy Distribution Type, aborting..." <<std::endl;
            std::exit(4);
        }
    }

    void Init_Lookup_Tables(int DISORDER_SEED, int N){
        std::mt19937 rng(DISORDER_SEED);
        std::uniform_real_distribution<double> real_dist(0,2*M_PI); 
        for(int i = 0;i<N;++i){
            if(fixed_angle)
                MinimumAtAngle.push_back(angle);
            else
                MinimumAtAngle.push_back(real_dist(rng));
        }
        switch(distribution){
            case distribution_t::fixed_p:
                p=std::vector<int>(N,p_max);
                break;
            case distribution_t::uniform_int:
                std::uniform_int_distribution<int> int_dist(1,p_max);
                for(int i = 0; i<N;++i){
                    p.push_back(int_dist(rng));
                }
                break;
        }
        MinimumAtAngle.shrink_to_fit();
        p.shrink_to_fit();
    }
};

#endif //MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H

