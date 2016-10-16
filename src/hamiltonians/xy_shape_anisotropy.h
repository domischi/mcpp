#ifndef MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H
#define MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H

#include <vector>
#include <random>
#include <cmath>

#include "xy_hamiltonian.h"

// Describes a local shape anisotropy, such that a specific directions are more favorable in terms of an angle
// The parameter A describes how much this is actually the case
// The parameter p describes how many minima the actual shape anisotropy has.
class XY_Shape_Anisotropy : public XY_Hamiltonian{
public:
    XY_Shape_Anisotropy (alps::Parameters const& params) : 
    XY_Hamiltonian(params),
    A(params.value_or_default("A",1.)),
    p_max(params["p_max"])
    {
        Init_Lookup_Tables(params["DISORDER_SEED"], alps::graph_helper<>(params).num_sites()); 
    }
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i) const {
        double e=0.;
        e -= A * std::cos(p[i]*(spins[i]-MinimumAtAngle[i]));
        return e;
    }
    
private:
    int p_max;
    double A;

    std::vector<double> MinimumAtAngle;
    std::vector<int> p;

    void Init_Lookup_Tables(int DISORDER_SEED, int N){
        std::mt19937 rng(DISORDER_SEED);
        std::uniform_int_distribution<int> int_dist(1,p_max); //TODO this might not be the most meaningful distribution for p, as this randomly draws an isotropy, however if most of the dots are circular, then this would enforce a useless constraint
        std::uniform_real_distribution<double> real_dist(0,2*M_PI); 
        for(int i = 0;i<N;++i){
            MinimumAtAngle.push_back(real_dist(rng));
            p.push_back(int_dist(rng));
        }
        MinimumAtAngle.shrink_to_fit();
        p.shrink_to_fit();
    }
};

#endif //MCPP_XY_HAMILTONIAN_SHAPE_ANISOTROPY_H

