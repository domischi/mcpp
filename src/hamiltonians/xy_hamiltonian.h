#ifndef MCPP_XY_HAMILTONIAN_H
#define MCPP_XY_HAMILTONIAN_H

#include <alps/config.h>
#include <alps/lattice.h>
#include <alps/model.h>

// This class serves as an interface for all the interesting Hamiltonians, one wants to define for XY systems, these could for example include dipolar interaction, exchange, random field, local anisotropy,...
template<typename DERIVED>
class XY_Hamiltonian {
public:
    typedef double spin_t;
    XY_Hamiltonian (alps::Parameters const&) {}
    double SingleSiteEnergy(std::vector<double> const& spins, int i){
        return derived().SingleSiteEnergy(spins,i);  
    }
    double Energy( std::vector<double> const& spins) {
        double E=0.;
        for(int i=0;i<spins.size();++i) {
            E+=derived().SingleSiteEnergy(spins,i);  
        } 
        return E;
    }
private:
    //Convenience function for CRTP
    DERIVED& derived() {
        return *static_cast<DERIVED*>(this);
    }
};

#endif //MCPP_XY_HAMILTONIAN_H
