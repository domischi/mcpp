#ifndef MCPP_XY_HAMILTONIAN_J1J2_H
#define MCPP_XY_HAMILTONIAN_J1J2_H

#include "xy_hamiltonian.h"

// Apply the convention that J>0 -> Ferromagnetic
class XY_J1J2 : public XY_Hamiltonian<XY_J1J2>{
public:
    XY_J1J2 (alps::Parameters const& params) : 
    XY_Hamiltonian<XY_J1J2>(params),
    J1(params.value_or_default("J1",1)),
    J2(params.value_or_default("J2",1)),
    L(params["L"]),
    gh(params){
    }
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i) const {
        double e=0.;
        neighbour_bond_iterator b, b_end;
        for(std::tie(b,b_end)=gh.neighbor_bonds(i); b!= b_end;++b){
            e-=(gh.bond_type(*b) ? J2 : J1)* std::cos(spins[i]-spins[gh.target(*b)]);
        }
        return e/2.;
    }
    
    typedef alps::graph_helper<>::neighbor_iterator neighbour_iterator;
    typedef alps::graph_helper<>::neighbor_bond_iterator neighbour_bond_iterator;

private:
    const int L;
    const double J1,J2;
    alps::graph_helper<> gh;
};

#endif //MCPP_XY_HAMILTONIAN_J1J2_H
