#ifndef MCPP_XY_HAMILTONIAN_EXCHANGE_H
#define MCPP_XY_HAMILTONIAN_EXCHANGE_H

#include "xy_hamiltonian.h"


// Apply the convention that J>0 -> Ferromagnetic
class XY_Exchange : public XY_Hamiltonian<XY_Exchange>{
public:
    //typedef typename XY_Hamiltonian<G>::spin_t spin_t;
    XY_Exchange (alps::Parameters const& params) : 
    XY_Hamiltonian<XY_Exchange>(params),
    J(params.value_or_default("J",1)),
    L(params["L"]),
    gh(params){
    }
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i){
        double e=0.;
        neighbour_iterator j, j_end;
        for(std::tie(j,j_end)=gh.neighbors(i); j!= j_end;++j){
            e -= J * std::cos(spins[i]-spins[*j]);
        }
        return e;
    }
    
    typedef alps::graph_helper<>::neighbor_iterator neighbour_iterator;

private:
    int L;
    double J;
    alps::graph_helper<> gh;
};

#endif //MCPP_XY_HAMILTONIAN_EXCHANGE_H

