#ifndef MCPP_COMPONENT_OBSERVABLES_H_
#define MCPP_COMPONENT_OBSERVABLES_H_

#include <alps/parameter.h> 
#include "../observable.h"
#include "../../utilities.h"
class component_observables : public observable{
public:
    typedef double spin_t;

    component_observables(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
        observable(p,gh_,hl_)
    {
    }
                
    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) {
        double mxs,mys, ms;
        std::tie(mxs,mys)=mcpp::mxs_and_mys(spins,L);
        ms=std::sqrt(mxs*mxs+mys*mys);
        obs["M striped"    ]<< std::abs((mxs*mxs-mys*mys)/ms); 
        obs["M microvortex"]<< std::abs(2*mxs*mys/ms);
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs <<alps::RealObservable("M striped");
        obs <<alps::RealObservable("M microvortex");
    }

    // Intentionally left empty as only const values are inside the class (initialization due to constructor)
    void save(alps::ODump &dump) const{ }
private:

};
#endif //MCPP_COMPONENT_OBSERVABLES_H_
