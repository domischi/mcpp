#ifndef MCPP_BASIC_OBSERVABLES_H_
#define MCPP_BASIC_OBSERVABLES_H_

#include <alps/parameter.h> 
#include "../observable.h"
#include "../../utilities.h"
class basic_observables : public observable{
public:
    typedef double spin_t;

    basic_observables(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
        observable(p,gh_,hl_)
    {
    }
                
    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) {
        double M, Ms;
        std::tie(M,Ms)=mcpp::magnetization_and_staggered_magnetization(spins,L);
        obs["M"            ] << M; 
        obs["M^2"          ] << M*M;
        obs["M^4"          ] << M*M*M*M;
        obs["M staggered"  ] << Ms;
        obs["M staggered^2"] << Ms*Ms;
        obs["M staggered^4"] << Ms*Ms*Ms*Ms;
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealObservable("M");
        obs << alps::RealObservable("M^2");
        obs << alps::RealObservable("M^4");
        obs << alps::RealObservable("M staggered");
        obs << alps::RealObservable("M staggered^2");
        obs << alps::RealObservable("M staggered^4");
    }

    // Intentionally left empty as only const values are inside the class (initialization due to constructor)
    void save(alps::ODump &dump) const{ }
private:

};
#endif //MCPP_BASIC_OBSERVABLES_H_
