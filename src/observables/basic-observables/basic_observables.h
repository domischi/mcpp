#ifndef MCPP_BASIC_OBSERVABLES_H_
#define MCPP_BASIC_OBSERVABLES_H_

#include <alps/parameter.h>
#include "../../utilities.h"
class basic_observables {
public:
    typedef double spin_t;

    basic_observables(const alps::Parameters& p) :
    L(p["L"]),
    N(mcpp::init_N(p))
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
    const int L,N;

};
#endif //MCPP_BASIC_OBSERVABLES_H_
