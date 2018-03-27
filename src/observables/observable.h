#ifndef MCPP_OBSERVABLE_H_
#define MCPP_OBSERVABLE_H_

// This severes as an interface class for the various observables
#include <memory>
#include <alps/parameter.h>
#include <alps/lattice.h>
#include "../utilities.h"
class observable {
public:
    typedef double spin_t;

    observable(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_ ) :
        graph_helper_(gh_),
        hamiltonian_list_(hl_),
        L(p["L"])
    {
    }
                
    virtual void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) = 0; 
    virtual void init_observables(alps::ObservableSet& obs) const = 0; 
    virtual void save(alps::ODump &dump) const{ }
    virtual void load(alps::IDump &dump) const{ }
protected:
    std::shared_ptr<alps::graph_helper<>> graph_helper_;
    std::shared_ptr<Hamiltonian_List> hamiltonian_list_;
    const int L;
};
#endif //MCPP_OBSERVABLE_H_
