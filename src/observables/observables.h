#ifndef MCPP_OBSERVABLES_H_
#define MCPP_OBSERVABLES_H_

#include <vector>
#include <memory>

#include "observable.h"
#include "mcrg/mcrg.h"
#include "basic-observables/basic_observables.h"
#include "spin-autocorrelation/spin_autocorrelation.h"
#include "structure-factor/structure_factor.h"
#include "llg/llg.h"

std::vector<std::shared_ptr<observable>> construct_observables(const alps::Parameters& p, std::shared_ptr<Hamiltonian_List> hl_ ) {
    std::vector<std::shared_ptr<observable>> observables;
    std::shared_ptr<alps::graph_helper<>> gh_=std::make_shared<alps::graph_helper<>>(p);
    if(static_cast<bool>(p.value_or_default("basic_observables",true))) {
        observables.push_back(std::make_shared<basic_observables>(p,gh_,hl_));
    }
    if(static_cast<bool>(static_cast<int>(p.value_or_default("mcrg_iteration_depth",-1))>0)) {
        int mcrg_it_depth=p["mcrg_iteration_depth"];
        std::cout << "Initialize MCRG measurement with iteration depth "<<mcrg_it_depth<<"..."<<std::flush;
        observables.push_back(std::make_shared<mcrg>(p,0,mcrg_it_depth,gh_,hl_));
        std::cout << "done"<<std::endl;;
    }
    if(static_cast<bool>(p.value_or_default("structure_factor",false))){
        std::cout << "Initialize Structure Factor measurement..."<<std::flush;
        observables.push_back(std::make_shared<structure_factor>(p,gh_,hl_));
        std::cout << "done"<<std::endl;;
    }
    if(static_cast<bool>(static_cast<int>(p.value_or_default("Spin autocorrelation analysis length",-1))>0)){
        std::cout << "Initialize Spin autocorrelation measurement..."<<std::flush;
        observables.push_back(std::make_shared<spin_autocorrelation>(p,gh_,hl_));
        std::cout << "done"<<std::endl;;
    }
    if(static_cast<bool>(p.value_or_default("llg",false))) {
        std::cout << "Initialize LLG measurement..."<<std::flush;
        observables.push_back(std::make_shared<llg>(p,gh_,hl_));
        std::cout << "done"<<std::endl;;
    }
    return observables;
}


#endif //MCPP_OBSERVABLES_H_
