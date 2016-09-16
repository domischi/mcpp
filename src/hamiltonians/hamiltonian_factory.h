#ifndef MCPP_HAMILTONIAN_FACTORY_H
#define MCPP_HAMILTONIAN_FACTORY_H

#include <vector>
#include <memory> //unique_ptr

#include "xy_dipolar.h"

std::vector<std::shared_ptr<XY_Hamiltonian>> HamiltonianFactory(alps::Parameters const& params){
    std::vector<std::shared_ptr<XY_Hamiltonian>> HamiltonianList;
    if(params.defined("D") && static_cast<double>(params["D"])!=0.){
        HamiltonianList.push_back(std::make_shared<XY_Dipole>(params)); 
    }
    return HamiltonianList;
}

#endif //MCPP_HAMILTONIAN_FACTORY_H
