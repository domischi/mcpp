#ifndef MCPP_HAMILTONIAN_FACTORY_H
#define MCPP_HAMILTONIAN_FACTORY_H

#include <vector>
#include <memory> //unique_ptr
#include <iostream>
#include "xy_dipolar.h"
#include "xy_exchange.h"
#include "xy_shape_anisotropy.h"

std::vector<std::shared_ptr<XY_Hamiltonian>> HamiltonianFactory(alps::Parameters const& params){
    std::vector<std::shared_ptr<XY_Hamiltonian>> HamiltonianList;
     if(params.defined("D") && static_cast<double>(params["D"])!=0.){
        if((params.defined("Dilution Rate") && static_cast<double>(params["Dilution Rate"])!=0.) || (params.defined("Position Disorder") && static_cast<double>(params["Position Disorder"])!=0.)){
            HamiltonianList.push_back(std::make_shared<XY_Dipole<true>>(params));
        }
        else
            HamiltonianList.push_back(std::make_shared<XY_Dipole<false>>(params));
    }
    if(params.defined("J") && static_cast<double>(params["J"])!=0.){
        HamiltonianList.push_back(std::make_shared<XY_Exchange>(params)); 
    }
    if(params.defined("p_max") && static_cast<double>(params["A"])!=0.){
        HamiltonianList.push_back(std::make_shared<XY_Shape_Anisotropy>(params)); 
    }
    if(HamiltonianList.size()==0) {
        std::cerr << "NO HAMILTONIAN LOADED";
        std::exit(56);
    }
    return HamiltonianList;
}

#endif //MCPP_HAMILTONIAN_FACTORY_H
