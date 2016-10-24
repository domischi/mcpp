#ifndef MCPP_HAMILTONIAN_LIST_H
#define MCPP_HAMILTONIAN_LIST_H

#include <vector>
#include <memory> //unique_ptr
#include <iostream>
#include "xy_dipolar.h"
#include "xy_exchange.h"
#include "xy_shape_anisotropy.h"

struct Hamiltonian_List {
    
    std::unique_ptr<XY_Dipole<true>>     p_dipolar_disordered;
    std::unique_ptr<XY_Dipole<false>>    p_dipolar;
    std::unique_ptr<XY_Exchange>         p_exchange;
    std::unique_ptr<XY_Shape_Anisotropy> p_shape_anisotropy;

    Hamiltonian_List(alps::Parameters const& params) {
         if(params.defined("D") && static_cast<double>(params["D"])!=0.){
            if((params.defined("Dilution Rate") && static_cast<double>(params["Dilution Rate"])!=0.) || (params.defined("Position Disorder") && static_cast<double>(params["Position Disorder"])!=0.)){
                p_dipolar_disordered=std::unique_ptr<XY_Dipole<true>>(new XY_Dipole<true>(params)); 
            }
            else
                p_dipolar=std::unique_ptr<XY_Dipole<false>>(new XY_Dipole<false>(params)); 
        }
        if(params.defined("J") && static_cast<double>(params["J"])!=0.){
            p_exchange=std::unique_ptr<XY_Exchange>(new XY_Exchange(params)); 
        }
        if( params.defined("A")     && static_cast<double>(params["A"])    !=0.
            params.defined("p_max") && static_cast<double>(params["p_max"])!=0.){
            p_shape_anisotropy=std::unique_ptr<XY_Shape_Anisotropy>(new XY_Shape_Anisotropy(params)); 
        }
        if(!p_dipolar_disordered && !p_dipolar && !p_exchange && !p_shape_anisotropy) {
            std::cerr << "NO HAMILTONIAN LOADED";
            std::exit(56);
        }
    }
    double Energy(std::vector<double> const& spins) const{
        double E=0.; 
        SingleHamiltonianAllSites(p_dipolar_disordered, spins, E);
        SingleHamiltonianAllSites(p_dipolar,            spins, E);
        SingleHamiltonianAllSites(p_exchange,           spins, E);
        SingleHamiltonianAllSites(p_shape_anisotropy,   spins, E);
        return E;
    }
    double SingleSiteEnergy(std::vector<double> const& spins, int i) const{
        double e=0.;
        SingleHamiltonianSingleSite(p_dipolar_disordered, spins, i, e);
        SingleHamiltonianSingleSite(p_dipolar,            spins, i, e);
        SingleHamiltonianSingleSite(p_exchange,           spins, i, e);
        SingleHamiltonianSingleSite(p_shape_anisotropy,   spins, i, e);
        return e;
    }

    template<typename p_t>
    inline void SingleHamiltonianSingleSite(p_t const& p, std::vector<double> const& spins, int i, double& e) const {
        if(p) e+=p->SingleSiteEnergy(spins,i);
    }
    template<typename p_t>
    inline void SingleHamiltonianAllSites(p_t const& p, std::vector<double> const& spins, double& e) const {
        if(p) e+=p->Energy(spins);
    }
};

#endif //MCPP_HAMILTONIAN_LIST_H
