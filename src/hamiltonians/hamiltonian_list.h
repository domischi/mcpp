#ifndef MCPP_HAMILTONIAN_LIST_H
#define MCPP_HAMILTONIAN_LIST_H

#include <vector>
#include <memory> //unique_ptr
#include <iostream>
#include "xy_dipolar.h"
#include "xy_exchange.h"
#include "xy_j1j2.h"
#include "xy_shape_anisotropy.h"
#include "../utilities.h"
struct Hamiltonian_List {
    std::unique_ptr<XY_Dipole<true>>     p_dipolar_disordered;
    std::unique_ptr<XY_Dipole<false>>    p_dipolar;
    std::unique_ptr<XY_J1J2>             p_J1J2;
    std::unique_ptr<XY_Exchange>         p_exchange;
    std::unique_ptr<XY_Shape_Anisotropy> p_shape_anisotropy;
    Hamiltonian_List(alps::Parameters const& params) {

        if(params.defined("D") && static_cast<double>(params["D"])!=0.){
            if(mcpp::is_disordered(params) || static_cast<bool>(params["LATTICE"]!="square lattice")){
                p_dipolar_disordered=std::unique_ptr<XY_Dipole<true>>(new XY_Dipole<true>(params)); 
            }
            else{
                p_dipolar=std::unique_ptr<XY_Dipole<false>>(new XY_Dipole<false>(params));
            }
        }
        if(params.defined("J") && static_cast<double>(params["J"])!=0.){
            p_exchange=std::unique_ptr<XY_Exchange>(new XY_Exchange(params)); 
        }
        if(params.defined("J1") && static_cast<double>(params["J1"])!=0. &&
           params.defined("J2")
                ){
            p_J1J2=std::unique_ptr<XY_J1J2>(new XY_J1J2(params)); 
        }
        if( params.defined("Shape Anisotropy Strength") && static_cast<double>(params["Shape Anisotropy Strength"]) !=0. &&
            params.defined("Shape Anisotropy p")        && static_cast<double>(params["Shape Anisotropy p"])        !=0.){
            p_shape_anisotropy=std::unique_ptr<XY_Shape_Anisotropy>(new XY_Shape_Anisotropy(params)); 
        }
        if(!p_dipolar_disordered && !p_dipolar && !p_J1J2 && !p_exchange && !p_shape_anisotropy) {
            std::cerr << "NO HAMILTONIAN LOADED";
            std::exit(56);
        }
    }
    double Energy(std::vector<double> const& spins) const{
        double E=0.; 
        SingleHamiltonianAllSites(p_dipolar_disordered, spins, E);
        SingleHamiltonianAllSites(p_dipolar,            spins, E);
        SingleHamiltonianAllSites(p_exchange,           spins, E);
        SingleHamiltonianAllSites(p_J1J2,               spins, E);
        SingleHamiltonianAllSites(p_shape_anisotropy,   spins, E);
        return E;
    }
    double SingleSiteEnergy(std::vector<double> const& spins, int i) const{
        double e=0.;
        SingleHamiltonianSingleSite(p_dipolar_disordered, spins, i, e);
        SingleHamiltonianSingleSite(p_dipolar,            spins, i, e);
        SingleHamiltonianSingleSite(p_exchange,           spins, i, e);
        SingleHamiltonianSingleSite(p_J1J2,               spins, i, e);
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

std::shared_ptr<Hamiltonian_List> get_hl(mcpp::parameter_type p) {
    return std::make_shared<Hamiltonian_List>(p);
}

#endif //MCPP_HAMILTONIAN_LIST_H
