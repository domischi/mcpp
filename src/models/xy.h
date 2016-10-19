#ifndef MCPP_MODELS_XY_H
#define MCPP_MODELS_XY_H

#include <alps/parapack/worker.h>
#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris.h>
#include <alps/osiris/dump.h>
#include <alps/expression.h>
#include <alps/lattice.h>

#include <thread> //sleep_for
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm> //sort
#include <utility> //pair

#include "../hamiltonians/hamiltonian_list.h"

#include "../special-observables/special_observables.h"

class xy_worker : public alps::parapack::lattice_mc_worker<>{
public :
    enum init_t {GS, Random}; 

    xy_worker(const alps::Parameters& params) :
        alps::parapack::lattice_mc_worker<>(params), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",100)),
        Measure_Sweeps(params.value_or_default("SWEEPS",5000)),
        L(params["L"]),
        N(num_sites()),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        HamiltonianList(params),
        D(params.value_or_default("D",1.)),
        cutoff_distance(static_cast<double>(params.value_or_default("cutoff_distance",3.))*static_cast<double>(params.value_or_default("a",1.))),
        Step_Number(0),
        mx(0.),
        my(0.),
        mx_stag(num_sites()),//stagered in stripes
        my_stag(0.),
        accepted(0),
        mcrg_it_depth(params.value_or_default("mcrg_iteration_depth",-1)),
        measure_mcrg(static_cast<bool>(static_cast<int>(params.value_or_default("mcrg_iteration_depth",-1))>0)),
        measure_structure_factor(static_cast<bool>(params.value_or_default("structure_factor",false))),
        Each_Measurement(params.value_or_default("Each_Measurement",15)),
        targeted_acc_ratio(params.value_or_default("Targeted Acceptance Ratio",0.5)),
        angle_dev(0.1*M_PI)
        {
            init_spins(params);
            En=Energy();
            if(measure_mcrg){
                std::cout << "\tInitialize MCRG with iteration depth "<<mcrg_it_depth<<"..."<<std::flush;
                mcrg_=std::make_shared<mcrg>(params,0,mcrg_it_depth);
                std::cout << "\tdone"<<std::endl;
            }
            if(measure_structure_factor){
                std::cout << "\tInitialize Structure Factor Measurement..."<<std::flush;
                structure_factor_=std::unique_ptr<structure_factor>(new structure_factor(params));
                std::cout << "\tdone"<<std::endl;
            }
        }
    
    void init_observables(alps::Parameters const&, alps::ObservableSet& obs){
            obs << alps::RealObservable("Energy");
            obs << alps::RealObservable("Energy^2");
            obs << alps::RealObservable("M");
            obs << alps::RealObservable("M^2");
            obs << alps::RealObservable("M^4");
            obs << alps::RealObservable("Mx");
            obs << alps::RealObservable("Mx^2");
            obs << alps::RealObservable("M staggered");
            obs << alps::RealObservable("M staggered^2");
            obs << alps::RealObservable("M staggered^4");
            obs << alps::RealObservable("Mx staggered");
            obs << alps::RealObservable("Mx staggered^2");
            obs << alps::RealObservable("Acceptance Ratio"); //Probably very useful for debugging
            if(measure_mcrg){
                mcrg_->init_observables(obs);
            }
            if(measure_structure_factor){
                structure_factor_->init_observables(obs);
            }
    }

    void save(alps::ODump &dump) const{
        dump 
        //<< L 
        //<< N
        //<< T
        << Step_Number 
        << spins 
        //<< Thermalization_Sweeps
        //<< Measure_Sweeps
        //<< D
        //<< cutoff_distance
        //<< Each_Measurement
        //<< targeted_acc_ratio
        << angle_dev
        //<< dist3_
        //<< phi
        //<< neighbour_list
        << En
        << mx
        << my
        << mx_stag
        << my_stag
        << accepted
        //<< measure_mcrg
        //<< mcrg_it_depth
        ;
    }
    void load(alps::IDump &dump){
        dump 
        //>> L 
        //>> N
        //>> T
        >> Step_Number 
        >> spins 
        //>> Thermalization_Sweeps
        //>> Measure_Sweeps
        //>> D
        //>> cutoff_distance
        //>> Each_Measurement
        //>> targeted_acc_ratio
        >> angle_dev
        //>> dist3_
        //>> phi
        //>> neighbour_list
        >> En
        >> mx
        >> my
        >> mx_stag
        >> my_stag
        >> accepted
        //>> measure_mcrg
        //>> mcrg_it_depth
        ;
    }
    void run(alps::ObservableSet& obs){
        using namespace alps::alea;
        ++Step_Number;
        int n_steps = N;
        accepted=0;
        for(int i = 0;i<n_steps;++i){// Do a typewriter sweep (possibly avoiding cache reload every step)
            update(i);
        }

        if(Step_Number&(1<<10)){ //Every 1024 steps do a random site lattice sweep to avoid ergodicity problems at 0K
            for(int i = 0;i<n_steps;++i){
                update();
            }
        }
        obs["Acceptance Ratio"]<<update_angle_deviation(accepted,n_steps);// This is also of interest before thermalization
        if(is_thermalized()&&(!(Step_Number%Each_Measurement))){
            measure(obs);
            accepted=0;
        }
    }
    
    bool is_thermalized() const {
        return Step_Number >= Thermalization_Sweeps;   
    }
    double progress() const {
        return Step_Number/(static_cast<double>(Measure_Sweeps*Each_Measurement+Thermalization_Sweeps));
    }
private:
    //System parameters
    const int Thermalization_Sweeps;
    const int Measure_Sweeps;
    int Step_Number;
    const int L;
    const int N; //L^2
    init_t init_type;
    std::vector<double> spins; //saves the spins
    const double T;
    const double D;
    const double cutoff_distance;
    const int Each_Measurement;
    
    const double targeted_acc_ratio;
    double angle_dev;

    // struct of Hamiltonian pointers -> allows CRTP
    Hamiltonian_List HamiltonianList;

    //Easy observables
    double En;
    double mx;
    double my;
    double mx_stag;
    double my_stag;
    int accepted;
    //Complex observables
    const bool measure_mcrg;
    std::shared_ptr<mcrg> mcrg_;
    const int mcrg_it_depth; //to which depth the mcrg is done
    const bool measure_structure_factor;
    std::unique_ptr<structure_factor> structure_factor_;

    inline void update(){update(random_int(num_sites()));}
    void update(int site){
        //propose a new state
        double old_state=spins[site];
        double old_energy=single_site_Energy(site); 
        double new_state=old_state+random_real_shifted(angle_dev);
        spins[site]=new_state;
        double new_energy=single_site_Energy(site);
        if(random_real()<=std::exp(-(new_energy-old_energy)/T)){
            //update variables due to local change
            mx+=std::cos(new_state)-std::cos(old_state);
            my+=std::sin(new_state)-std::sin(old_state);

            //update with the corresponding prefactor
            int prefactor_x=1;
            int prefactor_y=1;
            if(((site%L)%2)) prefactor_x=-1; //for even y sites -1
            if(((site/L)%2)) prefactor_y=-1; //for even x sites -1
            mx_stag+=prefactor_x*(std::cos(new_state)-std::cos(old_state));
            my_stag+=prefactor_y*(std::sin(new_state)-std::sin(old_state));
            En+=(new_energy-old_energy)/2;
            ++accepted;
        }
        else{ //switch back
            spins[site]=old_state;
        }
    }
    //this updates the range for the update
    double update_angle_deviation(int accepted, int n_steps){
        double acc_ratio=accepted*1./n_steps;
        if(targeted_acc_ratio>acc_ratio && angle_dev>1e-3*M_PI)
            angle_dev/=1.2;
        if(targeted_acc_ratio<acc_ratio && angle_dev<2*M_PI)
            angle_dev*=1.2;
        accepted=0;
        return acc_ratio; 
    }
    void measure(alps::ObservableSet& obs){
        obs["Energy"]<<En/num_sites();
        obs["Energy^2"]<<std::pow(En/num_sites(),2);
        double M= std::sqrt(mx*mx+my*my)/num_sites();
        obs["M"]<<M;
        obs["M^2"]<<M*M;
        obs["M^4"]<<M*M*M*M;
        double Mx=mx/num_sites();
        obs["Mx"]<<Mx; 
        obs["Mx^2"]<<Mx*Mx;
        M= std::sqrt(mx_stag*mx_stag+my_stag*my_stag)/num_sites();
        obs["M staggered"]<<M;
        obs["M staggered^2"]<<M*M;
        obs["M staggered^4"]<<M*M*M*M;
        Mx=mx_stag/num_sites();
        obs["Mx staggered"]<<Mx; 
        obs["Mx staggered^2"]<<Mx*Mx;
        if(measure_mcrg) 
            mcrg_->measure(spins, obs);
        if(measure_structure_factor) 
            structure_factor_->measure(spins, obs);
    }
    inline void init_spins(const alps::Parameters& params){
        init_type=(params.value_or_default("Initialization","GS")=="Random" ? init_t::Random : init_t::GS);
        spins.resize(N, 0.);
        switch(init_type) {
            case init_t::GS:
                if(params["LATTICE"]=="square lattice")
                    for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter){
                        if(((*s_iter)%2)){//odd y site
                            spins[*s_iter]=M_PI;
                        }
                    }
                else {
                    std::cerr <<"Ground state not explicitly defined, for this lattice, implement this or assume all spins to point in the x direction"<<std::endl;
                }
                break;
            case init_t::Random:
                for(auto& s: spins) {
                    s=random_real(0,2*M_PI);
                }
                break;
            default:
                std::cerr << "Smth went terribly wrong as this line should never be hit"<<std::endl;
        }
            

    }
    inline double beta() const { 
        return 1./T;
    }
    inline int random_int(int j){
        return static_cast<int>(j*uniform_01());
    }
    inline double random_real(double a=0,double b=1){
        return a+(b-a)*uniform_01();
    }
    inline double random_real_shifted(double s){
        return -s/2+s*uniform_01();
    }
    double Energy() const {
        return HamiltonianList.Energy(spins);
    }
    double single_site_Energy(int i) const{
        return HamiltonianList.SingleSiteEnergy(spins, i);
    }
};


class xy_evaluator : public alps::parapack::simple_evaluator {
private:
    double T;
    int L;
    int N;
    double beta() const {
        return 1./T;
    }
public:
    xy_evaluator(alps::Parameters const& params) : 
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        L(params["L"]),
        N(L*L)
        {
        }
    void evaluate(alps::ObservableSet& obs) const {
        // Binder cumulant
        if(obs.has("M^4")&&obs.has("M^2")){
            alps::RealObsevaluator m2 = obs["M^2"];
            alps::RealObsevaluator m4 = obs["M^4"];
            alps::RealObsevaluator binder("BinderCumulant"); 
            binder = m4/(m2*m2);
            obs.addObservable(binder); 
        } else std::cerr << "Binder cumulant will not be calculated"<<std::endl;
        if(obs.has("M staggered^4")&&obs.has("M staggered^2")){
            alps::RealObsevaluator m2 = obs["M staggered^2"];
            alps::RealObsevaluator m4 = obs["M staggered^4"];
            alps::RealObsevaluator binder("BinderCumulant staggered"); 
            binder = 1.-m4/(m2*m2*3);
            obs.addObservable(binder); 
        } else std::cerr << "Binder stag cumulant will not be calculated"<<std::endl;
        // c_V  
        if(obs.has("Energy")&&obs.has("Energy^2")){
            alps::RealObsevaluator E = obs["Energy"];
            alps::RealObsevaluator E2 = obs["Energy^2"];
            alps::RealObsevaluator c_V("c_V");
            c_V = beta()*beta() * (E2-E*E) * N;

            obs.addObservable(c_V); 
        } else std::cerr << "c_V will not be calculated"<<std::endl;
        // susceptibility 
        if(obs.has("Mx")&&obs.has("Mx^2")){
            alps::RealObsevaluator Mx = obs["Mx"];
            alps::RealObsevaluator Mx2 = obs["Mx^2"];
            alps::RealObsevaluator chi("susceptibility");
            chi = beta() * (Mx2-Mx*Mx);
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility will not be calculated"<<std::endl;
        if(obs.has("M staggered")&&obs.has("M staggered^2")){
            alps::RealObsevaluator M = obs["M staggered"];
            alps::RealObsevaluator M2 = obs["M staggered^2"];
            alps::RealObsevaluator chi("susceptibility staggered");
            chi = N*beta() * (M2-M*M); 
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility staggered will not be calculated"<<std::endl;
    }

};
#endif //MCPP_MODELS_XY_H
