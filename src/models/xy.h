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
#include "../observables/observables.h"

class xy_worker : public alps::parapack::lattice_mc_worker<>{
public :
    enum init_t {GS, Random, Ferro, Vortex}; 

    xy_worker(const alps::Parameters& params) :
        alps::parapack::lattice_mc_worker<>(params), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",100)),
        Measure_Sweeps(params.value_or_default("SWEEPS",5000)),
        L(params["L"]),
        N(num_sites()),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        HamiltonianList(params),
        ising(static_cast<bool>(params.value_or_default("Ising", false))),
        Step_Number(0),
        accepted(0),
        measure_basic_observables(static_cast<bool>(params.value_or_default("basic_observables",true))),
        measure_last_configuration(static_cast<bool>(params.value_or_default("measure last configuration",false))),
        mcrg_it_depth(params.value_or_default("mcrg_iteration_depth",-1)),
        measure_mcrg(static_cast<bool>(static_cast<int>(params.value_or_default("mcrg_iteration_depth",-1))>0)),
        measure_structure_factor(static_cast<bool>(params.value_or_default("structure_factor",false))),
        measure_spin_autocorrelation(static_cast<bool>(static_cast<int>(params.value_or_default("Spin autocorrelation analysis length",-1))>0)),
        Each_Measurement(params.value_or_default("Each_Measurement",15)),
        targeted_acc_ratio(params.value_or_default("Targeted Acceptance Ratio",0.5)),
        angle_dev(    params.value_or_default("Angle Deviation Start", 0.1*M_PI)),
        angle_dev_fac(params.value_or_default("Angle Deviation Factor",1.2)),
        angle_dev_min(params.value_or_default("Angle Deviation Min",   1e-3*M_PI)),
        angle_dev_max(params.value_or_default("Angle Deviation Max",   2*M_PI))
        {
            init_spins(params);
            En=Energy();
            if(measure_basic_observables){
                basic_observables_=std::unique_ptr<basic_observables>(new basic_observables(params));
            }
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
            if(measure_spin_autocorrelation){
                std::cout << "\tInitialize Spin Autocorrelation Measurement..."<<std::flush;
                spin_autocorrelation_=std::unique_ptr<spin_autocorrelation>(new spin_autocorrelation(params));
                std::cout << "\tdone"<<std::endl;
            }
        }
    
    void init_observables(alps::Parameters const&, alps::ObservableSet& obs){
        obs << alps::RealObservable("Energy");
        obs << alps::RealObservable("Energy^2");
        obs << alps::RealObservable("Acceptance Ratio"); //Probably very useful for debugging
        if(measure_last_configuration){
            obs<<alps::RealVectorObservable("Last Configuration");
        }
        if(measure_mcrg){
            mcrg_->init_observables(obs);
        }
        if(measure_basic_observables){
            basic_observables_->init_observables(obs);
        }
        if(measure_structure_factor){
            structure_factor_->init_observables(obs);
        }
        if(measure_spin_autocorrelation){
            spin_autocorrelation_->init_observables(obs);
        }
    }

    void save(alps::ODump &dump) const{
        dump 
        << Step_Number 
        << spins 
        << angle_dev
        << En
        << accepted;
        if(measure_basic_observables){
            basic_observables_->save(dump);
        }
        if(measure_mcrg){
            mcrg_->save(dump);
        }
        if(measure_structure_factor){
            structure_factor_->save(dump);
        }
        if(measure_spin_autocorrelation){
            spin_autocorrelation_->save(dump);
        }
    }
    void load(alps::IDump &dump){
        dump 
        >> Step_Number 
        >> spins 
        >> angle_dev
        >> En
        >> accepted
        ;
        if(measure_spin_autocorrelation)
            spin_autocorrelation_->load(dump);
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
        if(is_thermalized()&&(!((Step_Number-Thermalization_Sweeps)%Each_Measurement))){
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
    bool ising;
    init_t init_type;
    std::vector<double> spins; //saves the spins
    const double T;
    const int Each_Measurement;
    
    const double targeted_acc_ratio;
    double angle_dev;
    const double angle_dev_min, angle_dev_max, angle_dev_fac;

    // struct of Hamiltonian pointers -> allows CRTP
    Hamiltonian_List HamiltonianList;

    //Easy observables
    double En;
    int accepted;
    const bool measure_last_configuration;
    //Observables
    const bool measure_basic_observables;
    std::shared_ptr<basic_observables> basic_observables_;
    const bool measure_mcrg;
    std::shared_ptr<mcrg> mcrg_;
    const int mcrg_it_depth; //to which depth the mcrg is done
    const bool measure_structure_factor;
    std::unique_ptr<structure_factor> structure_factor_;
    const bool measure_spin_autocorrelation;
    std::unique_ptr<spin_autocorrelation> spin_autocorrelation_;

    inline void update(){update(random_int(num_sites()));}
    void update(int site){
        //propose a new state
        double old_state=spins[site];
        double old_energy=single_site_Energy(site);
        double new_state=mod2Pi(old_state+random_real_shifted(angle_dev));
        if(ising) {
            new_state= (old_state>0.5*M_PI ? 0 : M_PI); //Only allow Spin Flip updates
        }
        spins[site]=new_state;
        double new_energy=single_site_Energy(site);
        if(random_real()<=std::exp(-(new_energy-old_energy)/T)){
            En+=new_energy-old_energy;
            ++accepted;
        }
        else{ //switch back
            spins[site]=old_state;
        }
    }
    //this updates the range for the update
    double update_angle_deviation(int accepted, int n_steps){
        double acc_ratio=accepted*1./n_steps;
        if(targeted_acc_ratio>acc_ratio && angle_dev>angle_dev_min) 
            angle_dev/=angle_dev_fac;
        if(targeted_acc_ratio<acc_ratio && angle_dev<angle_dev_max)
            angle_dev*=angle_dev_fac;
        accepted=0;
        return acc_ratio; 
    }

    inline bool is_last_MC_step(){
        return Step_Number==Measure_Sweeps*Each_Measurement+Thermalization_Sweeps;
    }
    
    double inline mod2Pi(double const& s) const {
        return s-std::floor(s/(2*M_PI))*2*M_PI;
    }
    void measure(alps::ObservableSet& obs){
        obs["Energy"]<<En/num_sites();
        obs["Energy^2"]<<std::pow(En/num_sites(),2);
        if(measure_last_configuration && is_last_MC_step()){
            std::valarray<double> s(spins.data(),spins.size());
            obs["Last Configuration"]<<s;
        }
        if(measure_basic_observables) 
            basic_observables_->measure(spins, obs);
        if(measure_mcrg) 
            mcrg_->measure(spins, obs);
        if(measure_structure_factor)
            structure_factor_->measure(spins, obs);
        if(measure_spin_autocorrelation)
            spin_autocorrelation_->measure(spins, obs);
    }
    inline void init_spins(const alps::Parameters& params){
        if(params.value_or_default("Initialization","GS")=="Random"){
            init_type=init_t::Random;
        } else if (params.value_or_default("Initialization","GS")=="GS") {
            init_type=init_t::GS;
        } else if (params.value_or_default("Initialization","GS")=="Ferro"){
            init_type=init_t::Ferro;
        } else if (params.value_or_default("Initialization","GS")=="Vortex"){
            init_type=init_t::Vortex;
        } else {
            std::cerr<< "Did not recognise the initialization type, typo? Abort now...";
            std::exit(3);
        }
        spins.resize(N, 0.);
        switch(init_type) {
            case init_t::GS:
                if((params["LATTICE"]=="square lattice" || params["LATTICE"]=="anisotropic square lattice") && !(L%2))
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
            case init_t::Ferro: //Already initialized as Ferro due to resize
                break;
            case init_t::Vortex: //Already initialized as Ferro due to resize
                if((params["LATTICE"]=="square lattice" || params["LATTICE"]=="anisotropic square lattice") && !(L%2))
                    for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter){
                        if((*s_iter)%2){//odd y site
                            if((*s_iter/L)%2)//odd x site
                                spins[*s_iter]=7./4*M_PI;
                            else//even x site
                                spins[*s_iter]=1./4*M_PI;
                        }
                        else{//even y site
                            if((*s_iter/L)%2)//odd x site
                                spins[*s_iter]=5./4*M_PI;
                            else//even x site
                                spins[*s_iter]=3./4*M_PI;
                        }
                    }
                else {
                    std::cerr <<"Vortex state not explicitly defined, for this lattice, implement this or assume all spins to point in the x direction"<<std::endl;
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
            binder = 1.-m4/(m2*m2*3);
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
        if(obs.has("M")&&obs.has("M^2")){
            alps::RealObsevaluator M = obs["M"];
            alps::RealObsevaluator M2 = obs["M^2"];
            alps::RealObsevaluator chi("susceptibility");
            chi = N*beta() * (M2-M*M); 
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
