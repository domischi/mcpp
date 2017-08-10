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
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm> //sort
#include <utility> //pair

#include "../hamiltonians/hamiltonian_list.h"
#include "../observables/observables.h"
#include "../utilities.h"
class xy_worker : public alps::parapack::lattice_mc_worker<>{
public :
    enum init_t {GS, Random, Ferro, Vortex}; 

    xy_worker(const alps::Parameters& params) :
        alps::parapack::lattice_mc_worker<>(params), 
        is_exmc(static_cast<bool>(params["ALGORITHM"]!="xy")),
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",100)),
        Each_Measurement(params.value_or_default("Each_Measurement",10)),
        Sweeps(static_cast<int>(params.value_or_default("SWEEPS",5000))),
        L(params["L"]),
        N(num_sites()),
        T(mcpp::init_T(params)),
        ising(static_cast<bool>(params.value_or_default("Ising", false))),
        Step_Number(0),
        accepted(0),
        measure_last_configuration(static_cast<bool>(params.value_or_default("measure last configuration",false))),
        HamiltonianList(get_hl(params)),
        observables(construct_observables(params, HamiltonianList)),
        targeted_acc_ratio(params.value_or_default("Targeted Acceptance Ratio",0.5)),
        angle_dev(    params.value_or_default("Angle Deviation Start", 0.1*M_PI)),
        angle_dev_fac(params.value_or_default("Angle Deviation Factor",1.2)),
        angle_dev_min(params.value_or_default("Angle Deviation Min",   1e-3*M_PI)),
        angle_dev_max(params.value_or_default("Angle Deviation Max",   2*M_PI)),
        shape_anisotropy_p(params.value_or_default("Shape Anisotropy p",-1)),
        print_debug_information(static_cast<bool>(params.value_or_default("print debug information",false)))
        {
            if(Sweeps<Each_Measurement){
                std::cerr<< "Less Sweeps than measurement frequency set, this does not make sense. Abort..."<<std::endl;
                std::exit(4);
            }
            init_spins(params);
            En=Energy();
            if(print_debug_information)
                print_debug();
        }
    
    void init_observables(alps::Parameters const&, alps::ObservableSet& obs){
        obs << alps::RealObservable("Energy");
        obs << alps::RealObservable("Energy^2");
        obs << alps::RealObservable("Acceptance Ratio"); //Probably very useful for debugging
        if(measure_last_configuration){
            obs<<alps::RealVectorObservable("Last Configuration");
        }
        for (auto& o : observables)
            o->init_observables(obs);
    }
    void save(alps::ODump &dump) const{
        dump 
        << Step_Number 
        << spins 
        << angle_dev
        << En
        << accepted;
        for(auto & o: observables)
            o->save(dump);
    }
    void load(alps::IDump &dump){
        dump 
        >> Step_Number 
        >> spins 
        >> angle_dev
        >> En
        >> accepted
        ;
        for(auto& o : observables)
            o->load(dump);
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
        return Step_Number > Thermalization_Sweeps;
    }
    double progress() const {
        return Step_Number/(static_cast<double>(Sweeps+Thermalization_Sweeps));
    }
    
    //for exmc
    typedef double weight_parameter_type;
    void set_beta(double beta_) { beta=beta_; T=1./beta_;}
    weight_parameter_type weight_parameter() const {return -Energy();} 
    static double log_weight(weight_parameter_type gw, double beta) {return beta*gw;}
private:
    //System parameters
    const bool is_exmc;
    const int Thermalization_Sweeps;
    const int Each_Measurement;
    const int Sweeps;
    int Step_Number;
    const int L;
    const int N; //L^2
    bool ising;
    init_t init_type;
    std::vector<double> spins; //saves the spins
    const bool print_debug_information;

    double T, beta; // non-const because of exmc

    const double targeted_acc_ratio;
    double angle_dev;
    const double angle_dev_min, angle_dev_max, angle_dev_fac;
    const int shape_anisotropy_p;

    // struct of Hamiltonian pointers -> allows CRTP
    std::shared_ptr<Hamiltonian_List> HamiltonianList;
    std::vector<std::shared_ptr<observable>> observables;

    //Easy observables
    double En;
    int accepted;
    const bool measure_last_configuration;

    inline void update(){update(random_int(num_sites()));}
    void update(int site){
        //propose a new state
        double old_state=spins[site];
        double old_energy=single_site_Energy(site);
        double new_state=mod2Pi(old_state+random_real_shifted(angle_dev));
        if(shape_anisotropy_p>0) {
            new_state=mod2Pi(new_state+random_int(shape_anisotropy_p)*2*M_PI/shape_anisotropy_p); 
        }
        if(ising) {
            new_state= (old_state>0.5*M_PI ? 0 : M_PI); //Only allow Spin Flip updates
        }
        spins[site]=new_state;
        double new_energy=single_site_Energy(site);
        if(random_real()<=std::exp(-(new_energy-old_energy)/T)){
            En+=(new_energy-old_energy);
            ++accepted;
        }
        else{ //switch back
            spins[site]=old_state;
        }
    }

    std::string pprint_vector(vector_type v){
        std::stringstream ret;
        ret<<"(";
        for(int d=0;d<dimension()-1;++d)
            ret<<std::setw(8)<<v[d]<<",";
        ret<<std::setw(8)<<v[dimension()-1];
        ret<<")";
        return ret.str();
    }

    void print_debug(){
        std::vector<vector_type> basis;
        basis_vector_iterator v, v_end;
        for(std::tie(v,v_end)=basis_vectors();v!= v_end;++v){
            basis.push_back(*v);
        }
        std::vector<vector_type> periodic_translations;
        int dim = dimension();
        assert(dim==basis[0].size());
        periodic_translations.push_back(vector_type(dim,0.));
        vector_type pmb(dimension()),ppb(dimension()); //periodic_translation+-L*basis
        for(auto& actual_basis_vector : basis){
            int size_periodic_translations=periodic_translations.size();
            for(int i=0;i<size_periodic_translations;++i){//to avoid double counting a specific vector, and to not have a segfault
                vector_type p = periodic_translations[i];
                for(int d =0;d<dimension();++d){
                    pmb[d]=(p[d]-L*(actual_basis_vector[d]));
                    ppb[d]=(p[d]+L*(actual_basis_vector[d]));
                }
                periodic_translations.push_back(pmb);
                periodic_translations.push_back(ppb);
            }
        }
        for(int i=0;i<periodic_translations.size();++i){
            std::cout <<"Periodic Translation No " <<std::setw(2) <<i<<" with vector "<<pprint_vector(periodic_translations[i]) <<std::endl;
        }
        for(site_iterator s_iter = sites().first; s_iter !=sites().second; ++s_iter){
            vector_type c(coordinate(*s_iter));
            std::cout <<"Site: " <<std::setw(3) <<*s_iter<<" with vector "<<pprint_vector(c);
            std::cout<<std::endl;
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
        return Step_Number>Sweeps+Thermalization_Sweeps-Each_Measurement;
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
        for(auto& o: observables) 
            o->measure(spins, obs);
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
    //inline double beta() const { 
    //    return 1./T;
    //}

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
        return HamiltonianList->Energy(spins);
    }
    double single_site_Energy(int i) const{
        return HamiltonianList->SingleSiteEnergy(spins, i);
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
    std::shared_ptr<field_histogram_evaluator> fhe;
public:
    xy_evaluator(alps::Parameters const& params) : 
        T(mcpp::init_T(params)),
        L(params["L"]),
        N(mcpp::init_N(params))
        {
            if(params.value_or_default("Field Histogram",false)) {
                fhe=std::make_shared<field_histogram_evaluator>(params);
            }
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
            c_V = beta()*beta() * (E2-E*E) * 4 * N;

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
        if(fhe)
            fhe->evaluate(obs);
    }

};
#endif //MCPP_MODELS_XY_H
