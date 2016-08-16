#ifndef MCPP_SSF_H
#define MCPP_SSF_H

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
#include <map>
#include <cassert>
#include <algorithm> //sort
#include <utility> //pair

#include "../special-observables/mcrg.h"

class ssf_worker : public alps::parapack::lattice_mc_worker<>{
public :
    ssf_worker(const alps::Parameters& params) :
        alps::parapack::lattice_mc_worker<>(params), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",100)),
        Measure_Sweeps(params.value_or_default("SWEEPS",5000)),
        L(params["L"]),
        N(L*L),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        D(params.value_or_default("D",1.)),
        Step_Number(0),
        mx(0.),
        my(0.),
        mx_stag(num_sites()),//stagered in stripes
        my_stag(0.),
        accepted(0),
        //ac_measured(false),
        //safety_factor(params.value_or_default("safety_factor",3)),
        //autocorrelation(0),
        //ac_N(0),
        //ac_obs("autocorr obs"),
        mcrg_it_depth(params.value_or_default("mcrg_iteration_depth",-1)),
        Each_Measurement(params.value_or_default("Each_Measurement",15))
        {
            measure_mcrg=(mcrg_it_depth>0);
            cutoff_distance=params.value_or_default("cutoff_distance",3.);
            double a=params.value_or_default("a",1.);
            cutoff_distance*=a;
            //Initialize Spins and the local observables
            spins.resize(N, 0.);
            if(is_bipartite()&&true)//TODO implement check if ground state is striped as used below
                for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter){
                    if(((*s_iter)%2)){//odd y site
                        spins[*s_iter]=M_PI;
                    }
                }
            else {
                std::cerr <<"Ground state not explicitly defined, for this lattice, implement this or assume all spins to point in the x direction"<<std::endl;
            }
            //Init the lookup tables
            Init_Lookup_Tables(); 

            En=Energy();
            //if(T<=1e-3){ //at T=0 there is a problem with this, as the system doesn't move at all (0/0 problem) 
            //    ac_measured=true;
            //    autocorrelation=1;
            //}
            if(measure_mcrg){
                std::cout << "\tInitialize MCRG with iteration depth "<<mcrg_it_depth<<"..."<<std::flush;
                mcrg_=std::make_shared<mcrg>(params,0,mcrg_it_depth);
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
            obs << alps::RealObservable("Acceptance Rate"); //Probably very useful for debugging
            if(measure_mcrg){
                mcrg_->init_observables(obs);
            }
    }

    //static void print_copyright(std::ostream & out){
    //    out << "You are using mc++"<<std::endl
    //        << "copyright (c) by Dominik Schildknecht"<<std::endl
    //        << "if you reuse this project, please mention the ALPS project and me as a fair user"<<std::endl;
    //}

    void save(alps::ODump &dump) const{
        dump << L << N<< T<< Step_Number << spins;
    }
    void load(alps::IDump &dump){
        dump >> L >> N>> T>> Step_Number >> spins;
    }
    void run(alps::ObservableSet& obs){
        using namespace alps::alea;
        ++Step_Number;
        //int MAX_STEPS_AC = 1e6;
        //int n_steps = ac_measured ? autocorrelation*safety_factor+1 : 0;
        int n_steps = N;
        for(int i = 0;i<n_steps;++i){// Do a typewriter sweep (possibly avoiding cache reload every step)
            update(i);
        }
        if(Step_Number&(1<<10)){ //Every 1024 steps do a random site lattice sweep to avoid ergodicity problems at 0K
            for(int i = 0;i<n_steps;++i){
                update();
            }
        }
        //if(!ac_measured && is_thermalized()) { //the ac_time still needs to be determined
        //    for(int i = 0;i<MAX_STEPS_AC;++i){
        //        update();
        //        ac_obs<<En;
        //        if(!(i % 512)) { //only with full bins the ac time can be measured
        //            measure_ac_time();
        //            if(ac_measured) break;
        //        }
        //        if(i==MAX_STEPS_AC-1) {
        //            std::cerr << "didn't find any meaningful autocorrelation time!!"<<std::endl;
        //            std::exit(13);
        //        }
        //    }
        //}
        //if(ac_measured && is_thermalized()){
        if(is_thermalized()&&(Step_Number%Each_Measurement)){
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
    int Thermalization_Sweeps;
    int Measure_Sweeps;
    int Step_Number;
    int L;
    int N; //L^2
    std::vector<double> spins; //saves the spins
    double T;
    double D;
    double cutoff_distance;
    int Each_Measurement;
    ////autocorrelation time parameters
    //alps::RealObservable ac_obs;
    //bool ac_measured;
    //int ac_N;
    //double safety_factor;
    //double autocorrelation;

    //Lookup tables 
    std::map<std::pair<int,int>,double> dist_3;
    std::map<std::pair<int,int>,double> phi;
    std::vector<std::vector<int>> neighbour_list; //saves which neighbours are relevant (as not only nearest neighbours count in dipole )

    //Easy observables
    double En;
    double mx;
    double my;
    double mx_stag;
    double my_stag;
    int accepted;
    //Complex observables
    bool measure_mcrg;
    std::shared_ptr<mcrg> mcrg_;
    int mcrg_it_depth; //to which depth the mcrg is done

    //double measure_ac_time(){
    //    if(is_thermalized()&&!ac_measured) {
    //        if(std::isfinite(N*ac_obs.tau()) && N*ac_obs.tau()>0&&N*ac_obs.tau()>autocorrelation){ 
    //            autocorrelation=N*ac_obs.tau();
    //            ac_N=0;
    //        }
    //        else
    //            ++ac_N;

    //        if(ac_N>=10){//autocorrelation time is stable 
    //            ac_measured=true;
    //        }
    //        ac_obs.reset(true);
    //    }
    //    return autocorrelation;
    //}

    void update(){update(random_int(num_sites()));}
    void update(int site){
        //propose a new state
        double new_state=random_real(0.,2*M_PI);
        double old_energy=single_site_Energy(site);
        double old_state=spins[site];
        
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
        //obs["Acceptance Rate"] << (1.0*accepted)/(autocorrelation*safety_factor+1);
        if(measure_mcrg) 
            mcrg_->measure(spins, obs);
    }
    inline double beta(){ 
        return 1./T;
    }
    void Init_Lookup_Tables(){ 
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
        for(int i=0;i<num_sites();++i) neighbour_list.push_back(std::vector<int>());
        std::map<std::pair<int,int>,double> dist_map;
        for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter)
        for(site_iterator s_iter2= s_iter; s_iter2!=sites().second; ++s_iter2)
        for(auto& p : periodic_translations)
            if(*s_iter!=*s_iter2){
                vector_type c1(coordinate(*s_iter));
                vector_type c2(coordinate(*s_iter2));
                double dist=distance(c1,c2,p);
                std::pair<int,int> pair_ = std::make_pair(*s_iter,*s_iter2);
                if(dist<=cutoff_distance && (dist_map[pair_]==0. || dist<dist_map[pair_])){ //new value is smaller than the previous calculated for this pair
                    dist_map[pair_]=dist;
                    std::pair<int,int> pair_inverse=std::make_pair(pair_.second,pair_.first);
                    dist_3[pair_]=std::pow(dist,-3);
                    dist_3[pair_inverse]=std::pow(dist,-3);
                    phi[pair_]=std::atan2(-(c1[1]+p[1]-c2[1]),-(c1[0]+p[0]-c2[0]));
                    phi[pair_inverse]=std::atan2((c1[1]+p[1]-c2[1]),(c1[0]+p[0]-c2[0]));
                    neighbour_list[*s_iter].push_back(*s_iter2);
                    neighbour_list[*s_iter2].push_back(*s_iter);
                }
            }
        //Shrink the neighbour_list
        for(auto& nl : neighbour_list) nl.shrink_to_fit();
        neighbour_list.shrink_to_fit();
    }
    inline double distance(vector_type& x, vector_type& y, vector_type& periodic){
        return std::sqrt(std::pow(x[0]-y[0]+periodic[0],2)+std::pow(x[1]-y[1]+periodic[1],2));
    }
    inline double inv_distance_cubed(int i,int j) {
        return inv_distance_cubed(std::make_pair(i,j));
    }
    inline double inv_distance_cubed(std::pair<int,int> pair_){
        double ret_val=dist_3[pair_];
        assert(ret_val>0.);
        return ret_val;
    }
    inline double angle_w_x(int i, int j) { 
        return angle_w_x(std::make_pair(i,j));
    }
    inline double angle_w_x(std::pair<int,int> pair_){
        return phi[pair_];
    }
    inline int random_int(int j){
        return static_cast<int>(j*uniform_01());
    }
    inline double random_real(double a=0,double b=1){
        return a+(b-a)*uniform_01();
    }
    double Energy(){
        double e=0.;
        for(int i =0;i<N;++i){
            e+=single_site_Energy(i);
        }
        return e/2;
    }
    double single_site_Energy(int i){
        double e=0.;
        for(int j : neighbour_list[i]){
            e-=D*inv_distance_cubed(i,j)*(1.5*std::cos(spins[i]+spins[j]-2*angle_w_x(i,j))+0.5*std::cos(spins[i]-spins[j]));
        }
        return e;
    }
};


class ssf_evaluator : public alps::parapack::simple_evaluator {
private:
    double T;
    int L;
    int N;
    int mcrg_it_depth;
    bool measure_mcrg;
    double beta() const {
        return 1./T;
    }
public:
    ssf_evaluator(alps::Parameters const& params) : 
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        L(params["L"]),
        N(L*L),
        mcrg_it_depth(params.value_or_default("mcrg_iteration_depth",-1)),
        measure_mcrg(mcrg_it_depth>0)
        {}
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
            chi = beta() * (Mx2-Mx*Mx); //TODO divide by num_sites()?
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility will not be calculated"<<std::endl;
        if(obs.has("M staggered")&&obs.has("M staggered^2")){
            alps::RealObsevaluator M = obs["M staggered"];
            alps::RealObsevaluator M2 = obs["M staggered^2"];
            alps::RealObsevaluator chi("susceptibility staggered");
            chi = N*beta() * (M2-M*M); 
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility staggered will not be calculated"<<std::endl;
        //if(measure_mcrg){
        //    for(int it=1;it<=mcrg_it_depth;++it){
        //        if (!(obs.has("MCRG S_alpha even " + std::to_string(it)) &&
        //              obs.has("MCRG S_alpha even "+ std::to_string(it-1)) &&
        //              obs.has("MCRG S_alpha S_beta same iteration even "+ std::to_string(it)) &&
        //              obs.has("MCRG S_alpha S_beta next iteration even "+ std::to_string(it)) &&
        //              obs.has("MCRG S_alpha odd " + std::to_string(it)) &&
        //              obs.has("MCRG S_alpha odd "+ std::to_string(it-1)) &&
        //              obs.has("MCRG S_alpha S_beta same iteration odd "+ std::to_string(it)) &&
        //              obs.has("MCRG S_alpha S_beta next iteration odd "+ std::to_string(it)))) 
        //            std::cerr << "NOT ALL INFORMATION FOR A MCRG ANALYSIS WAS GIVEN"<<std::endl;
        //        else{ //calculate the dS^(n)/dK^(n-1) and dS^(n)/dK^(n)
        //            //even
        //            std::valarray<double> S_a_n_e = alps::RealVectorObsevaluator(obs["MCRG S_alpha even "+std::to_string(it)]).mean();
        //            std::valarray<double> S_a_nm1_e = alps::RealVectorObsevaluator(obs["MCRG S_alpha even "+std::to_string(it-1)]).mean();
        //            std::valarray<double> S_a_n_S_b_n_e = alps::RealVectorObsevaluator(obs["MCRG S_alpha S_beta same iteration even "+std::to_string(it)]).mean(); 
        //            std::valarray<double> S_a_n_S_b_nm1_e = alps::RealVectorObsevaluator(obs["MCRG S_alpha S_beta next iteration even "+std::to_string(it)]).mean();
        //            std::vector<double> m_S_a_m_S_b_st_e;
        //            std::vector<double> m_S_a_m_S_b_nt_e;
        //            int n_alpha_e=mcrg::n_interactions_even();
        //            for(int i=0;i<n_alpha_e;++i)
        //                for(int j=0;j<n_alpha_e;++j){
        //                    m_S_a_m_S_b_st_e.push_back((S_a_n_e[i]*S_a_n_e[j]));
        //                    m_S_a_m_S_b_nt_e.push_back((S_a_nm1_e[i]*S_a_nm1_e[j]));
        //                }
        //            std::valarray<double> mSamSb_st_e(m_S_a_m_S_b_st_e.data(),m_S_a_m_S_b_st_e.size());
        //            std::valarray<double> mSamSb_nt_e(m_S_a_m_S_b_nt_e.data(),m_S_a_m_S_b_nt_e.size());
        //            //odd
        //            std::valarray<double> S_a_n_o = alps::RealVectorObsevaluator(obs["MCRG S_alpha odd "+std::to_string(it)]).mean();
        //            std::valarray<double> S_a_nm1_o = alps::RealVectorObsevaluator(obs["MCRG S_alpha odd "+std::to_string(it-1)]).mean();
        //            std::valarray<double> S_a_n_S_b_n_o = alps::RealVectorObsevaluator(obs["MCRG S_alpha S_beta same iteration odd "+std::to_string(it)]).mean(); 
        //            std::valarray<double> S_a_n_S_b_nm1_o = alps::RealVectorObsevaluator(obs["MCRG S_alpha S_beta next iteration odd "+std::to_string(it)]).mean();
        //            std::vector<double> m_S_a_m_S_b_st_o;
        //            std::vector<double> m_S_a_m_S_b_nt_o;
        //            int n_alpha_o=mcrg::n_interactions_odd();
        //            for(int i=0;i<n_alpha_o;i+=2)
        //                for(int j=0;j<n_alpha_o;j+=2){
        //                    m_S_a_m_S_b_st_o.push_back((S_a_n_o[i]*S_a_n_o[j]+S_a_n_o[i+1]*S_a_n_o[j+1]));
        //                    m_S_a_m_S_b_nt_o.push_back((S_a_nm1_o[i]*S_a_nm1_o[j]+S_a_nm1_o[i+1]*S_a_nm1_o[j+1]));
        //                }
        //            std::valarray<double> mSamSb_st_o(m_S_a_m_S_b_st_o.data(),m_S_a_m_S_b_st_o.size());
        //            std::valarray<double> mSamSb_nt_o(m_S_a_m_S_b_nt_o.data(),m_S_a_m_S_b_nt_o.size());
        //            if(true){
        //                std::valarray<double> dS_n_dK_n_e = S_a_n_S_b_n_e - mSamSb_st_e;
        //                std::valarray<double> dS_n_dK_nm1_e = S_a_n_S_b_nm1_e - mSamSb_nt_e;
        //                alps::RealVectorObservable dS_n_dK_n_e_obs("MCRG dSdK same iteration even "+ std::to_string(it));
        //                alps::RealVectorObservable dS_n_dK_nm1_e_obs("MCRG dSdK next iteration even "+ std::to_string(it));
        //                dS_n_dK_n_e_obs<<dS_n_dK_n_e;
        //                dS_n_dK_nm1_e_obs<<dS_n_dK_nm1_e;
        //                obs.addObservable(dS_n_dK_n_e_obs);
        //                obs.addObservable(dS_n_dK_nm1_e_obs);
        //                std::valarray<double> dS_n_dK_n_o =S_a_n_S_b_n_o - mSamSb_st_o;
        //                std::valarray<double> dS_n_dK_nm1_o=S_a_n_S_b_nm1_o - mSamSb_nt_o;
        //                alps::RealVectorObservable dS_n_dK_n_o_obs("MCRG dSdK same iteration odd "+ std::to_string(it));
        //                alps::RealVectorObservable dS_n_dK_nm1_o_obs("MCRG dSdK next iteration odd "+ std::to_string(it));
        //                dS_n_dK_n_o_obs<<dS_n_dK_n_o;  
        //                //dS_n_dK_nm1_o_obs<< dS_n_dK_nm1_o;
        //                obs.addObservable(dS_n_dK_n_o_obs);
        //                obs.addObservable(dS_n_dK_nm1_o_obs);
        //            } else 
        //                std::cerr << "NO MEASUREMENTS WERE MADE IN THE OBSERVABLES NEEDED FOR MCRG";
        //        }
        //    }
        //}
    }

};
#endif //MCPP_SSF_H
