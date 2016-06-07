#ifndef MCPP_SSF_H
#define MCPP_SSF_H

#include <alps/scheduler.h>
#include <alps/parapack/worker.h>
#include <alps/alea.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris.h>
#include <alps/osiris/dump.h>
#include <alps/expression.h>
#include <alps/lattice.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm> //sort
#include <utility> //pair
//#include "../spins/xy.h"

class ssf_worker : public alps::parapack::lattice_mc_worker<>{
public :
    //ssf_worker(const alps::ProcessList & where,const alps::Parameters& params, int node) :
    ssf_worker(const alps::Parameters& params) :
        alps::parapack::lattice_mc_worker<>(params), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",100)),
        Measure_Sweeps(params.value_or_default("SWEEPS",5000)),
        Each_Measurement(params.value_or_default("Each_Measurement",1)),
        L(params["L"]),
        N(L*L),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        D(params.value_or_default("D",1.)),
        Step_Number(0),
        mx(0.),
        my(0.),
        mx_stag(num_sites()),//stagered in stripes
        my_stag(0.),
        accepted(0){
            cutoff_distance=params.value_or_default("cutoff_distance",3.);
            double a=params.value_or_default("a",1.);
            cutoff_distance*=a;
            //Initialize Spins and the local observables
            spins.resize(N, 0.);
            if(is_bipartite()&&true)//TODO implement check if ground state is striped as used below
                for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter){
                    if((((*s_iter)/L)%2)){//odd y site
                        //spins[*s_iter]=1.5*M_PI;
                        spins[*s_iter]=M_PI;
                    }
                    //else{
                    //    spins[*s_iter]=M_PI/2;
                    //}
                }
            else {
                std::cerr <<"Ground state not explicitly defined, for this lattice, implement this or assume all spins to point in the x direction";
            }
            //Init the lookup tables
            Init_Lookup_Tables(); 

            En=Energy();
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
    }

    static void print_copyright(std::ostream & out){
        out << "You are using mc++"<<std::endl
            << "copyright (c) by Dominik Schildknecht"<<std::endl
            << "if you reuse this project, please mention the ALPS project and me as a fair user"<<std::endl;
    }

    void save(alps::ODump &dump) const{
        dump << L << N<< T<< Step_Number << spins;
    }
    void load(alps::IDump &dump){
        dump >> L >> N>> T>> Step_Number >> spins;
    }
    void run(alps::ObservableSet& obs){
        ++Step_Number;
        for(int i = 0;i<N;++i){
            update();
        }
        if(!(Step_Number%Each_Measurement) && Step_Number>Thermalization_Sweeps){
            measure(obs);
            accepted=0;
        }
    }
    bool is_thermalized() const {
        return Step_Number >= Thermalization_Sweeps;   
    }
    double progress() const {
        return Step_Number/(static_cast<double>(Measure_Sweeps+Thermalization_Sweeps));
    }
private:
    int Thermalization_Sweeps;
    int Measure_Sweeps;
    int Each_Measurement; //do Each_Measurement steps, then measure
    int Step_Number;
    int L;
    int N; //L^2
    std::vector<double> spins; //saves the spins
    double T;
    double D;
    double cutoff_distance;

    std::map<std::pair<int,int>,double> dist_3;
    std::map<std::pair<int,int>,double> phi;

    std::vector<std::vector<int>> neighbour_list; //saves which neighbours are relevant (as not only nearest neighbours count in dipole )

    double En;
    double mx;
    double my;
    double mx_stag;
    double my_stag;
    int accepted;
    
    void update(){
        //Choose a site
        int site=random_int(num_sites());
        //int site=L;
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
            if(((site/L)%2)) prefactor_x=-1; //for even y sites -1
            if(((site%L)%2)) prefactor_y=-1; //for even x sites -1
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
        obs["Acceptance Rate"] << (1.0*accepted)/Each_Measurement;
    }
    //TODO decide in which units of measurement to measure...(natural, eV and K, SI?)
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
                    phi[pair_]=std::atan2((c1[1]+p[1]-c2[1]),(c1[0]+p[0]-c2[0]));
                    phi[pair_inverse]=std::atan2(-(c1[1]+p[1]-c2[1]),-(c1[0]+p[0]-c2[0]));
                    neighbour_list[*s_iter].push_back(*s_iter2);
                    neighbour_list[*s_iter2].push_back(*s_iter);
                }
            }
    }

    inline double distance(vector_type& x, vector_type& y, vector_type& periodic){
        return std::sqrt(std::pow(x[0]-y[0]+periodic[0],2)+std::pow(x[1]-y[1]+periodic[1],2));
    }

    double inv_distance_cubed(int i,int j) {return inv_distance_cubed(std::make_pair(i,j));}

    double inv_distance_cubed(std::pair<int,int> pair_){
        double ret_val=dist_3[pair_];
        assert(ret_val>0.);
        return ret_val;
    }

    double angle_w_x(int i, int j) {return angle_w_x(std::make_pair(i,j));}

    double angle_w_x(std::pair<int,int> pair_){
        return phi[pair_];
    }

    int random_int(int j){
        return static_cast<int>(j*uniform_01());
    }
    double random_real(double a=0,double b=1){
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
            e+=D*inv_distance_cubed(i,j)*(1.5*std::cos(spins[i]+spins[j]-2*angle_w_x(i,j))+0.5*std::cos(spins[i]-spins[j]));
        }
        return e;
    }
};


class ssf_evaluator : public alps::parapack::simple_evaluator {
private:
    double T;
    int L;
    int N;
    double beta() const {
        return 1./T;
    }
public:
    ssf_evaluator(alps::Parameters const& params) : T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),L(params["L"]),N(L*L) {}
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
            binder = m4/(m2*m2);
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
            chi = beta()* (Mx2-Mx*Mx); //TODO divide by num_sites()?
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility will not be calculated"<<std::endl;
        if(obs.has("M staggered")&&obs.has("M staggered^2")){
            alps::RealObsevaluator Mx = obs["M staggered"];
            alps::RealObsevaluator Mx2 = obs["M staggered^2"];
            alps::RealObsevaluator chi("susceptibility staggered");
            chi = beta() * (Mx2-Mx*Mx) * N; 
            obs.addObservable(chi); 
        } else std::cerr << "susceptibility staggered will not be calculated"<<std::endl;
    }

};
#endif //MCPP_SSF_H
