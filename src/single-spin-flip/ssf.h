#ifndef MCPP_SSF_H
#define MCPP_SSF_H

#include <alps/scheduler.h>
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
#include <utility> //pair
//#include "../spins/xy.h"

class ssf : public alps::scheduler::LatticeMCRun<>{
public :
    ssf(const alps::ProcessList & where,const alps::Parameters& params, int node) :
        alps::scheduler::LatticeMCRun<>(where,params,node), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",10000)),
        Measure_Sweeps(params.value_or_default("SWEEPS",50000)),
        Each_Measurement(params.value_or_default("Each_Measurement",100)),
        L(params["L"]),
        N(L*L),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        D(params.value_or_default("D",1.)),
        cutoff_distance(params.value_or_default("cutoff_distance",3.)),
        Step_Number(0),
        mx(num_sites()),
        my(0.) {
            //Initialize Spins and the local observables
            spins.resize(N, 0.); 
            
            //Initialize the measurements
            measurements << alps::RealObservable("Energy");
            measurements << alps::RealObservable("M");
            measurements << alps::RealObservable("Mx");
            
            Init_Lookup_Tables(); 

            En=Energy();
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
    void dostep(){
        ++Step_Number;
        update();
        if(Step_Number%Each_Measurement && Step_Number>0){
            measure();
        }
        if(Step_Number==Measure_Sweeps+Thermalization_Sweeps){
            evaluate();
        }
    }
    bool is_thermalized() const {
        return Step_Number >= Thermalization_Sweeps;   
    }
    double work_done() const {
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
    
    void update(){
        //Choose a site
        int site=random_int(num_sites());
     
        //propose a new state
        double new_state=random_real(0.,2*M_PI);
        double old_energy=single_site_Energy(site);
        double old_state=spins[site];

        spins[site]=new_state;
        double new_energy=single_site_Energy(site);
        //std::cout << "old:"<<std::setw(10)<<old_energy<<"\tnew:"<<std::setw(10)<<new_energy<<std::endl;
        if(random_real()<=std::exp(-(new_energy-old_energy)/T)){//check this line, as this is not yet checked
            //update variables due to local change
            mx+=std::cos(new_state)-std::cos(old_state);
            my+=std::sin(new_state)-std::sin(old_state);
            En+=(new_energy-old_energy)/2;
        }
        else{ //switch back
            spins[site]=old_state;
        }
    }

    void measure(){
        measurements["Energy"]<<En/num_sites();
        double M= std::sqrt(mx*mx+my*my)/num_sites();
        measurements["M"]<<M;
        measurements["M^2"]<<M*M;
        measurements["M^4"]<<M*M*M*M;
        double Mx = mx/num_sites();
        measurements["Mx"]<<Mx; 
        measurements["Mx^2"]<<Mx*Mx;  
    }
    //Calculate the properties like the Binder Cumulant, the susceptibility and so on.
    void evaluate(){
        // Binder cumulant
        if(measurements.has("M^4")&&measurements.has("M^2")){
            alps::RealObsevaluator m2 = measurements["M^2"];
            alps::RealObsevaluator m4 = measurements["M^4"];
            alps::RealObsevaluator binder = m2*m2/m4;
            measurements.addObservable(binder); 
        } else std::cout << "Binder cumulant will not be calculated"<<std::endl;
        // c_V 
        if(measurements.has("Energy")&&measurements.has("Energy^2")){
            alps::RealObsevaluator E = measurements["Energy"];
            alps::RealObsevaluator E2 = measurements["Energy^2"];
            alps::RealObsevaluator c_V= beta()*beta() * (E2-E*E); //TODO divide by num_sites()?
            measurements.addObservable(c_V); 
        } else std::cout << "c_V will not be calculated"<<std::endl;
        // susceptibility 
        if(measurements.has("Mx")&&measurements.has("Mx^2")){
            alps::RealObsevaluator Mx = measurements["Mx"];
            alps::RealObsevaluator Mx2 = measurements["Mx^2"];
            alps::RealObsevaluator chi= beta() * (Mx2-Mx*Mx); //TODO divide by num_sites()?
            measurements.addObservable(chi); 
        } else std::cout << "susceptibility will not be calculated"<<std::endl;
    }

    //TODO decide in which units of measurement to measure...(natural, eV and K, SI?)
    inline double beta(){ 
        return 1./T;
    }

    void Init_Lookup_Tables(){ 
        std::cout <<"\tInitialize Lookup Tables..."<<std::flush;
        //TODO i need to find the graph_helper function which returns the vector which is describing the super lattice
        std::vector<vector_type> basis;
        basis_vector_iterator v, v_end;
        for(std::tie(v,v_end)=basis_vectors();v!= v_end;++v){
            basis.push_back(*v);
        }
        std::vector<vector_type> periodic_translations;
        //TODO generalize to more than 2D, and probably also make it nicer (maybe with a valarray)
        for(int i=-1;i<=1;++i)
        for(int j=-1;j<=1;++j){
            vector_type vec, b1, b2;
            b1=basis[0];
            b2=basis[1];
            vec.push_back(b1[0]*L*i+b2[0]*L*j);
            vec.push_back(b1[1]*L*i+b2[1]*L*j);
            periodic_translations.push_back(vec);
        }
        for(int i=0;i<num_sites();++i) neighbour_list.push_back(std::vector<int>());
        std::map<std::pair<int,int>,double> dist_map;
        for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter)
        for(site_iterator s_iter2= sites().first; s_iter2!=sites().second; ++s_iter2)
        for(auto& p : periodic_translations)
            if(*s_iter!=*s_iter2){
                vector_type c1(coordinate(*s_iter));
                vector_type c2(coordinate(*s_iter2));
                double dist=distance(c1,c2,p);
                std::pair<int,int> pair_ = std::make_pair(*s_iter,*s_iter2);
                if(dist<cutoff_distance && (dist_map[pair_]==0. || dist<dist_map[pair_])){ //new value is smaller than the previous calculated for this pair
                    dist_map[pair_]=dist;
                    std::pair<int,int> pair_inverse=std::make_pair(pair_.second,pair_.first);
                    dist_3[pair_]=std::pow(dist,-3);
                    dist_3[pair_inverse]=std::pow(dist,-3);
                    phi[pair_]=std::atan2((c1[1]-c2[1]),(c1[0]-c2[0]));
                    phi[pair_inverse]=std::atan2((c1[1]-c2[1]),(c1[0]-c2[0]));
                    neighbour_list[*s_iter].push_back(*s_iter2);
                    neighbour_list[*s_iter2].push_back(*s_iter);
                }
        }
        std::cout <<"done"<<std::endl;
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

typedef alps::scheduler::SimpleMCFactory<ssf> ssf_factory;
#endif //MCPP_SSF_H
