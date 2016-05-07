#ifndef MCPP_SSF_H
#define MCPP_SSF_H

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris.h>
#include <alps/osiris/dump.h>
#include <alps/expression.h>
#include <alps/lattice.h>

#include <iostream>
#include <vector>
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
        mx(1.),
        my(0.) {
            //Initialize Spins and the local observables
            spins.resize(N, 0.); 
            
            //coordinate(site_iterator) is a function provided by the graph_helper, 
            //this class is inherited by it, therefore it is here usable,
            //it returns the coordinate as a vector<double> type 
            std::cout << "lattice vectors:" <<std::endl<< /*graph_helper*/coordinate(site(0)).size()<<std::endl;
            //Initialize the measurements
            measurements << alps::RealObservable("Energy");
            measurements << alps::RealObservable("|Magnetization|");
            measurements << alps::RealObservable("Mx");
            
            //TODO get the distance right and implement an angle function between the two points
            //TODO make a table with r_ij^3
            Init_Lookup_Tables(); 
            //TODO also move the neighbour list generation in this lookup_table function 
            //Initialize the neighbour list
            for(int i=0;i<N;++i){
                std::vector<int> tmp;
                for(int j=0;j<N;++j){
                    if(i!=j && distance(i,j)<cutoff_distance) tmp.push_back(j);
                    //std::cout << "i="<<std::setw(4)<<i<<"\tj="<<std::setw(4)<<j<<"\tdist(i,j)="<<std::setw(4)<<distance(site(i),site(j))<<std::endl;
                }
                neighbour_list.push_back(tmp);
            }

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

    typedef std::pair<int,int> spin_pair;

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

    void Init_Lookup_Tables(){//TODO implement this in a general fashion, this is just rediculous....
        for(site_iterator s_iter= sites().first; s_iter!=sites().second; ++s_iter){
            std::cout << "The coods of site "<<*s_iter<<" are: (" <<coordinate(*s_iter)[0]<<","<<coordinate(*s_iter)[1]<<")"<<std::endl;
        }
    }

    double distance(int i, int j){
        return 1.;
    }

    double angle_w_x(int i, int j){
        return 0.;
    }

    void measure(){
        measurements["Energy"]<<En/num_sites();    
        measurements["|Magnetization|"]<<std::sqrt(mx*mx+my*my)/num_sites();    
        measurements["Mx"]<<mx/num_sites();    
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
            e-=D/std::pow(distance(i,j),3)*(1.5*std::cos(spins[i]+spins[j]-2*angle_w_x(i,j))+0.5*std::cos(spins[i]-spins[j]));
        }
        return e;
    }
};

typedef alps::scheduler::SimpleMCFactory<ssf> ssf_factory;
#endif //MCPP_SSF_H
