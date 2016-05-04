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

//#include "../spins/xy.h"

class ssf : public alps::scheduler::LatticeMCRun<>{
public :
    ssf(const alps::ProcessList & where,const alps::Parameters& params, int node) :
        alps::scheduler::LatticeMCRun<>(where,params,node), 
        Thermalization_Sweeps(params.value_or_default("THERMALIZATION",10000)),
        Measure_Sweeps(params.value_or_default("SWEEPS",50000)),
        L(params["L"]),
        N(L*L),
        T(params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"])),
        Step_Number(0) {
           spins.resize(N, 0.); 
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
        if(Step_Number%Each_Measurement){
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

    void update(){}
    void measure(){}
};

typedef alps::scheduler::SimpleMCFactory<ssf> ssf_factory;
#endif //MCPP_SSF_H
