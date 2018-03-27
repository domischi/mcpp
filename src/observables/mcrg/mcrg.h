#ifndef MCPP_MCRG_H_
#define MCPP_MCRG_H_

#include <alps/lattice.h>
#include <cassert>
#include <valarray>
#include <alps/parameter.h>
#include <tuple>
#include "interactions.h"
#include "reductions.h"
#include "../../utilities.h"
#include "../observable.h"
class mcrg : public observable {
public:
    typedef double spin_t;
    typedef mcpp::mcrg::shift_t shift_t; //(dx,dy, component)
    typedef mcpp::mcrg::ReductionTechnique ReductionTechnique;
    typedef mcpp::mcrg::LatticeType LatticeType;

    mcrg(const alps::Parameters& p, const int MCRG_It_, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
    observable(p,gh_,hl_),
    iteration(0),
    max_iterations(MCRG_It_),
    L(p["L"]),
    L_adapted(L),
    N_last(mcpp::init_N(p)),
    lattice_type(mcpp::mcrg::init_lattice_type(p)),
    reduction_type(mcpp::mcrg::init_reduction_technique(p)),
    scale_factor_b(init_b(reduction_type)),
    interactions(mcpp::mcrg::init_interactions(p)),
    entry_point(N_last-1)
    {
    }

    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs){
        ++entry_point%=N_last;//loop over the sites as entry points
        std::vector<std::vector<double>> Se,So; 
        auto reduced_spins=spins;
        //measure the S_alpha in the not reduced lattice
        for(int iteration=0;iteration<=max_iterations;++iteration){
            L_adapted= L/static_cast<int>(std::pow(scale_factor_b,iteration));
            std::vector<double> oute, outo;
            std::tie(oute,outo)=mcpp::mcrg::all_S_alphas(reduced_spins, interactions, lattice_type);
            Se.push_back(oute);
            So.push_back(outo);
            if(!is_last_iteration(iteration))
                reduced_spins=mcpp::mcrg::reduce(reduced_spins, entry_point,L_adapted, reduction_type, is_first_iteration(iteration));
        }
        for(int iteration=0;iteration<=max_iterations;++iteration){
            // Save the correlators <S_alpha>
            std::valarray<double> OUTe(Se[iteration].data(), Se[iteration].size());
            std::valarray<double> OUTo(So[iteration].data(), So[iteration].size());
            std::valarray<double> save_outo=OUTo;//Massive cheat, and I don't exactly know why
            std::valarray<double> save_oute=OUTe;
            obs["MCRGe S_alpha"+ std::to_string(iteration)]<<save_oute;
            obs["MCRGo S_alpha"+ std::to_string(iteration)]<<save_outo;
            //measure <S_alpha n S_beta n>
            std::valarray<double> outouto=mcpp::outer(OUTo,OUTo);
            std::valarray<double> outoute=mcpp::outer(OUTe,OUTe);
            obs["MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outouto;
            obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outoute;
            if(!is_last_iteration(iteration)){//This is not the last instantation of the mcrg
                std::valarray<double> INe(Se[iteration+1].data(), Se[iteration+1].size());
                std::valarray<double> INo(So[iteration+1].data(), So[iteration+1].size());
                //calculate <S_alpha n-1 S_beta n>
                //and save it into obs
                std::valarray<double> outine=mcpp::outer(OUTe, INe);
                std::valarray<double> outino=mcpp::outer(OUTo, INo);
                obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outine;
                obs["MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outino;
            }
        }
        return;
    }
    void init_observables(alps::ObservableSet& obs) const {
        if(is_first_iteration(iteration)) 
            for(int i = 0; i<=max_iterations;++i) {
                obs << alps::RealVectorObservable("MCRGe S_alpha"+ std::to_string(i));
                obs << alps::RealVectorObservable("MCRGe S_alpha"+std::to_string(i) +" S_beta"+std::to_string(i));
                if(!is_last_iteration(i))
                    obs << alps::RealVectorObservable("MCRGe S_alpha"+std::to_string(i) +" S_beta"+std::to_string(i+1));
                obs << alps::RealVectorObservable("MCRGo S_alpha"+ std::to_string(i));
                obs << alps::RealVectorObservable("MCRGo S_alpha"+std::to_string(i) +" S_beta"+std::to_string(i));
                if(!is_last_iteration(i))
                    obs << alps::RealVectorObservable("MCRGo S_alpha"+std::to_string(i) +" S_beta"+std::to_string(i+1));
            }
    }
    void save(alps::ODump &dump) const{
        dump
            << entry_point
            << L_adapted;
    }
    void load(alps::IDump &dump) {
        dump 
            >> entry_point
            >> L_adapted;
    }
private:
    const int iteration;
    const int max_iterations;
    const int L,N_last;
    int L_adapted;
    const LatticeType lattice_type;
    const ReductionTechnique reduction_type;
    //first index: which interaction
    //second index: listing the vectors of the shift
    //third index (shift_t): a pair of ints denoting the shift wrt to the first one (always (0,0))
    const std::vector<std::vector<shift_t>> interactions;
    const int scale_factor_b;
    int entry_point;
    void update_entry_point(){
        entry_point=mcpp::mcrg::new_entry_point(entry_point, L,scale_factor_b);
    }
    bool inline is_last_iteration(const int i) const{
        return max_iterations<=i;
    }
    bool inline is_first_iteration(const int i) const {
        return !i; //if iteration==0 this is true
    }
};
#endif //MCPP_MCRG_H_
