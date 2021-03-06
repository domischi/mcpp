#ifndef MCPP_SPIN_AUTOCORRELATION_H_
#define MCPP_SPIN_AUTOCORRELATION_H_

#include <alps/lattice.h>
#include <alps/parameter.h>
#include "../../utilities.h"
#include "../observable.h"
class spin_autocorrelation : public observable {
public:
    typedef double spin_t;
    typedef std::vector<spin_t> configuration_t;

    spin_autocorrelation(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
    observable(p,gh_,hl_),
    filled(false),
    actual_index(0),
    autocorr_analysis_depth(p["Spin autocorrelation analysis length"])
    {
        old_configurations.resize(autocorr_analysis_depth); //need to store the last n configurations to calculate the autocorrelation
    }
                
    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) {
        old_configurations[actual_index]=spins;
        if(filled){
            std::valarray<double> acs(autocorr_analysis_depth);
            for(int i =1;i<autocorr_analysis_depth;++i){
                configuration_t old_configuration=old_configurations[(actual_index-i+autocorr_analysis_depth)%autocorr_analysis_depth];
                acs[i]=mcpp::spin_config_overlap(old_configuration,spins);
            }
            acs[0]=1; //the spins are perfectly self correlated 
            obs["Spin Autocorrelation"]<<acs; //TODO do this in terms of the mcpp utility function
        }
        else{
            filled=(actual_index==autocorr_analysis_depth-1);
        }
        ++actual_index%=autocorr_analysis_depth;
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealVectorObservable("Spin Autocorrelation");
    }

    void save(alps::ODump &dump) const{
        dump
            << actual_index 
            << filled 
            << old_configurations;
    }

    void load(alps::IDump &dump){
        dump
            >> actual_index 
            >> filled 
            >> old_configurations;
    }
private:
    //const int L,N;
    const int autocorr_analysis_depth;
    int actual_index;
    bool filled;
    //As I dont want to enforce boost (not really fast, some issues over versions) I implemented my own version of a circular buffer typish thing
    std::vector<configuration_t> old_configurations;

};
#endif //MCPP_SPIN_AUTOCORRELATION_H_
