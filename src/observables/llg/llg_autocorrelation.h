#ifndef MCPP_LLG_AUTOCORRELATION_H_
#define MCPP_LLG_AUTOCORRELATION_H_

#include "llg_timeseries.h"

class llg_autocorrelation : public llg_timeseries{
public:
    typedef mcpp::spin_t spin_t;
    
    llg_autocorrelation(const int length_) : 
    llg_timeseries("LLG Autocorrelation", length_),
    first_spins(0){
    }
    virtual void measure(std::vector<spin_t> const& spins){
        if(current_index==0) {
            first_spins=spins;
            llg_timeseries::measure(1.); //perfect autocorrelation time at t=0
        } 
        else 
            llg_timeseries::measure(mcpp::spin_config_overlap(first_spins,spins));
    }
    virtual void save(alps::ODump &dump) const{
        llg_timeseries::save(dump);
        dump
            << first_spins;
    }
    virtual void load(alps::IDump &dump){
        llg_timeseries::load(dump);
        dump
            >> first_spins;
    }
private:
    std::vector<spin_t> first_spins;
};
#endif//MCPP_LLG_AUTOCORRELATION_H_
