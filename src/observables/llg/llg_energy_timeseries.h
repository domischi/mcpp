#ifndef MCPP_LLG_ENERGY_TIMESERIES_H_
#define MCPP_LLG_ENERGY_TIMESERIES_H_

#include "llg_timeseries.h"

class llg_energy_timeseries : public llg_timeseries{
public:
    typedef mcpp::spin_t spin_t;
    
    llg_energy_timeseries(const int length_, std::shared_ptr<Hamiltonian_List> hl_ ) : 
    llg_timeseries("LLG Energy Timeseries", length_),
    hamiltonian_list_(hl_)
    {
    }
    virtual void measure(std::vector<spin_t> const& spins){
        llg_timeseries::measure(hamiltonian_list_->Energy(spins));
    }
private:
    std::shared_ptr<Hamiltonian_List> hamiltonian_list_;
};
#endif//MCPP_LLG_ENERGY_TIMESERIES_H_
