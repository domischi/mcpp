#ifndef MCPP_LLG_H_
#define MCPP_LLG_H_

#include <alps/parameter.h>
#include <alps/lattice.h>
#include "../observable.h"
#include "llg_energy_timeseries.h"
#include "llg_autocorrelation.h"
#include "muon_depolarization.h"
class llg : public observable{
public:
    typedef std::vector<double> vector_t;

    llg(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
    observable(p,gh_,hl_),
    dt(p.value_or_default("dt",1)),
    T(p.defined("T") ? static_cast<double>(p["T"]) : 1./static_cast<double>(p["beta"])),
    alpha(p.value_or_default("alpha",1)),
    gamma(p.value_or_default("gamma",1)),
    D(p.value_or_default("D",1)),
    thermalization(p.value_or_default("LLG thermalization",0)),
    measure_frequency(p.value_or_default("LLG measure frequency",1)),
    measures_per_measure(p["LLG measures per measure"]),
    descending_ascending(false),
    nd(0.,1.),
    rng(static_cast<int>(p["SEED"]))
    {
        std::tie(nl, difference_vectors) = mcpp::get_neighbours(*graph_helper_, p);
        std::tie(is_deleted, coordinates) = mcpp::get_coordinates(*graph_helper_, p);
        if(p.value_or_default("LLG Measure Energy Timeseries",false))
            energy_timeseries=std::make_shared<llg_energy_timeseries>(measures_per_measure, hl_);
        if(p.value_or_default("LLG Measure Autocorrelation",false))
            autocorrelation=std::make_shared<llg_autocorrelation>(measures_per_measure);
        if(p.value_or_default("LLG Measure Muon",false))
            muon=std::make_shared<muon_depolarization>(p,gh_);
    }

    void measure(const std::vector<spin_t>& s, alps::ObservableSet& obs) {
        spins=s;
        for(int j = 0; j< thermalization; ++j) {
            update_spins();
        }
        for(int i = 0; i < measures_per_measure; ++i) {
            for(int j = 0; j< measure_frequency; ++j) {
                update_spins();
            }
            double M, Ms, E;
            std::tie(M,Ms)=mcpp::magnetization_and_staggered_magnetization(spins,L);
            E=hamiltonian_list_->Energy(spins)/graph_helper_->num_sites();
            obs["LLG M" ] << M; 
            obs["LLG Ms"] << Ms;
            obs["LLG E" ] << E;
            if(energy_timeseries)
                energy_timeseries->measure(spins);
            if(muon)
                muon->measure(spins);
            if(autocorrelation) 
                autocorrelation->measure(spins);
        }
        if(energy_timeseries)
            energy_timeseries->write(obs);
        if(muon)
            muon->write(obs);
        if(autocorrelation)
            autocorrelation->write(obs);
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealObservable("LLG M");
        obs << alps::RealObservable("LLG Ms");
        obs << alps::RealObservable("LLG E");
        if(energy_timeseries) 
            energy_timeseries->init_observable(obs); 
        if(muon) 
            muon->init_observable(obs); 
        if(autocorrelation)
            autocorrelation->init_observable(obs);
    }

    void save(alps::ODump &dump) const{
        dump
            << descending_ascending
            << spins;
        if(autocorrelation)
            autocorrelation->save(dump);
        if(muon)
            muon->save(dump);
        if(energy_timeseries)
            energy_timeseries->save(dump);
    }
    void load(alps::IDump &dump){
        dump
            >> descending_ascending
            >> spins;
        std::cerr << "WARNING: reload a simulation of LLG, however the internal RNG is not backed up\nThe results are therefore not bitwise reproducible"<<std::endl; //TODO
        if(autocorrelation)
            autocorrelation->load(dump);
        if(muon)
            muon->load(dump);
        if(energy_timeseries)
            energy_timeseries->load(dump);
    }
private:
    const int measure_frequency;
    const int measures_per_measure;
    const int thermalization;
    const double dt, alpha, gamma;
    const double T;
    const double D;
    std::mt19937 rng;
    std::normal_distribution<double> nd;
    typedef mcpp::spin_t spin_t;
    typedef mcpp::site_iterator site_iterator;
    typedef mcpp::vector_type vector_type;
    std::vector<spin_t> spins;
    bool descending_ascending;
    std::vector<bool> is_deleted;//TODO this is not yet implemented as a use of
    std::vector<vector_type> coordinates;
    mcpp::neighbour_list_type nl;
    mcpp::difference_vector_list_type difference_vectors;

    std::shared_ptr<llg_energy_timeseries> energy_timeseries;
    std::shared_ptr<llg_autocorrelation> autocorrelation;
    std::shared_ptr<muon_depolarization> muon;
    void update_spins() {
        if(descending_ascending) for(int i=0  ;i< N;++i) update_one_spin(i);
        else                     for(int i=N-1;i>=0;--i) update_one_spin(i);
        if(muon)
            muon->update_muon_spin(spins);
        descending_ascending=!descending_ascending;
    }
    void update_one_spin(int i) {
        update_one_spin_stochastic   (i, dt/2);
        update_one_spin_deterministic(i, dt  );
        update_one_spin_stochastic   (i, dt/2);
    }
    void update_one_spin_stochastic(int i, double dt){
        spins[i]+=-gamma*sqrt_mu()*dt*nd(rng);
    }
    inline double sqrt_mu() const {
        return std::sqrt(2*alpha*T/(gamma*dt));
    }
    void update_one_spin_deterministic(int i, double dt){
        auto h = get_field_at_site(i);
        spins[i]+=dt*alpha*gamma*(std::sin(spins[i])*h[0]-std::cos(spins[i])*h[1]);
    }
    vector_t get_field_at_site(int i) const {
        vector_t h{0.,0.};
        for(int j=0;j<nl[i].size();++j){
            vector_type tmp=mcpp::get_one_dipolar_field(difference_vectors[i][j],spins[nl[i][j]], D);
            h[0]+=tmp[0];
            h[1]+=tmp[1];
        }
        return h; 
    }
};
#endif //MCPP_LLG_H_
