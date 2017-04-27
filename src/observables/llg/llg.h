#ifndef MCPP_LLG_H_
#define MCPP_LLG_H_

#include <alps/parameter.h>
#include <alps/lattice.h>
#include "../observable.h"
#include "llg_autocorrelation.h"

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
    measure_muon(p.value_or_default("LLG Measure Muon",false)),
    muon_in_plane_polarized(p.value_or_default("LLG Muon in plane",false)),
    fixed_xy(p.defined("LLG Muon X") || p.defined("LLG Muon Y")),
    fixed_z( p.defined("LLG Muon Z")),
    x_fix(p.value_or_default("LLG Muon X",p.value_or_default("LLG Muon Y",-1))),
    y_fix(p.value_or_default("LLG Muon Y",p.value_or_default("LLG Muon X",-1))),
    z_fix(p.value_or_default("LLG Muon Z",-1)),
    z_mean(p.value_or_default("LLG Muon Z mean",1)),
    z_std (p.value_or_default("LLG Muon Z std" ,0.5)),
    z_cutoff(p.value_or_default("LLG Muon Z cutoff" ,0.05)),
    gamma_mu(p.value_or_default("LLG Muon gamma" ,1)),
    muon_spin(3),
    muon_location(3),
    descending_ascending(false),
    nd(0.,1.),
    urd(cutoff_distance_muon,static_cast<double>(p.value_or_default("a",1.))*static_cast<double>(p["L"])-cutoff_distance_muon),
    rng(static_cast<int>(p["SEED"])),
    depolarization(measures_per_measure),
    debug(false),
    cutoff_distance_muon(p.value_or_default("LLG Muon cutoff",3.))
    {
        std::tie(nl, difference_vectors) = mcpp::get_neighbours(*graph_helper_, p);
        std::tie(is_deleted, coordinates) = mcpp::get_coordinates(*graph_helper_, p);
        if(p.value_or_default("LLG Measure Autocorrelation",false))
            autocorrelation=std::make_shared<llg_autocorrelation>(measures_per_measure);
    }
                
    void measure(const std::vector<spin_t>& s, alps::ObservableSet& obs) {
        spins=s;
        if(measure_muon)
            initialize_muon();
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
            if(measure_muon) {
                if(muon_in_plane_polarized) // in plane
                    depolarization[i]=muon_spin[0];
                else // out of plane
                    depolarization[i]=muon_spin[2];
            }
            if(std::isnan(muon_spin[0])&&!debug){
                debug=true;
                std::cout 
                    << std::setw(15)<< muon_spin[0]
                    << std::setw(15)<< muon_spin[1]
                    << std::setw(15)<< muon_spin[2]
                    << std::setw(4)<< i
                    << std::setw(15) << sqrt_mu()
                    << std::setw(15) << muon_location[2]
                    << std::endl;
            }
            if(autocorrelation) {
                autocorrelation->measure(spins);
            }
        }
        if(measure_muon){
            obs["LLG Depolarization"]<< depolarization;
        }
        if(autocorrelation)
            autocorrelation->write(obs);
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealObservable("LLG M");
        obs << alps::RealObservable("LLG Ms");
        obs << alps::RealObservable("LLG E");
        if(measure_muon) {
            obs << alps::RealVectorObservable("LLG Depolarization");
        }
        if(autocorrelation){
            autocorrelation->init_observable(obs);
        }
    }

    void save(alps::ODump &dump) const{
        dump
            << descending_ascending
            << muon_spin
            << muon_location
            << depolarization;
        if(autocorrelation)
            autocorrelation->save(dump);
    }
    void load(alps::IDump &dump){
        dump
            >> descending_ascending
            >> muon_spin
            >> muon_location
            >> depolarization;
        std::cerr << "WARNING: reload a simulation of LLG, however the internal RNG is not backed up\nThe results are therefore not bitwise reproducible"<<std::endl; //TODO
        if(autocorrelation)
            autocorrelation->load(dump);
    }
private:
    const int measure_frequency;
    const int measures_per_measure;
    const int thermalization;
    const double dt, alpha, gamma;
    const double T;
    const double D;
    const double cutoff_distance_muon;
    std::mt19937 rng;
    std::normal_distribution<double> nd;
    std::uniform_real_distribution<double> urd;
    typedef mcpp::spin_t spin_t;
    typedef mcpp::site_iterator site_iterator;
    typedef mcpp::vector_type vector_type;
    std::vector<spin_t> spins;
    bool descending_ascending;
    std::vector<bool> is_deleted;
    std::vector<vector_type> coordinates;
    mcpp::neighbour_list_type nl;
    mcpp::difference_vector_list_type difference_vectors;
    bool debug;
    const bool measure_muon;
    const double gamma_mu;
    const bool fixed_xy, fixed_z, muon_in_plane_polarized;
    const double x_fix, y_fix, z_fix, z_mean, z_std, z_cutoff;

    vector_t muon_spin, muon_location;
    std::valarray<double> depolarization;
    
    std::shared_ptr<llg_autocorrelation> autocorrelation;

    void update_spins() {
        if(descending_ascending) for(int i=0  ;i< N;++i) update_one_spin(i);
        else                     for(int i=N-1;i>=0;--i) update_one_spin(i);
        if(measure_muon)
            update_muon_spin();
        descending_ascending=!descending_ascending;
    }
    void update_muon_spin(){
        vector_t h(get_field_at_muon(muon_location));
        muon_spin[0]+=dt*gamma_mu*(muon_spin[1]*h[2]-muon_spin[2]*h[1]);
        muon_spin[1]+=dt*gamma_mu*(muon_spin[2]*h[0]-muon_spin[0]*h[2]);
        muon_spin[2]+=dt*gamma_mu*(muon_spin[0]*h[1]-muon_spin[1]*h[0]);
        double n=mcpp::norm(muon_spin);
        for(auto& x: muon_spin) x/=n;
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
    vector_t get_one_field(vector_type const& r_ij, double const& s) const {
        return get_one_field(r_ij, {std::cos(s),std::sin(s)});
    }
    vector_t get_one_field(vector_type const& r_ij, vector_type const& s_j) const {
        double norm=mcpp::norm(r_ij);
        double r_ij_dot_s_j=mcpp::dot(r_ij,s_j);
        vector_type h(r_ij.size());
        for(int i = 0; i<h.size();++i){ //this version works in 2D as well as in 3D
            h[i]=(D/(2.*std::pow(norm,3)))*(s_j[i]-3*r_ij[i]*r_ij_dot_s_j/(norm*norm));
        }
        return h;    
    }
    vector_t get_field_at_site(int i) const {
        vector_t h{0.,0.};
        for(int j=0;j<nl[i].size();++j){
            vector_type tmp=get_one_field(difference_vectors[i][j],spins[nl[i][j]]);
            h[0]+=tmp[0];
            h[1]+=tmp[1];
        }
        return h; 
    }
    inline double get_x(){
        if(fixed_xy)
            return x_fix;
        else{
            return urd(rng);
        }
    }
    inline double get_y(){
        if(fixed_xy)
            return y_fix;
        else{
            return urd(rng);
        }
    }
    inline double get_z(){
        if(fixed_z)
            return z_fix;
        else{
            double z=-1;
            while(z<=z_cutoff){//ignore the very last part close to the sample due to numerical instabilities
                z=z_mean+(nd(rng)*z_std);
            }
            return z;
        }
    }
    inline void initialize_muon() {
        debug=false;
        if(muon_in_plane_polarized)
            muon_spin=vector_t{1,0,0};
        else
            muon_spin=vector_t{0,0,1};
        muon_location[0]=get_x();
        muon_location[1]=get_y();
        muon_location[1]=get_z();
    }
    vector_t get_field_at_muon(vector_t xyz) const{
        vector_t h={0,0,0};
        for(site_iterator site=graph_helper_->sites().first; site!=graph_helper_->sites().second; ++site){
            vector_t c=coordinates[*site];
            if(std::pow(c[0]-muon_location[0],2)+std::pow(c[1]-muon_location[1],2)<cutoff_distance_muon*cutoff_distance_muon){
                vector_t tmp=get_one_field(muon_location, spins[*site]);
                h[0]+=tmp[0];
                h[1]+=tmp[1];
                h[2]+=tmp[2];
            }
        }
        return h;
    }
    vector_t get_field_at_muon(const double x, const double y, const double z) const{
        return get_field_at_muon(vector_t{x,y,z});
    }
};
#endif //MCPP_LLG_H_
