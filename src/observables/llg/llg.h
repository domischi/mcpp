#ifndef MCPP_LLG_H_
#define MCPP_LLG_H_

#include <alps/parameter.h>
#include "../observable.h"
#include <alps/lattice.h>
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
    measures_per_measure(p.value_or_default("LLG measures per measure",1)),
    descending_ascending(false),
    rnd(0.,1.),
    rng(static_cast<int>(p["SEED"]))
    {
        std::tie(nl, difference_vectors) = mcpp::get_neighbours(*graph_helper_, p);
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
        }
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealObservable("LLG M");
        obs << alps::RealObservable("LLG Ms");
        obs << alps::RealObservable("LLG E");
    }

    // Intentionally left empty as only const values are inside the class (initialization due to constructor)
    void save(alps::ODump &dump) const{ 
        dump << descending_ascending;
    }
private:
    const int measure_frequency;
    const int measures_per_measure;
    const int thermalization;
    const double dt, alpha, gamma;
    const double T;
    const double D;
    std::mt19937 rng;
    std::normal_distribution<double> rnd;
    typedef mcpp::spin_t spin_t;
    typedef mcpp::vector_type vector_type;
    std::vector<spin_t> spins;
    bool descending_ascending;
    mcpp::neighbour_list_type nl;
    mcpp::difference_vector_list_type difference_vectors;

    void update_spins() {
        if(descending_ascending) for(int i=0  ;i< N;++i) update_one_spin(i);
        else                     for(int i=N-1;i>=0;--i) update_one_spin(i);
        descending_ascending=!descending_ascending;
    }
    void update_one_spin(int i) {
        update_one_spin_stochastic   (i, dt/2);
        update_one_spin_deterministic(i, dt  );
        update_one_spin_stochastic   (i, dt/2);
    }
    void update_one_spin_stochastic(int i, double dt){
        spins[i]+=-gamma*sqrt_mu()*dt*rnd(rng); 
    }
    inline double sqrt_mu() const {
        return std::sqrt(2*alpha*T/(gamma*dt));
    }
    void update_one_spin_deterministic(int i, double dt){
        auto h = get_field_at_site(i);
        spins[i]+=dt*alpha*gamma*(std::sin(spins[i])*h[0]-std::cos(spins[i])*h[1]);
    }
    vector_t get_one_field(vector_type const& r_ij, double const& s) {
        return get_one_field(r_ij, {std::cos(s),std::sin(s)});
    }
    vector_t get_one_field(vector_type const& r_ij, vector_type const& s_j) {
        double norm=mcpp::norm2(r_ij);
        double r_ij_dot_s_j=mcpp::dot(r_ij,s_j);
        vector_type h{0.,0.};
        h[0]=(D/(2.*std::pow(norm,3)))*(s_j[0]-3*r_ij[0]*r_ij_dot_s_j/(norm*norm));
        h[1]=(D/(2.*std::pow(norm,3)))*(s_j[1]-3*r_ij[1]*r_ij_dot_s_j/(norm*norm));
        return h;    
    }
    vector_t get_field_at_site(int i) {
        vector_t h{0.,0.};
        for(int j=0;j<nl[i].size();++j){
            vector_type tmp=get_one_field(difference_vectors[i][j],spins[nl[i][j]]);
            h[0]+=tmp[0];
            h[1]+=tmp[1];
        }
        return h; 
    }
};
#endif //MCPP_LLG_OBSERVABLES_H_
