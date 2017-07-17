#ifndef MCPP_LLG_MUON_DEPOLARIZATION_H_
#define MCPP_LLG_MUON_DEPOLARIZATION_H_

//TODO make this work with more than one muon
//TODO make use of the new function in utilities that does this in a fully fashioned way
class muon_depolarization : public llg_timeseries{
public:
    typedef std::vector<double> vector_t;
    typedef mcpp::spin_t spin_t;
    typedef mcpp::site_iterator site_iterator;
    typedef mcpp::vector_type vector_type;

    muon_depolarization(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_) :
    llg_timeseries("LLG Muon Depolarization", p["LLG measures per measure"]),
    graph_helper_(gh_),
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
    dt(p.value_or_default("dt",1)),
    D(p.value_or_default("D" ,1)),
    muon_spin(3),
    muon_location(3),
    nd(0.,1.),
    urd(cutoff_distance_muon,static_cast<double>(p.value_or_default("a",1.))*static_cast<double>(p["L"])-cutoff_distance_muon),
    rng(static_cast<int>(p["SEED"])), 
    cutoff_distance_muon(p.value_or_default("LLG Muon cutoff",3.)) {
        std::tie(is_deleted, coordinates) = mcpp::get_coordinates(*graph_helper_, p);
    }

    virtual void measure(std::vector<spin_t> const& spins){
        if(current_index==0) 
            initialize_muon();
        llg_timeseries::measure(muon_in_plane_polarized ? muon_spin[0] : muon_spin[2]);
        
    }
    virtual void save(alps::ODump &dump) const{
        llg_timeseries::save(dump);
        dump
            << muon_spin
            << muon_location;
    }
    virtual void load(alps::IDump &dump){
        llg_timeseries::load(dump);
        dump
            >> muon_spin
            >> muon_location;
    }

    void update_muon_spin(const std::vector<spin_t>& s){ 
        vector_t h(get_field_at_muon(s, muon_location));
        muon_spin[0]+=dt*gamma_mu*(muon_spin[1]*h[2]-muon_spin[2]*h[1]);
        muon_spin[1]+=dt*gamma_mu*(muon_spin[2]*h[0]-muon_spin[0]*h[2]);
        muon_spin[2]+=dt*gamma_mu*(muon_spin[0]*h[1]-muon_spin[1]*h[0]);
        double n=mcpp::norm(muon_spin);
        for(auto& x: muon_spin) x/=n;
    }
private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> urd;
    std::normal_distribution<double> nd; 
    const bool muon_in_plane_polarized, fixed_xy, fixed_z;
    const double x_fix, y_fix, z_fix, z_mean, z_std, z_cutoff;
    const double gamma_mu, dt, D;
    const double cutoff_distance_muon;
    vector_t muon_spin;
    vector_t muon_location;
    
    std::shared_ptr<alps::graph_helper<>> graph_helper_;
    std::vector<vector_type> coordinates;
    std::vector<bool> is_deleted;

    inline void initialize_muon() {
        if(muon_in_plane_polarized)
            muon_spin=vector_t{1,0,0};
        else
            muon_spin=vector_t{0,0,1};
        muon_location[0]=get_x();
        muon_location[1]=get_y();
        muon_location[2]=get_z();
    }
    vector_t get_field_at_muon(const std::vector<spin_t>& s, vector_t xyz) const{
        vector_t h={0,0,0};
        for(site_iterator site=graph_helper_->sites().first; site!=graph_helper_->sites().second; ++site){
            vector_t c=coordinates[*site];
            if(std::pow(c[0]-muon_location[0],2)+std::pow(c[1]-muon_location[1],2)+std::pow(muon_location[2],2)<cutoff_distance_muon*cutoff_distance_muon){
                vector_t coor{c[0],c[1],0.};//need a 3d vector here
                vector_t tmp=mcpp::get_one_dipolar_field(coor, muon_location, s[*site], D);
                h[0]+=tmp[0];
                h[1]+=tmp[1];
                h[2]+=tmp[2];
            }
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
};

#endif//MCPP_LLG_MUON_H_
