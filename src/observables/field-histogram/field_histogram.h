#ifndef MCPP_FIELD_HISTOGRAM_H_
#define MCPP_FIELD_HISTOGRAM_H_

#include "../../utilities.h"
#include "../datatypes/histogram.h"
class field_histogram : public observable{
public:
    typedef mcpp::vector_type vector_type;
    
    field_histogram(const alps::Parameters& p, std::shared_ptr<alps::graph_helper<>> gh_, std::shared_ptr<Hamiltonian_List> hl_) :
    observable(p,gh_,hl_),
    graph_helper_(gh_),
    measures_per_measure(p.value_or_default("Field Histogram Measures per Measure",1024)),
    fixed_z(p.value_or_default("Field Histogram fixed z", true)),
    z_mean(p.value_or_default("Field Histogram z", 0.75)),
    z_std (p.value_or_default("Field Histogram z stddev", 0.5)),
    z_cutoff(p.value_or_default("Field Histogram z cutoff", 0.1)),
    cutoff(p.value_or_default("Field Histogram cutoff" ,3.)),
    diameter_split(p.value_or_default("Field Histogram diameter split" ,1)),
    radius(static_cast<double>(p.value_or_default("a",1.))*static_cast<double>(p.value_or_default("Dot Radius" ,0.35))), // one of the samples has this configuration
    M(p.value_or_default("M" ,1.)),
    urd(0.,static_cast<double>(p.value_or_default("a",1.))*static_cast<double>(p["L"])),
    nd(0.,1.),
    rng(static_cast<int>(p["SEED"]))
    {
        std::tie(is_deleted, coordinates) = mcpp::get_coordinates(*graph_helper_, p);
        periodic_translations_=mcpp::get_periodic_translations(*graph_helper_, p);
        
        maxfield_=maxfield();
        bool log=p.value_or_default("Field Histogram Log Scale", false);
        if(log) {
            histogram_abs=std::make_shared<histogram>("Field Histogram Absolute Value", maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::log);
            histogram_xy =std::make_shared<histogram>("Field Histogram xy"            , maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::symmetric_log);
            histogram_z  =std::make_shared<histogram>("Field Histogram z"             , maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::symmetric_log);
        }                                                                               
        else{
            histogram_abs=std::make_shared<histogram>("Field Histogram Absolute Value", 0.        , maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
            histogram_xy =std::make_shared<histogram>("Field Histogram xy"            , -maxfield_, maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
            histogram_z  =std::make_shared<histogram>("Field Histogram z"             , -maxfield_, maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
        }
    }

    virtual void measure(std::vector<spin_t> const& spins, alps::ObservableSet& obs){
        add_constant(obs["Field Histogram MaxField"], maxfield_);
        for(int i=0;i<measures_per_measure;++i){
            vector_type pos={get_x(),get_y(),get_z()};
            vector_type field=get_field_at_position(spins, pos);
            histogram_abs->measure(mcpp::norm(field));
            histogram_xy ->measure(field[0]);
            histogram_xy ->measure(field[1]);
            histogram_z  ->measure(field[2]);
        }
        histogram_abs->write(obs);
        histogram_xy ->write(obs);
        histogram_z  ->write(obs);
    }
    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::SimpleRealObservable("Field Histogram MaxField");
        histogram_abs->init_observable(obs);
        histogram_xy ->init_observable(obs);
        histogram_z  ->init_observable(obs);
    }
    virtual void save(alps::ODump &dump) const{
        histogram_abs->save(dump);
        histogram_xy ->save(dump);
        histogram_z  ->save(dump);
    }
    virtual void load(alps::IDump &dump){
        histogram_abs->load(dump);
        histogram_xy ->load(dump);
        histogram_z  ->load(dump);
        std::cerr << "WARNING: reload a simulation of Field Histogram, however the internal RNG is not backed up\nThe results are therefore not bitwise reproducible"<<std::endl; //TODO
    }

private:
    std::shared_ptr<alps::graph_helper<>> graph_helper_;
    std::vector<vector_type> coordinates;
    std::vector<bool> is_deleted;
    std::vector<vector_type> periodic_translations_;
   
    const int measures_per_measure;
     
    const double cutoff;
    const int diameter_split;
    const double radius;
    const double M;

    const bool fixed_z;
    const double z_mean, z_std, z_cutoff;
    std::mt19937 rng;
    std::normal_distribution<double> nd; 
    std::uniform_real_distribution<double> urd;
    double maxfield_;
    std::shared_ptr<histogram> histogram_abs, histogram_xy, histogram_z;

    vector_type get_field_at_position(const std::vector<spin_t>& s, vector_type xyz) const{
        return mcpp::get_field_at_position(graph_helper_, s, xyz, coordinates, is_deleted, periodic_translations_, cutoff, diameter_split, radius, M);
    }
    inline double maxfield() const {
        return 3*mcpp::norm(mcpp::get_field_one_dot(vector_type {0.,0.,0.} , vector_type{0.,0.,z_mean}, 0., M, diameter_split,radius)); // The 3 in front seems more to be of a phenomenological structure taking into account the neighbours withour really calculating them
    }

    inline double get_x(){
        return urd(rng);
    }
    inline double get_y(){
        return urd(rng);
    }
    inline double get_z(){
        if(fixed_z)
            return z_mean;
        else{
            double z=-1;
            while(z<=z_cutoff){//ignore the very last part close to the sample due to numerical instabilities
                z=z_mean+(nd(rng)*z_std);
            }
            return z;
        }
    }
};

#endif//MCPP_FIELD_HISTOGRAM_H_
