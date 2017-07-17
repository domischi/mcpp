#ifndef MCPP_FIELD_HISTOGRAM_EVALUATOR_H_
#define MCPP_FIELD_HISTOGRAM_EVALUATOR_H_

class field_histogram_evaluator{
public: 
    field_histogram_evaluator(alps::Parameters const& p)
    :
    z(p.value_or_default("Field Histogram z", 0.75)),
    real_z(p.value_or_default("Field Histogram real z", 0.75)),
    diameter_split(p.value_or_default("Field Histogram diameter split" ,1)),
    radius(static_cast<double>(p.value_or_default("a",1.))*static_cast<double>(p.value_or_default("Dot Radius" ,0.35))),
    M(p.value_or_default("M" ,1.)),
    use_si(p.value_or_default("Field Histogram Use SI", false)),
    maxfield_(mcpp::maxfield(M,real_z,radius,diameter_split))
    {
        if (use_si){
            if(!p.defined("M")||!p.defined("Field Histogram real z")){
                std::cerr<<"\"M\" or \"Field Histogram real z\" not defined. No Field Histogram X value will be calculated"<<std::endl;
                return;
            }
            maxfield_*=1e-7/* mu_0/4pi */;
        } 
    bool log=p.value_or_default("Field Histogram Log Scale", false);
    if(log) {
        he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::log);
        he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::symmetric_log);
        he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , maxfield_/1024, 256*maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::symmetric_log);
    }
    else{
        he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", 0.        , maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
        he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , -maxfield_, maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
        he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , -maxfield_, maxfield_, p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128)),histogram::hist_type::linear);
    }
    }
    virtual void evaluate(alps::ObservableSet& obs) const {
        he_abs->evaluate(obs);
        he_xy ->evaluate(obs);
        he_z  ->evaluate(obs);
    }
private:
    const double z; 
    const double real_z; 
    const int diameter_split;
    const double radius;
    const double M;
    double maxfield_;
    const bool use_si;
    std::shared_ptr<histogram_evaluator> he_abs,he_xy,he_z;
};
#endif//MCPP_FIELD_HISTOGRAM_EVALUATOR_H_
