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
    Br(p.value_or_default("Field Histogram Br" ,1.)),
    V(p.value_or_default("Field Histogram Volume Dot" ,1.)),
    prefactor(1.),
    use_si(p.value_or_default("Field Histogram Use SI", false)),
    auto_detect_range(p.value_or_default("Field Histogram auto detect range",true)),
    maxfield_(mcpp::maxfield(Br,real_z,radius,diameter_split))
    {
        if (use_si){
            if(!p.defined("Field Histogram Br")||!p.defined("Field Histogram Volume Dot")||!p.defined("Field Histogram real z")){
                std::cerr<<"\"M\" or \"Field Histogram real z\" not defined. No Field Histogram X value will be calculated"<<std::endl;
                return;
            }
            prefactor*=Br*V/std::pow(real_z,3)/* mu_0/4pi */;
        } 
        bool log=p.value_or_default("Field Histogram Log Scale", false);
        const int n_bins= p.value_or_default("Field Histogram number of bins",p.value_or_default("Field Histogram n_bins", 128));
        if(auto_detect_range){
            if(log) {
                he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", n_bins,histogram::hist_type::log, prefactor);
                he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , n_bins,histogram::hist_type::symmetric_log, prefactor);
                he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , n_bins,histogram::hist_type::symmetric_log, prefactor);
            }
            else{
                he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", n_bins,histogram::hist_type::linear, prefactor);
                he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , n_bins,histogram::hist_type::linear, prefactor);
                he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , n_bins,histogram::hist_type::linear, prefactor);
            }

        }
        else{
            if(log) {
                const double lower_factor =p.value_or_default("Field Histogram lower histogram truncation" ,4096);
                const double higher_factor=p.value_or_default("Field Histogram higher histogram truncation",(log?64:3));
                he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", maxfield_/lower_factor, higher_factor*maxfield_, n_bins,histogram::hist_type::log, prefactor);
                he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , maxfield_/lower_factor, higher_factor*maxfield_, n_bins,histogram::hist_type::symmetric_log, prefactor);
                he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , maxfield_/lower_factor, higher_factor*maxfield_, n_bins,histogram::hist_type::symmetric_log, prefactor);
            }
            else{
                he_abs=std::make_shared<histogram_evaluator>("Field Histogram Absolute Value", 0.        , maxfield_, n_bins,histogram::hist_type::linear, prefactor);
                he_xy =std::make_shared<histogram_evaluator>("Field Histogram xy"            , -maxfield_, maxfield_, n_bins,histogram::hist_type::linear, prefactor);
                he_z  =std::make_shared<histogram_evaluator>("Field Histogram z"             , -maxfield_, maxfield_, n_bins,histogram::hist_type::linear, prefactor);
            }
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
    const double Br, V;
    double maxfield_;
    double prefactor;
    const bool auto_detect_range;
    const bool use_si;
    std::shared_ptr<histogram_evaluator> he_abs,he_xy,he_z;
};
#endif//MCPP_FIELD_HISTOGRAM_EVALUATOR_H_
