#ifndef MCPP_HISTOGRAM_H_
#define MCPP_HISTOGRAM_H_

#include <string>
#include <alps/parameter.h>

class histogram{
public:
    histogram(std::string const& name_, const double min_, const double max_, const int n_bins=128) :
    bins(n_bins),
    min(min_),
    max(max_),
    name(name_){
    }
    void init_observable(alps::ObservableSet& obs){
        obs<<alps::RealVectorObservable(name);
    }
    virtual void measure(const double value) {
        if(std::isnan(value)){
            std::cerr << "WARNING: Tried to save NaN value, Skipping this entry."<<std::endl;
            return;
        }
        if(value <=min)
            bins[0]+=1.;
        else if(value > max)
            bins[bins.size()-1]+=1.;
        else 
            bins[int((value-min)/(max-min)*(bins.size()-2))]+=1.;
    }
    virtual void write(alps::ObservableSet& obs){
        const double normalization=bins.sum();
        if(normalization <1.){
            std::cerr<< "Histogram variable "+name+" wants to save before being measured. Not saving anything..."<<std::endl;
            return;
        }
        bins/=normalization;
        obs[name]<<bins;
        bins=std::valarray<double>(bins.size());
    }
    virtual void save(alps::ODump &dump) const{
        dump
            << bins;
    }
    virtual void load(alps::IDump &dump){
        dump
            >> bins;
    }
protected:
    std::valarray<double> bins;
    const std::string name;
    const double min, max;
};
#endif//MCPP_HISTOGRAM_H_
