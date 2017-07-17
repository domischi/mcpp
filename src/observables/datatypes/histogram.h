#ifndef MCPP_HISTOGRAM_H_
#define MCPP_HISTOGRAM_H_

#include <string>
#include <alps/parameter.h>
// If linear, log min means mean
// If symmetric_log, min means smallest value above 0 to get its own bin

class histogram{
public:
    enum hist_type {linear, log, symmetric_log};
    histogram(std::string const& name_, const double min_, const double max_, const int n_bins=128, const hist_type type_=hist_type::linear) :
    bins(n_bins),
    min(min_),
    max(max_),
    name(name_),
    type(type_){
    }
    void init_observable(alps::ObservableSet& obs){
        obs<<alps::RealVectorObservable(name);
    }
    virtual void measure(const double value) {
        if(std::isnan(value)){
            std::cerr << "WARNING: Tried to save NaN value, Skipping this entry."<<std::endl;
            return;
        }
        switch(type){
            case linear:
                if(value <=min)
                    bins[0]+=1.;
                else if(value > max)
                    bins[bins.size()-1]+=1.;
                else 
                    bins[int((value-min)/(max-min)*(bins.size()-2))]+=1.;
                break;
            case log:
                if (value<=0. || std::log(value)<std::log(min))
                    bins[0]+=1.;
                else if (std::log(value)>std::log(max))
                    bins[bins.size()-1]+=1.;
                else
                    bins[int((std::log(value)-std::log(min))/(std::log(max)-std::log(min))*(bins.size()-2))+1]+=1.;
                break;
            case symmetric_log:
                if(value<=-max)
                   bins[0]+=1.;
                else if(value>=max)
                   bins[bins.size()-1]+=1.;
                else if(std::abs(value) < min)
                    if (value>0) //positive and small
                        bins[bins.size()/2]+=1.;
                    else //negative and small
                        bins[bins.size()/2-1]+=1.;
                else 
                    if (value>0) //positive and falling in a bin
                        bins[int((1+(std::log(value)-std::log(min))/(std::log(max)-std::log(min)))*(bins.size()/2-2)+3)]+=1.;
                    else
                        bins[-std::floor((std::log(-value)-std::log(max))/(std::log(max)-std::log(min))*(bins.size()/2-2))]+=1.;
                break;
        }
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
    const hist_type type;
};

#include "histogram_evaluator.h"
#endif//MCPP_HISTOGRAM_H_
