#ifndef MCPP_HISTOGRAM_H_
#define MCPP_HISTOGRAM_H_

#include <string>
#include <alps/parameter.h>
// If linear, log min means mean
// If symmetric_log, min means smallest value above 0 to get its own bin

class histogram{
public:
    enum hist_type {linear, log, symmetric_log};
    //fixed range histogram
    histogram(std::string const& name_, const double min_, const double max_, const int n_bins=128, const hist_type type_=hist_type::linear) :
    bins(n_bins),
    min(min_),
    max(max_),
    name(name_),
    auto_detect_range(false),
    n_first_values(0),
    writeable(true),
    type(type_){
    }
    // auto-detect-range histogram
    histogram(std::string const& name_, const int n_bins=128, const int n_first_values_=1024, const hist_type type_=hist_type::linear) :
    bins(n_bins),
    min(0),
    max(0),
    name(name_),
    auto_detect_range(true),
    n_first_values(n_first_values_),
    writeable(false),
    type(type_){
    }
    void init_observable(alps::ObservableSet& obs){
        obs<<alps::RealVectorObservable(name);
        obs<<alps::SimpleRealObservable(name+" Min");
        obs<<alps::SimpleRealObservable(name+" Max");
    }
    virtual void measure(const double value) {
        if(std::isnan(value)){
            std::cerr << "WARNING: Tried to save NaN value, Skipping this entry."<<std::endl;
            return;
        }
        if(auto_detect_range && !writeable){ //fill first_values
            first_values.push_back(value);
            if(first_values.size()>=n_first_values){
                if(type==hist_type::symmetric_log){// special care has to be taken for log scaled ones
                    min=std::abs(*std::min_element(std::begin(first_values),std::end(first_values), [](double a, double b){return std::abs(a)<std::abs(b);}));
                }
                else{
                    min=*std::min_element(std::begin(first_values),std::end(first_values));
                }
                max=*std::max_element(std::begin(first_values),std::end(first_values));
                for(auto & x : first_values) fill_bins(x);
                first_values=std::vector<double>();
                writeable=true;
            }
        }
        else{
            fill_bins(value);
        }
    }
    virtual void write(alps::ObservableSet& obs){
        if(writeable){
            const double normalization=bins.sum();
            if(normalization <1.){
                std::cerr<< "Histogram variable "+name+" wants to save before being measured. Not saving anything..."<<std::endl;
                return;
            }
            bins/=normalization;
            obs[name]<<bins;
            add_constant(obs[name+ " Min"], min);
            add_constant(obs[name+ " Max"], max);
            bins=std::valarray<double>(bins.size());
        }
    }
    virtual void save(alps::ODump &dump) const{
        dump
            << min
            << max
            << writeable
            << first_values
            << bins;
    }
    virtual void load(alps::IDump &dump){
        dump
            >> min
            >> max
            >> writeable
            >> first_values
            >> bins;
    }
protected:
    std::valarray<double> bins;
    const std::string name;
    double min, max;
    const hist_type type;
    const bool auto_detect_range;
    const int n_first_values;
    bool writeable;
    std::vector<double> first_values;

    void fill_bins(const double value){
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
};

#include "histogram_evaluator.h"
#endif//MCPP_HISTOGRAM_H_
