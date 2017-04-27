#ifndef MCPP_LLG_TIMESERIES_H_
#define MCPP_LLG_TIMESERIES_H_

#include <string>
#include <alps/parameter.h>

class llg_timeseries{
public:
    llg_timeseries(std::string const& name_, const int length) :
    current_index(0),
    values(length),
    name(name_){
    }
    void init_observable(alps::ObservableSet& obs){
        obs<<alps::RealVectorObservable(name);
    }
    virtual void measure(double value) {
        if(current_index>values.size()){
            std::cerr<<"LLG timeseries encountered an error.\nYou try to store more data into the observable "+name+" then there are allowed. Aborting..."<<std::endl;
            std::exit(43);
        }
        values[current_index]=value;
        ++current_index;
    }
    virtual void write(alps::ObservableSet& obs){
        if(current_index!=values.size()) {
            std::cerr<<"LLG timeseries encountered an error.\nYou try to write before filling up in observable "+name+". Aborting..."<<std::endl;
            std::exit(45);
        }
        obs[name]<<values;
        values=std::valarray<double>(values.size());
        current_index=0;
    }
    virtual void save(alps::ODump &dump) const{
        dump
            << current_index
            << values;
    }
    virtual void load(alps::IDump &dump){
        dump
            >> current_index
            >> values;
    }

protected:
    int current_index;
    std::valarray<double> values;
    const std::string name;
};
#endif//MCPP_LLG_TIMESERIES_H_
