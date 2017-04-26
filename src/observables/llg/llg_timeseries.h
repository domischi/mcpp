#ifndef MCPP_LLG_TIMESERIES_H_
#define MCPP_LLG_TIMESERIES_H_

#include <string>
#include <alps/parameter.h>

class llg_timeseries{
public:
    llg_timeseries(std::string const& name_, const int length) :
    current_index(0),
    values(length),
    name(name){
    }
    void init_observable_set(alps::ObservableSet& obs){
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

protected:
    int current_index;
    std::valarray<double> values;
    std::string name;
};
#endif//MCPP_LLG_TIMESERIES_H_
