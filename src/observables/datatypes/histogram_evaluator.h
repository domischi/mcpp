#ifndef MCPP_HISTOGRAM_EVALUATOR_H_
#define MCPP_HISTOGRAM_EVALUATOR_H_

#include <string>
#include <alps/parameter.h>
class histogram_evaluator {
private:
    const int nbins;
    const std::string name;
    const double min, max;
    const histogram::hist_type type;

    std::valarray<double> lin_x() const{
        std::valarray<double> ret(nbins);
        for(int i=1;i<=nbins-2;++i){
            double xi=(max-min)*i/nbins+min;
            double xip1=(max-min)*(i+1.)/nbins+min;
            double xm=(xi+xip1)/2;
            ret[i]=xm;
        }
        ret[0      ]=min-(max-min)/(nbins/2.);
        ret[nbins-1]=max+(max-min)/(nbins/2.);
        return ret;
    }
    std::valarray<double> log_x() const {
        std::valarray<double> ret(nbins);
        for(int i=1;i<=nbins-2;++i){
            //i==(std::log(xi)-std::log(min))/(std::log(max)-std::log(min))*(nbins.size()-2)+1
            //<=> xi=min*(max/min)^{(k-1)/(n-2)}
            double xi  =min*std::pow(max/min,(i-1.)/(nbins-1.));
            double xip1=min*std::pow(max/min,(i   )/(nbins-1.));
            ret[i]=(xi+xip1)/2.;
        }
        ret[0      ]=min/2;
        ret[nbins-1]=max+(max-min)/(nbins/2.);
        return ret;
    }
    std::valarray<double> slg_x() const{
        std::valarray<double> ret(nbins);
        for(int i=1;i<=nbins/2-2;++i){
            //i==-(std::log(-xi)-std::log(max))/(std::log(max)-std::log(min))*(nnbins/2-2)
            //<=> xi=-max*(min/max)^{i/(nbins/2.-2)}
            double xi  =-max*std::pow((min/max),(i   )/(nbins/2.-2));
            double xip1=-max*std::pow((min/max),(i+1.)/(nbins/2.-2));
            double xm=(xi+xip1)/2.;
            ret[i]=xm;
        }
        for(int i=nbins/2+1;i<=nbins-2;++i){
            //i==(1+(std::log(xi)-std::log(min))/(std::log(max)-std::log(min)))*(nbins.size()/2-2)+3
            //<=> xi=min*(max/min)^{(k-3.)/(n/2.-2)-1}
            double xi  =min*std::pow(max/min,(i-3.)/(nbins/2.-2)-1);
            double xip1=min*std::pow(max/min,(i-2.)/(nbins/2.-2)-1);
            double xm=(xi+xip1)/2.;
            ret[i]=xm;
        }
        ret[0        ]=-(max+(max-min)/((nbins-2)/2/2.));
        ret[nbins-1  ]=  max+(max-min)/((nbins-2)/2/2.) ;
        ret[nbins/2  ]= min/2;
        ret[nbins/2-1]=-min/2;
        return ret;
    }

    std::valarray<double> init_x() const{
        switch(type){
            case histogram::hist_type::linear:
                return lin_x();
            case histogram::hist_type::log:
                return log_x();
            case histogram::hist_type::symmetric_log:
                return slg_x();
        }
    }
public:
    histogram_evaluator(std::string const& name_, const double min_, const double max_, const int n_bins=128, const histogram::hist_type type_ = histogram::hist_type::linear) :
    nbins(n_bins),
    min(min_),
    max(max_),
    name(name_),
    type(type_){}

    virtual void evaluate(alps::ObservableSet& obs) const {
       if(obs.has(name)){
            alps::RealVectorObservable x_vals(name+" X");
            std::valarray<double> x=init_x();
            for (int i=0;i<4;++i) x_vals<<x;
            obs.addObservable(x_vals);
       } else{
           std::cerr<< "x-values for histogram "<< name << " will not be calculated"<<std::endl; 
       } 
    }
};
#endif//MCPP_HISTOGRAM_EVALUATOR_H_
