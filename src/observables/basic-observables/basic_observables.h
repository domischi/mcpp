#ifndef MCPP_BASIC_OBSERVABLES_H_
#define MCPP_BASIC_OBSERVABLES_H_

#include <alps/parameter.h>

class basic_observables {
public:
    typedef double spin_t;

    basic_observables(const alps::Parameters& p) :
    L(p["L"]),
    N(L*L)
    {
    }
                
    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) {
        double M, Ms;
        std::tie(M,Ms)=magnetization_and_staggered_magnetization(spins);
        obs["M"            ] << M; 
        obs["M^2"          ] << M*M;
        obs["M^4"          ] << M*M*M*M;
        obs["M staggered"  ] << Ms;
        obs["M staggered^2"] << Ms*Ms;
        obs["M staggered^4"] << Ms*Ms*Ms*Ms;
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealObservable("M");
        obs << alps::RealObservable("M^2");
        obs << alps::RealObservable("M^4");
        obs << alps::RealObservable("M staggered");
        obs << alps::RealObservable("M staggered^2");
        obs << alps::RealObservable("M staggered^4");
    }

    // Intentionally left empty as only const values are inside the class (initialization due to constructor)
    void save(alps::ODump &dump) const{ }
private:
    const int L,N;

    std::pair<double,double> magnetization_and_staggered_magnetization(const std::vector<spin_t>& spins) const{
        double mx =0.;
        double my =0.;
        double mxs=0.;
        double mys=0.;
        //at least theoretically vectorizable, however should not be the bottleneck
        for(int i=0;i<N;++i) {
            spin_t sp=spins[i];
            double c=std::cos(sp);
            double s=std::sin(sp);
            mx+=c;
            my+=s;
            int prefactor_x=1;
            int prefactor_y=1;
            mxs+= ((i%L)%2 ? -1 : 1) * c; //for even y sites -1
            mys+= ((i/L)%2 ? -1 : 1) * s; //for even x sites -1
        }
        return std::make_pair(std::sqrt(mx*mx+my*my)/N,std::sqrt(mxs*mxs+mys*mys)/N);
    }
};
#endif //MCPP_BASIC_OBSERVABLES_H_
