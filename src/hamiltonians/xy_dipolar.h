#ifndef MCPP_XY_HAMILTONIAN_DIPOLE_H
#define MCPP_XY_HAMILTONIAN_DIPOLE_H

#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "xy_hamiltonian.h"
#include "../utilities.h"
template<bool DISORDERED>
class XY_Dipole : public XY_Hamiltonian<XY_Dipole<DISORDERED>>{
public:
    XY_Dipole (alps::Parameters const& params) : 
    XY_Hamiltonian<XY_Dipole<DISORDERED>>(params),
    D(params.value_or_default("D",1)),
    L(params["L"]),
    N(mcpp::init_N(params)),
    gh(params),
    print_debug_information(static_cast<bool>(params.value_or_default("debug dipolar",false))){
        if(DISORDERED){
            dilution_rate=params.value_or_default("Dilution Rate", 0.);
            position_std_dev=static_cast<double>(params.value_or_default("Position Disorder", 0.))*static_cast<double>(params.value_or_default("a",1.));
        }
        else{
            dilution_rate=0;
            position_std_dev=0;
        }
        cutoff_distance=static_cast<double>(params.value_or_default("cutoff_distance",3));
        if (cutoff_distance>=L/2.){
            std::cerr<<"Too small lattice for such a cutoff, there is an ambiguity and you should check what you are doing... Aborting..."<<std::endl;
            std::exit(3);
        }
        cutoff_distance*=std::max(static_cast<double>(params.value_or_default("a",1.)),static_cast<double>(params.value_or_default("b",1.)));
        Init_Lookup_Tables(params);
    } 
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i) const {
        double e=0.;
        for(int j : neighbour_list[i]){
            e-=D*inv_distance_cubed(i,j)*(1.5*std::cos(spins[i]+spins[j]-2*angle_w_x(i,j))+0.5*std::cos(spins[i]-spins[j]));
        }
        return e;
    }

private:
    int L,N;
    const bool print_debug_information;
    double D;
    double dilution_rate;
    double position_std_dev;
    double cutoff_distance;
    alps::graph_helper<> gh;
    //Lookup tables 
    std::vector<double> dist_3;
    std::vector<double> phi;
    std::vector<std::vector<int>> neighbour_list; //saves which neighbours are relevant (as not only nearest neighbours count in dipole )
                        
    typedef typename mcpp::vector_type vector_type;
    typedef typename mcpp::basis_vector_iterator basis_vector_iterator;
    typedef typename mcpp::site_iterator site_iterator;

    //This is a templated if, should be optimized away as it is constant in the class
    inline int reduced_index(int i, int j) const {
        if(DISORDERED)
            return i+N*j;
        else {
            int xi,yi,xj,yj,x,y;
            xi=i/L;
            yi=i%L;
            xj=j/L;
            yj=j%L;
            x=(xi-xj+L)%L;
            y=(yi-yj+L)%L;
            return L*x+y;
        }
    }
    inline int reduced_index(std::pair<int,int> p) const {
        return reduced_index(p.first,p.second);
    }
    void Init_Lookup_Tables(alps::Parameters const& p){
        const int DISORDER_SEED=1;
        if(DISORDERED) {
            dist_3=std::vector<double>(N*N);
            phi=std::vector<double>(N*N);
        }
        else {
            dist_3=std::vector<double>(N);
            phi=std::vector<double>(N);
        }
        std::vector<std::vector<vector_type>> difference_vector_list;
        std::tie(neighbour_list, difference_vector_list) = mcpp::get_neighbours(gh,p);
        for(int i = 0 ; i < neighbour_list.size();++i) {
            for(int index=0;index<neighbour_list[i].size();++index){
                const int j = neighbour_list[i][index];
                const int ri=reduced_index(i,j);
                const vector_type d=difference_vector_list[i][index];
                const double dist=mcpp::norm2(d);
                dist_3[ri]=std::pow(dist,-3);
                phi[ri]=std::atan2(-d[1],-d[0]);
            }
        }
    }
    inline double inv_distance_cubed(int i,int j) const {
        return dist_3[reduced_index(i,j)];
    }
    inline double inv_distance_cubed(std::pair<int,int> pair_) const {
        return inv_distance_cubed(pair_.first, pair_.second);
    }
    inline double angle_w_x(int i, int j) const {
        return phi[reduced_index(i,j)];
    }
    inline double angle_w_x(std::pair<int,int> pair_){
        return angle_w_x(pair_.first, pair_.second);
    }
};
#endif //MCPP_XY_HAMILTONIAN_DIPOLE_H
