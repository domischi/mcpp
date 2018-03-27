#ifndef MCPP_MCRG_UTILITIES_H_
#define MCPP_MCRG_UTILITIES_H_

#include <vector>
#include <cmath>
#include <string>
#include <alps/parameter.h>
#ifndef NARMADILLO
    #include <armadillo>
#endif //NARMADILLO
namespace mcpp{
namespace mcrg{
    typedef double spin_t;
    typedef std::tuple<int,std::vector<int>> shift_t; //component, (dx,dy,dz,...)
    enum LatticeType {SquareLatticeIsh, CubicLatticeIsh};
    spin_t inline mod2Pi(double s){
        return s-std::floor(s/(2*M_PI))*2*M_PI;
    }
    LatticeType init_lattice_type(std::string s) {
        if(s=="square lattice" or s=="anisotropic square lattice") {
            return LatticeType::SquareLatticeIsh;
        }
        else if(s=="simple cubic lattice" or s=="anisotropic simple cubic lattice") {
            return LatticeType::CubicLatticeIsh;
        }
        else {
            std::cerr << "Didn't find the latticetype for mcrg, check if your lattice is already implemented for MCRG"<<std::endl;
            std::exit(4);
        }
    }
    LatticeType init_lattice_type(alps::Parameters const& p) {
        return init_lattice_type(static_cast<std::string>(p["LATTICE"]));
    }
    std::vector<std::vector<shift_t>> init_interactions(std::string s) {
        std::vector<std::vector<shift_t>> i;
        //Choose Interaction Set
        if(s=="small"){
            i=mcpp::mcrg::small;
        } else if(s=="minimal"){
            i=mcpp::mcrg::minimal;
        } else if(s=="minimal2"){
            i=mcpp::mcrg::close_to_minimal;
        } else if(s=="very small"){
            i=mcpp::mcrg::very_small;
        } else if(s=="medium interaction"){
            i=mcpp::mcrg::medium_interaction;
        } else if(s=="medium range" || s=="medium" ){
            i=mcpp::mcrg::medium_range;
        } else if(s=="dXY" ){
            i=mcpp::mcrg::dXY_handcrafted;
        } else if(s=="dIsing" ){
            i=mcpp::mcrg::dIsing;
        } else if(s=="wilson"){
            i=mcpp::mcrg::wilson1975;
        } else if(s=="ising"){
            i=mcpp::mcrg::small_Ising;
        } else if(s=="xy-3d-small"){
            i=mcpp::mcrg::xy_3d_small;
        } else if(s=="xy-3d-very-small"){
            i=mcpp::mcrg::xy_3d_very_small;
        } else if(s=="xy-3d-medium"){
            i=mcpp::mcrg::xy_3d_medium;
        } else if(s=="xy-3d-massive"){
            i=mcpp::mcrg::xy_3d_massive;
        } else if(s=="xy-3d-swendsen"){
            i=mcpp::mcrg::xy_3d_swendsen;
        } else {
            std::cerr<<"Did not recognise the interaction set to use for MCRG, aborting..."<<std::endl;
            std::exit(20);
        }
        for(auto& in : i) in.shrink_to_fit();
        i.shrink_to_fit();
        return i;
    }
    std::vector<std::vector<shift_t>> init_interactions(alps::Parameters const& p){
        return init_interactions(static_cast<std::string>(p["MCRG Interactions"]));
    }
   	inline int index_o_neighbour_square (int i, int dx,int dy, int L){
        int x,y;
		int N=L*L;
        i+=N; dx+=N; dy+=N;
        x=(i/L+dx)%L;
        y=(i%L+dy)%L;
        assert(x<L&&y<L);
        return x+L*y;
    }
   	inline int index_o_neighbour_cubic (int i, int dx,int dy, int dz, int L){
        int x,y,z;
		int N=L*L*L;
        i+=N; dx+=N; dy+=N; dz+=N;
        x=(i/(L*L)+dx)%L;
        y=(i/L+dy)%L;
        z=(i%L+dz)%L;
        assert(x<L&&y<L&&z<L);
        return z+L*(y+L*x);
    }

    inline int index_o_neighbour_L_helper(int const& size, LatticeType const& l){
        switch(l){
            case LatticeType::SquareLatticeIsh:
                return nearest_integer(std::pow(size,1./2));
            case LatticeType::CubicLatticeIsh:
                return nearest_integer(std::pow(size,1./3));
        }
    }

    inline int index_o_neighbour (const int& i, std::vector<int> const& dxdydz, LatticeType const& lattice_type, int const& L) {
        switch(lattice_type){
            case LatticeType::SquareLatticeIsh:
                return mcpp::mcrg::index_o_neighbour_square(i,dxdydz[0],dxdydz[1],L);
            case LatticeType::CubicLatticeIsh:
                return mcpp::mcrg::index_o_neighbour_cubic(i,dxdydz[0],dxdydz[1],dxdydz[2],L);
        }
    }

    //This generally calculates an the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    double S_alpha (const std::vector<spin_t>& spins, std::vector<shift_t> shifts, LatticeType const& lattice_type) {
        double S_a=0;
        const int L=index_o_neighbour_L_helper(spins.size(), lattice_type);
        for(int i =0;i<spins.size();++i){
            double tmp=1.;
            for(int j = 0; j <shifts.size();++j){ //Pairwise build up the scalar products
                shift_t s = shifts[j];
                const double angle=spins[index_o_neighbour(i,std::get<1>(s),lattice_type,L)];
                const int component=std::get<0>(s);
                if(component==0){
                    tmp*=std::cos(angle);
                } else if(component==1){
                    tmp*=std::sin(angle);
                } else {
                    std::cerr << "In S_alpha with a not recognised component, The values is "<< component<<"... Aborting..."<< std::endl;
                    std::exit(19);
                }
            }
            S_a+=tmp;
        }
        return S_a;
    }

    //This calculates all the S_alphas and returns them separately for the two sectors in vectors
    std::pair<std::vector<double>,std::vector<double>> all_S_alphas(std::vector<spin_t> spins, std::vector<std::vector<shift_t>> const& interactionSet, LatticeType const& lattice_type) {
        std::vector<double> oute, outo;
        for(int interaction=0;interaction<interactionSet.size();++interaction) {
            if(interactionSet[interaction].size()%2) //odd
                outo.push_back(S_alpha(spins, interactionSet[interaction], lattice_type));
            else //even
                oute.push_back(S_alpha(spins, interactionSet[interaction], lattice_type));
        }
        return std::make_pair(oute,outo);
    }

    inline int new_entry_point(int entry_point, const int L, const int scale_factor_b){
        if((entry_point%L)%scale_factor_b==scale_factor_b-1) {
            entry_point+=L-scale_factor_b;
        }
        ++entry_point;
        if(entry_point==scale_factor_b*L) {
            entry_point=0;
        }
        return entry_point;
    }

    #ifndef NARMADILLO
    inline arma::mat dSdK(std::valarray<double> const& Data_Sa, std::valarray<double> const& Data_Sb, std::valarray<double> const& Data_SaSb){
        arma::vec Sa   = arma::vec(&Data_Sa  [0],Data_Sa.size());
        arma::vec Sb   = arma::vec(&Data_Sb  [0],Data_Sa.size());
        arma::mat SaSb = arma::mat(&Data_SaSb[0],Data_Sa.size(),Data_Sa.size());
        SaSb-=Sa*Sb.t();
        return SaSb;
    }
    arma::mat get_T(std::valarray<double> const& Data_Sa, std::valarray<double> const& Data_Sap, std::valarray<double> const& Data_SaSb, std::valarray<double> const& Data_SaSbp){
        assert(Data_Sa.size()==Data_Sap.size() &&Data_Sa.size()*Data_Sa.size()==Data_SaSb.size() &&Data_Sa.size()*Data_Sa.size()==Data_SaSbp.size());
        arma::mat SaSb =dSdK(Data_Sa,Data_Sa ,Data_SaSb);
        arma::mat SaSbp=dSdK(Data_Sa,Data_Sap,Data_SaSbp);
        return SaSbp*arma::inv(SaSb);
    }
    #endif //NARMADILLO
}//mcrg
}//mcpp
#endif //MCPP_MCRG_UTILITIES_H_
