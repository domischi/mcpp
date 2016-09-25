#ifndef MCPP_XY_HAMILTONIAN_DIPOLE_H
#define MCPP_XY_HAMILTONIAN_DIPOLE_H

#include <vector>
#include <numeric>
#include <random>
#include <cmath>
#include <algorithm>
#include <iomanip>

#include "xy_hamiltonian.h"

class XY_Dipole : public XY_Hamiltonian{
public:
    //typedef typename XY_Hamiltonian<G>::spin_t spin_t;
    XY_Dipole (alps::Parameters const& params) : 
    XY_Hamiltonian(params),
    D(params.value_or_default("D",1)),
    L(params["L"]),
    N(L*L),
    dilution_rate(params.value_or_default("Dilution Rate", 0.)),
    gh(params){ 
        cutoff_distance=static_cast<double>(params.value_or_default("cutoff_distance",3))*static_cast<double>(params.value_or_default("a",1.));
        Init_Lookup_Tables(params["DISORDER_SEED"]);
        //for (int i = 0; i < N;++i)
        //for (int j = 0; j < N;++j)
        //    std::cout<< std::setw(3)<<i<<std::setw(3)<<j<<std::setw(3)<<reduced_index(i,j).first<<std::setw(3)<<reduced_index(i,j).second<<std::endl;
        //for(int i=0;i<N;++i)
        //    for(int j=0;j<N;++j)
        //        std::cout <<std::setw(4)<<i<<std::setw(4)<<j<<"\t("<<reduced_index(i,j).first<<", "<<reduced_index(i,j).second<<")\t("<<reduced_index(j,i).first<<","<<reduced_index(j,i).second<<")"<<std::endl;
    
    } 
    virtual double SingleSiteEnergy(std::vector<double> const& spins, int i){
        double e=0.;
        for(int j : neighbour_list[i]){
            e-=D*inv_distance_cubed(i,j)*(1.5*std::cos(spins[i]+spins[j]-2*angle_w_x(i,j))+0.5*std::cos(spins[i]-spins[j]));
        }
        return e;
    }

private:
    int L,N, DISORDER_SEED;
    double D;
    double dilution_rate;
    double cutoff_distance;
    alps::graph_helper<> gh;
    //Lookup tables 
    //std::unordered_map<int,double> dist_3; //TODO check if an unordered map is faster...
    //std::unordered_map<int,double> phi;
    std::vector<double> dist_3; //TODO check if an unordered map is faster...
    std::vector<double> phi;
    std::vector<std::vector<int>> neighbour_list; //saves which neighbours are relevant (as not only nearest neighbours count in dipole )
                        
    typedef typename alps::graph_helper<>::vector_type vector_type;
    typedef typename alps::graph_helper<>::basis_vector_iterator basis_vector_iterator;
    typedef typename alps::graph_helper<>::site_iterator site_iterator;
    
    inline int reduced_index(int i, int j){
        int xi,yi,xj,yj,x,y;
        xi=i/L;
        yi=i%L;
        xj=j/L;
        yj=j%L;
        x=(xi-xj+L)%L;
        y=(yi-yj+L)%L;
        return L*x+y;
    }
    inline int reduced_index(std::pair<int,int> p) {
        return reduced_index(p.first,p.second);
    }
    void Init_Lookup_Tables(int DISORDER_SEED){
        dist_3=std::vector<double>(N);
        phi=std::vector<double>(N);
        std::vector<vector_type> basis;
        basis_vector_iterator v, v_end;
        for(std::tie(v,v_end)=gh.basis_vectors();v!= v_end;++v){
            basis.push_back(*v);
        }
        std::vector<vector_type> periodic_translations;
        int dim = gh.dimension();
        assert(dim==basis[0].size());
        periodic_translations.push_back(vector_type(dim,0.));
        vector_type pmb(gh.dimension()),ppb(gh.dimension()); //periodic_translation+-L*basis
        for(auto& actual_basis_vector : basis){
            int size_periodic_translations=periodic_translations.size();
            for(int i=0;i<size_periodic_translations;++i){//to avoid double counting a specific vector, and to not have a segfault
                vector_type p = periodic_translations[i];
                for(int d =0;d<gh.dimension();++d){
                    pmb[d]=(p[d]-L*(actual_basis_vector[d]));
                    ppb[d]=(p[d]+L*(actual_basis_vector[d]));
                }
                periodic_translations.push_back(pmb);
                periodic_translations.push_back(ppb);
            }
        }
        for(int i=0;i<gh.num_sites();++i) neighbour_list.push_back(std::vector<int>());
        std::map<int,double> dist_map;
        for(site_iterator s_iter = gh.sites().first; s_iter !=gh.sites().second; ++s_iter )
        for(site_iterator s_iter2= gh.sites().first; s_iter2!=gh.sites().second; ++s_iter2)
        for(auto& p : periodic_translations)
            if(*s_iter!=*s_iter2){
                vector_type c1(gh.coordinate(*s_iter));
                vector_type c2(gh.coordinate(*s_iter2));
                double dist=distance(c1,c2,p);
                int ri = reduced_index(*s_iter, *s_iter2);
                //if(*s_iter==8 && *s_iter2==5) std::cout << dist<<std::setw(15)<<reduced_index(8,5).first<<" "<<reduced_index(8,5).second<<std::endl;
                if(dist<=cutoff_distance && (dist_map[ri]==0. || dist<=dist_map[ri])){ //new value is smaller than the previous calculated for this pair
                    dist_map[ri]=dist;
                    dist_3[ri]=std::pow(dist,-3);
                    phi[ri]=std::atan2(-(c1[1]+p[1]-c2[1]),-(c1[0]+p[0]-c2[0]));
                    neighbour_list[*s_iter].push_back(*s_iter2);
                }
            }
        //Dilution Disorder //TODO this is broken as with the reduced_index business...
        if(dilution_rate > 0.){
            std::cout << "WARNING: dilution is only implemented for the dipolar interaction!"<<std::endl;
            std::mt19937 rng(DISORDER_SEED);
            std::vector<int> site_list(N);
            iota(site_list.begin(),site_list.end(),0);// List of the sites
            std::shuffle(site_list.begin(),site_list.end(),rng);
            int Num_Diluted_Sites=static_cast<int>(std::round(dilution_rate*N));
            std::cout << "Remove "<<Num_Diluted_Sites<<"/"<<N<<" resulting in an effective dilution rate of "<<std::setprecision(3)<<Num_Diluted_Sites*1.0/N<< " ~ "<<dilution_rate<<std::endl;
            site_list.resize(Num_Diluted_Sites);
            for (auto& l : site_list) {
                neighbour_list[l]=std::vector<int>(); //make sure these sites have no neighbours 
            }
            for(auto& nl : neighbour_list){ //make sure that they are not neighbours of others
                for(auto& s : site_list)
                    nl.erase(std::remove(nl.begin(),nl.end(),s),nl.end());
            }
        }
        //Sort & Shrink the neighbour_list
        for(auto& nl : neighbour_list) {
            std::sort(nl.begin(),nl.end());
            nl.shrink_to_fit();
        }
        neighbour_list.shrink_to_fit();
        //for(auto& nl:neighbour_list) std::cout << nl.size()<<std::endl;// not all of them have the same number of neighbours...
        //int test=8;
        //for (auto& s : neighbour_list[test]){
        //    std::cout << test<< std::setw(4)<<s <<std::setw(15)<<std::pow(inv_distance_cubed(test,s),(-2./3))<<std::setw(15)<<angle_w_x(test,s)/M_PI<<std::endl;
        //}
    }
    inline double distance(vector_type& x, vector_type& y, vector_type& periodic){
        return std::sqrt(std::pow(x[0]-y[0]+periodic[0],2)+std::pow(x[1]-y[1]+periodic[1],2));
    }
    inline double inv_distance_cubed(int i,int j) {
        return inv_distance_cubed(std::make_pair(i,j));
    }
    inline double inv_distance_cubed(std::pair<int,int> pair_){// TODO make the packing and packing of a pair consistent
        double ret_val=dist_3[reduced_index(pair_)];
        assert(ret_val>0.);
        return ret_val;
    }
    inline double angle_w_x(int i, int j) {//it might be useful to check what it is for site 0 and then do the transformation for the second site if one assumes a shift first. This would allow only to save the interactions for the 0th site, and therefore reduce the map by a Factor of N which will probably be easier to save & load 
        return angle_w_x(std::make_pair(i,j));
    }
    inline double angle_w_x(std::pair<int,int> pair_){
        return phi[reduced_index(pair_)];
    }

};

#endif //MCPP_XY_HAMILTONIAN_DIPOLE_H

