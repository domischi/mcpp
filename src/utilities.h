#ifndef MCPP_UTILITIES_H_
#define MCPP_UTILITIES_H_

#include <alps/lattice.h>
#include <random>
#include "hamiltonians/hamiltonian_list.h"
namespace mcpp{
    typedef typename alps::Parameters parameter_type;
    typedef typename alps::graph_helper<> graph_helper_type;
    typedef typename alps::graph_helper<>::vector_type vector_type;
    typedef typename alps::graph_helper<>::basis_vector_iterator basis_vector_iterator;
    typedef typename alps::graph_helper<>::site_iterator site_iterator;
    typedef std::vector<std::vector<int>> neighbour_list_type;
    typedef std::vector<std::vector<vector_type>> difference_vector_list_type;
    typedef double spin_t;
    typedef std::vector<double> vector_t;
    
    inline double norm(vector_type const& x){
        double s=0;
        for (auto& e: x) {
            s+=e*e;
        }
        return std::sqrt(s);
    }
    inline double dot(vector_type const& x, vector_type const& y) {
        double s=0;
        for(int i = 0; i< x.size();++i)
            s+=x[i]*y[i];
        return s;
    } 
    std::pair<double,double> magnetization_and_staggered_magnetization(const std::vector<spin_t>& spins, const int L) {
        double mx =0.;
        double my =0.;
        double mxs=0.;
        double mys=0.;
        //at least theoretically vectorizable, however should not be the bottleneck
        for(int i=0;i<spins.size();++i) {
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
        return std::make_pair(std::sqrt(mx*mx+my*my)/spins.size(),std::sqrt(mxs*mxs+mys*mys)/spins.size());
    }

    double spin_config_overlap(std::vector<spin_t> const& spins1, std::vector<spin_t> const& spins2){
        int N=spins1.size();
        if(N!=spins2.size()){
            std::cerr<<"Trying to calculate autocorrelation between different sized arrays. Aborting..." <<std::endl;
            std::exit(46);
        }
        double acs=0.; 
        for(int j=0;j<N;++j)
            acs+=std::cos(spins1[j]-spins2[j]);
        return acs/N;
    }

    vector_t get_one_dipolar_field(vector_t const& r_ij, vector_t const& s_j, const double D) {
        double n=norm(r_ij);
        double r_ij_dot_s_j=dot(r_ij,s_j);
        vector_type h(r_ij.size());
        for(int i = 0; i<h.size();++i){ //this version works in 2D as well as in 3D
            h[i]=(D/(2.*std::pow(n,3)))*(-s_j[i]+3*r_ij[i]*r_ij_dot_s_j/(n*n));
        }
        return h;
    }
    inline vector_t get_one_dipolar_field(vector_t const& r_ij, double const& s, const double D) {
        return get_one_dipolar_field(r_ij, vector_t{std::cos(s),std::sin(s), 0.}, D);
    }
    inline vector_t get_one_dipolar_field(vector_t const& r_i, vector_t const& r_j, double const& s, const double D) {
        vector_t r_ij=r_i;
        for(int i =0;i<r_i.size();++i) r_ij[i]-=r_j[i];
        return get_one_dipolar_field(r_ij, s, D);
    }
   
    // Calculates the field of a dot under the assumption of a finite size, where each subsection however points in the same direction
    // diameter splits defines the number subparts considered if one crosses the diameter of a dot. I.e. diameter_split=1 considers only the dot itself
    // diameter_split=3 splits the magnetic moment into 5 sections arragned in a cross type configuration
    inline vector_t get_field_one_dot(const vector_t coor, const vector_t xyz, const double s, const double M, const int diameter_split, const double radius) {
        int counter=0;
        vector_t ret={0.,0.,0.};
        for(int i=0; i<diameter_split; ++i)
        for(int j=0; j<diameter_split; ++j){
            double offset_x=i*2.*radius/diameter_split+radius/diameter_split-radius;
            double offset_y=j*2.*radius/diameter_split+radius/diameter_split-radius;
            vector_t c= {coor[0]-offset_x, coor[1]-offset_y, coor[2]};
            vector_t tmp=get_one_dipolar_field(c, xyz, s, 1.); // as the magnetization gets split up, I will do this afterwards
            ret[0]+=tmp[0];
            ret[1]+=tmp[1];
            ret[2]+=tmp[2];
            ++counter;
        }
        ret[0]*=M/counter;
        ret[1]*=M/counter;
        ret[2]*=M/counter;
        return ret;
    }

    // Computes the field at a position xyz due to the spins on the lattice with PBC and finite sized dots 
    vector_t get_field_at_position(std::shared_ptr<alps::graph_helper<>> gh_, const std::vector<spin_t>& s, vector_t& xyz, std::vector<vector_t> const& coordinates, std::vector<bool> const& is_deleted, std::vector<vector_type> pt_, const double cutoff = 3., const int diameter_split=1, const double radius=0.3, const double M=1.) {
        vector_t h={0,0,0};
        for(site_iterator site=gh_->sites().first; site!=gh_->sites().second; ++site){
            if(!is_deleted[*site]){
                vector_t c_before_pb=coordinates[*site];
                for(const auto p : pt_){
                    vector_t c=c_before_pb;
                    for(int dim=0; dim<gh_->dimension();++dim) c[dim]+=p[dim]; 
                    if(std::pow(c[0]-xyz[0],2)+std::pow(c[1]-xyz[1],2)+std::pow(xyz[2],2)<cutoff*cutoff){
                        vector_t coor{c[0],c[1],0.};//need a 3d vector here
                        vector_t tmp=get_field_one_dot(coor, xyz, s[*site], M, diameter_split, radius);
                        h[0]+=tmp[0];
                        h[1]+=tmp[1];
                        h[2]+=tmp[2];
                    }
                }
            }
        }
        return h;
    }

    inline double maxfield(const double M, const double z, const double radius, const int diameter_split) {
        // The 3 in front seems more to be of a phenomenological structure taking into account the neighbours without really calculating them
        return 3*norm(get_field_one_dot(vector_type {0.,0.,0.} , vector_type{0.,0.,z}, 0., M, diameter_split,radius)); 
    }
    int init_N (const alps::Parameters& p){
        return alps::graph_helper<>(p).num_sites();
    }
    inline double init_T(const alps::Parameters& params) {
        if(params["ALGORITHM"]=="xy")
            return params.defined("T") ? static_cast<double>(params["T"]) : 1./static_cast<double>(params["beta"]);
        else if(params["ALGORITHM"]=="exmc")
            return 1.;// just set it to a value, it anyway gets set afterwards
        else {
            std::cerr<< "Cannot initalize T in model xy, because the algorithm is not detected. Aborting..."<<std::endl;
            std::exit(44);
        }
    }
    inline bool is_disordered(parameter_type const& p) {
        return 
            static_cast<double>(p.value_or_default("Position Disorder",0.))>0.||
            static_cast<double>(p.value_or_default("Dilution Rate"    ,0.))>0.||
            static_cast<bool>  (p.value_or_default("Enforce Disorder",false));
    }
    inline vector_type difference_vector(vector_type const& x, vector_type const& y, vector_type const& periodic) {
        vector_type diff=x;
        for(int i=0;i<x.size();++i) {
            diff[i]-=y[i]-periodic[i];
        }
        return diff; 
    }
    inline vector_type difference_vector(vector_type const& x, vector_type const& y) {
        return difference_vector(x,y,vector_type(x.size()));
    }
    inline double distance(vector_type const& x, vector_type const& y, vector_type const& periodic){
        return norm(difference_vector(x,y,periodic));
    }
    std::vector<double> positional_disorder(graph_helper_type const& gh, parameter_type p, std::mt19937& rng) {
        std::vector<double> dR(gh.dimension()*gh.num_sites());
        std::normal_distribution<double> dist(0,p["Position Disorder"]);
        for(int i=0;i<gh.dimension()*gh.num_sites();++i) dR[i]=dist(rng);
        return dR;
    }
    std::set<int> dilution_disorder(graph_helper_type const& gh, parameter_type p, std::mt19937& rng) {
        double dilution_rate=p["Dilution Rate"];
        std::cout << "WARNING: dilution is only implemented for the dipolar interaction!"<<std::endl;
        std::vector<int> site_list(gh.num_sites());
        iota(site_list.begin(),site_list.end(),0);// List of the sites
        std::shuffle(site_list.begin(),site_list.end(),rng);
        int Num_Diluted_Sites=static_cast<int>(std::round(dilution_rate*gh.num_sites()));
        std::cout << "Remove "<<Num_Diluted_Sites<<"/"<<gh.num_sites()<<" resulting in an effective dilution rate of "<<std::setprecision(3)<<Num_Diluted_Sites*1.0/gh.num_sites()<< " ~ "<<dilution_rate<<std::endl;
        site_list.resize(Num_Diluted_Sites);
        return std::set<int>(site_list.begin(),site_list.end()); 
    }
    std::vector<vector_type> get_periodic_translations(graph_helper_type const& gh, parameter_type const& p) {
        int L=p["L"];
        std::vector<vector_type> basis;
        basis_vector_iterator v, v_end;
        for(std::tie(v,v_end)=gh.basis_vectors();v!= v_end;++v){
            basis.push_back(*v);
        }
        std::vector<vector_type> periodic_translations;
        int dim = gh.dimension();
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
        return periodic_translations;
    }

    std::pair<neighbour_list_type, std::vector<std::vector<vector_type>>> sort_neighbours(neighbour_list_type const& neighbours_in, std::vector<std::vector<vector_type>> const& diff_vectors_in) {
        neighbour_list_type neighbour_out=neighbours_in;
        std::vector<std::vector<vector_type>> diff_vectors_out=diff_vectors_in;
        for(int i =0 ; i<neighbour_out.size();++i) {
            std::vector<std::pair<int,vector_type>> data;
            for(int j = 0; j<neighbour_out[i].size();++j) {
                data.push_back(std::make_pair(neighbours_in[i][j],diff_vectors_in[i][j]));
            }
            std::sort(data.begin(),data.end(),[](const std::pair<int,vector_type>& a, const std::pair<int,vector_type>& b) -> bool {return a.first<b.first; }); //c++14 could make this easier with the auto's
            for(int j = 0; j<neighbour_out[i].size();++j) {
                neighbour_out[i][j]=data[j].first;
                diff_vectors_out[i][j]=data[j].second;
            }
        }
        return std::make_pair(neighbour_out, diff_vectors_out);
}

    // Returns:
    // first: vector if site is deleted
    // second: vector describing the coordinate
    std::pair<std::vector<bool>,std::vector<vector_type>> get_coordinates(graph_helper_type const& gh, parameter_type const& p){
        const double dilution_rate=p.value_or_default("Dilution Rate", 0.);
        const double position_std_dev=static_cast<double>(p.value_or_default("Position Disorder", 0.))*static_cast<double>(p.value_or_default("a",1.));
        const int DISORDER_SEED=p["DISORDER_SEED"];
        const bool DISORDERED=is_disordered(p);
        const int L=p["L"];
        const int N=gh.num_sites();
        std::mt19937 rng(DISORDER_SEED);
        //Position disorder 
        std::vector<double> dR(gh.dimension()*gh.num_sites());
        if(DISORDERED&&position_std_dev>0.)
            dR=mcpp::positional_disorder(gh, p, rng);
        //Dilution disorder
        std::set<int> deleted_sites;
        if(DISORDERED&&dilution_rate>0)
            deleted_sites=mcpp::dilution_disorder(gh,p,rng);
        std::vector<bool> is_deleted(N);
        for(site_iterator s = gh.sites().first; s !=gh.sites().second; ++s){
            if(deleted_sites.find(*s)!=deleted_sites.end()) // *s is deleted
                is_deleted[*s]=true;
        }
        std::vector<vector_type> disordered_coordinates(N);
        vector_type tmp;
        for(site_iterator s = gh.sites().first; s !=gh.sites().second; ++s){
            tmp=gh.coordinate(*s);
            for(int d = 0 ; d<gh.dimension();++d) {
                tmp[d]+=dR[(*s)*gh.dimension()+d];
            }
            disordered_coordinates[*s]=tmp;
        }
        return std::make_pair(is_deleted, disordered_coordinates); 
    }

    std::pair<neighbour_list_type, difference_vector_list_type> get_neighbours(graph_helper_type const& gh, parameter_type const& p){
        const bool print_debug_information=static_cast<bool>(p.value_or_default("debug dipolar",false));
        const int L=p["L"];
        const int N=gh.num_sites();
        const double cutoff_distance = static_cast<double>(p.value_or_default("cutoff_distance",3))
                                      *std::max(static_cast<double>(p.value_or_default("a",1.)),static_cast<double>(p.value_or_default("b",1.))); 
        neighbour_list_type neighbour_list(N);
        for(int i=0;i<gh.num_sites();++i) neighbour_list[i]=std::vector<int>();
        std::vector<std::vector<vector_type>> difference_vector_list(N);
        const std::vector<vector_type> periodic_translations=get_periodic_translations(gh,p);
        std::map<int,double> dist_map;
        std::vector<bool> is_deleted;
        std::vector<vector_type> coordinates;
        std::tie(is_deleted, coordinates)=get_coordinates(gh, p);
        for(site_iterator s_iter = gh.sites().first; s_iter !=gh.sites().second; ++s_iter )
        for(site_iterator s_iter2= gh.sites().first; s_iter2!=gh.sites().second; ++s_iter2)
        for(auto& p : periodic_translations)
            if(*s_iter!=*s_iter2 && //same site
                !is_deleted[*s_iter ] && //s_iter  is not deleted
                !is_deleted[*s_iter2] ){ //s_iter2 is not deleted
                vector_type c1(coordinates[*s_iter ]);
                vector_type c2(coordinates[*s_iter2]);
                if(print_debug_information && *s_iter2==1 && p==periodic_translations[0]){
                    std::cout <<"Site: " <<std::setw(3) <<*s_iter<<" with vector ("
                        << std::setw(8) << c1[0]<<","
                        << std::setw(8) << c1[1]<<")"<<std::endl;
                }
                const vector_type diff_vec=mcpp::difference_vector(c1,c2,p);
                const double dist=norm(diff_vec);
                const int ri = (*s_iter)+N*(*s_iter2);
                if(dist<=cutoff_distance){
                    if(dist_map[ri]==0.){
                        dist_map[ri]=dist;
                        neighbour_list[*s_iter].push_back(*s_iter2);
                        difference_vector_list[*s_iter].push_back(diff_vec);
                    }
                    else {
                        //WRONG distmap has an entry then there is an interaction with a close particle in different fashions, therefore abort
                        std::cerr 
                            << "There is an error here! dist_map already has a value for this pair. "<<std::endl
                            << "s_iter =" <<std::setw(4)<<*s_iter
                            << "s_iter2=" <<std::setw(4)<<*s_iter2
                            << "old distance = "<<std::setw(10) << dist_map[ri]
                            << "new distance = "<<std::setw(10) << dist<<std::endl;
                            std::exit(34);
                    }
                }
            }
        return sort_neighbours(neighbour_list, difference_vector_list);
    }
}
#endif //MCPP_UTILITIES_H_
