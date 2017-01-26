#ifndef MCPP_UTILITIES_H_
#define MCPP_UTILITIES_H_

#include <alps/lattice.h>
#include <random>
namespace mcpp{
    typedef typename alps::Parameters parameter_type;
    typedef typename alps::graph_helper<> graph_helper_type;
    typedef typename alps::graph_helper<>::vector_type vector_type;
    typedef typename alps::graph_helper<>::basis_vector_iterator basis_vector_iterator;
    typedef typename alps::graph_helper<>::site_iterator site_iterator;
    typedef std::vector<std::vector<int>> neighbour_list_type;
   
    int init_N (const alps::Parameters& p){
        return alps::graph_helper<>(p).num_sites();
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
    inline double norm2(vector_type const& x){
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
    inline double distance(vector_type const& x, vector_type const& y, vector_type const& periodic){
        return norm2(difference_vector(x,y,periodic));
    }

    std::vector<double> positional_disorder(graph_helper_type const& gh, parameter_type p, std::mt19937 rng) {
        std::vector<double> dR(gh.dimension()*gh.num_sites());
        std::normal_distribution<double> dist(0,p["Position Disorder"]);
        for(int i=0;i<gh.dimension()*gh.num_sites();++i) dR[i]=dist(rng);
        return dR;
    }

    std::set<int> dilution_disorder(graph_helper_type const& gh, parameter_type p, std::mt19937 rng) {
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
            std::sort(data.begin(),data.end(),[](auto& a, auto& b) -> bool {return a.first<b.first; });
            for(int j = 0; j<neighbour_out[i].size();++j) {
                neighbour_out[i][j]=data[j].first;
                diff_vectors_out[i][j]=data[j].second;
            }
        }
        return std::make_pair(neighbour_out, diff_vectors_out);
}

    std::pair<neighbour_list_type, std::vector<std::vector<vector_type>>> get_neighbours(graph_helper_type const& gh, parameter_type p){
        const bool print_debug_information=static_cast<bool>(p.value_or_default("debug dipolar",false));
        const double dilution_rate=p.value_or_default("Dilution Rate", 0.);
        const double position_std_dev=static_cast<double>(p.value_or_default("Position Disorder", 0.))*static_cast<double>(p.value_or_default("a",1.));
        const int DISORDER_SEED=p["DISORDER_SEED"];
        const bool DISORDERED=is_disordered(p);
        const int L=p["L"];
        const int N=gh.num_sites();
        const double cutoff_distance = static_cast<double>(p.value_or_default("cutoff_distance",3))
                                      *std::max(static_cast<double>(p.value_or_default("a",1.)),static_cast<double>(p.value_or_default("b",1.))); 
        neighbour_list_type neighbour_list(N);
        for(int i=0;i<gh.num_sites();++i) neighbour_list[i]=std::vector<int>();
        std::vector<std::vector<vector_type>> difference_vector_list(N);
        const std::vector<vector_type> periodic_translations=get_periodic_translations(gh,p);
        std::map<int,double> dist_map;
        std::mt19937 rng(DISORDER_SEED);
        //Position disorder 
        std::vector<double> dR(gh.dimension()*gh.num_sites());
        if(DISORDERED&&position_std_dev>0.)
            dR=mcpp::positional_disorder(gh, p, rng);
        //Dilution disorder
        std::set<int> deleted_sites;
        if(DISORDERED&&dilution_rate>0)
            deleted_sites=mcpp::dilution_disorder(gh,p,rng);
        for(site_iterator s_iter = gh.sites().first; s_iter !=gh.sites().second; ++s_iter )
        for(site_iterator s_iter2= gh.sites().first; s_iter2!=gh.sites().second; ++s_iter2)
        for(auto& p : periodic_translations)
            if(*s_iter!=*s_iter2 && //same site
               deleted_sites.find(*s_iter) ==deleted_sites.end() && //s_iter  is not deleted
               deleted_sites.find(*s_iter2)==deleted_sites.end() ){ //s_iter2 is not deleted
                vector_type c1(gh.coordinate(*s_iter));
                vector_type c2(gh.coordinate(*s_iter2));
                for(int d=0;d<gh.dimension();++d){
                    c1[d]+=dR[*s_iter +d];
                    c2[d]+=dR[*s_iter2+d];
                }
                if(print_debug_information && *s_iter2==1 && p==periodic_translations[0]){
                    vector_type c(gh.coordinate(*s_iter));
                    std::cout <<"Site: " <<std::setw(3) <<*s_iter<<" with vector ("
                        << std::setw(8) << c1[0]<<","
                        << std::setw(8) << c1[1]<<")"<<std::endl;
                }
                const vector_type diff_vec=mcpp::difference_vector(c1,c2,p);
                const double dist=norm2(diff_vec);
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
