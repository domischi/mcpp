#ifndef MCPP_MCRG_H_
#define MCPP_MCRG_H_

#include <alps/lattice.h>
#include <cassert>
#include <valarray>
#include <alps/parameter.h>
#include <tuple>
#include "interactions.h"
#include "reductions.h"
#include "../../utilities.h"

class mcrg; //Needed as it is recursively called
class mcrg {
public:
    typedef double spin_t;
    typedef mcrg_utilities::shift_t shift_t; //(dx,dy, component)
    enum ReductionTechnique {Decimation,Blockspin, FerroBlockspin, Blockspin4x4, FerroBlockspin4x4, IsingMajority, IsingTieBreaker, BlockspinCubic, DecimationCubic};
    enum LatticeType {SquareLatticeIsh, CubicLatticeIsh};
    
    mcrg(const alps::Parameters& p, int Iteration_, int MCRG_It_) :
    iteration(Iteration_),
    max_iterations(MCRG_It_),
    L(p["L"]),
    N(mcpp::init_N(p)),
    lattice_type(init_lattice_type(static_cast<std::string>(p["LATTICE"]))),
    reduction_type(init_reduction_technique(static_cast<std::string>(p["MCRG Reduction Technique"]))),
    scale_factor_b(init_b(reduction_type)),
    entry_point(N-1),
    interactions(init_interactions(static_cast<std::string>(p["MCRG Interactions"])))
    {
        //If this is not the last iteration, construct the MCRG class for the next smaller lattice 
        if(!is_last_iteration()) {
            alps::Parameters p_descendant=p;
            p_descendant["L"]=int(L/scale_factor_b); //probably should check for the sqrt updates
            descendant=std::make_shared<mcrg>(p_descendant,iteration+1,MCRG_It_);
        }
    }
                
    std::tuple<std::valarray<double>,std::valarray<double>> measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs){
        ++entry_point%=N;//loop over the sites as entry points //TODO implement a true loop...
        //std::valarray<double> OUT(n_interactions());
        std::vector<double> oute, outo;
        int counto=0;
        int counte=0;
        //measure the S_alpha in the not reduced lattice
        for(int i=0;i<n_interactions();++i) {
            if(interactions[i].size()%2) { //odd
                outo.push_back(S_alpha(spins, interactions[i]));
                ++counto;
            }
            else{
                oute.push_back(S_alpha(spins, interactions[i]));
                ++counte;
            }
        }
        std::valarray<double> OUTe(oute.data(), oute.size()), OUTo(outo.data(), outo.size());
        
        std::valarray<double> save_outo=OUTo;
        //save_outo/=N;
        std::valarray<double> save_oute=OUTe;
        //save_oute/=N;
        obs["MCRGe S_alpha"+ std::to_string(iteration)]<<save_oute;
        obs["MCRGo S_alpha"+ std::to_string(iteration)]<<save_outo;
        //measure <S_alpha n S_beta n>
        std::valarray<double> outouto(OUTo.size()*OUTo.size());
        std::valarray<double> outoute(OUTe.size()*OUTe.size());
        
        for(int i=0;i<counto;++i){
            for(int j=0;j<counto;++j){
                outouto[i*counto+j]=OUTo[i]*OUTo[j];
            }
        }
        for(int i=0;i<counte;++i){
            for(int j=0;j<counte;++j){
                outoute[i*counte+j]=OUTe[i]*OUTe[j];
            }
        }
        //outoute/=N;
        //outouto/=N;
        obs["MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outouto;
        obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outoute;
        
        if(!is_last_iteration()){//This is not the last instantation of the mcrg
            //reduce the lattice and call the descendant on the renormalized system
            std::vector<spin_t> reduced_spins = reduce(spins, entry_point);
            std::valarray<double> INo(OUTo.size()),INe(OUTe.size());
            std::tie(INo, INe) =  descendant->measure(reduced_spins,obs);
            //calculate <S_alpha n-1 S_beta n>
            //and save it into obs
            std::valarray<double> outine(INe.size()*INe.size());
            for(int i=0;i<INe.size();++i){
                for(int j=0;j<INe.size();++j){
                    outine[i*INe.size()+j]=(OUTe[i]*INe[j]);
                }
            }
            //outine/=N;
            obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outine;
            std::valarray<double> outino(INo.size()*INo.size());
            for(int i=0;i<INo.size();++i){
                for(int j=0;j<INo.size();++j){
                    outino[i*INo.size()+j]=(OUTo[i]*INo[j]);
                }
            }
            //outino/=N;
            obs["MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outino;
        }
        return std::tie(OUTo,OUTe);
    }

    void init_observables(alps::ObservableSet& obs){
        obs << alps::RealVectorObservable("MCRGe S_alpha"+ std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1));
        obs << alps::RealVectorObservable("MCRGo S_alpha"+ std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1));
        if(!is_last_iteration()){
            descendant->init_observables(obs);
        }
    }
    void save(alps::ODump &dump) const{
        dump
            << entry_point;
    }
    void load(alps::IDump &dump) {
        dump 
            >> entry_point;
    }
private:
    std::shared_ptr<mcrg> descendant;
    const int iteration;
    const int max_iterations; 
    const int L,N;
    const LatticeType lattice_type;
    const ReductionTechnique reduction_type;
    //first index: which interaction
    //second index: listing the vectors of the shift
    //third index (shift_t): a pair of ints denoting the shift wrt to the first one (always (0,0))
    const std::vector<std::vector<shift_t>> interactions;
    const int scale_factor_b;
    int entry_point;
    ReductionTechnique init_reduction_technique(std::string s) const {
        if(s=="Decimation"){
            if(lattice_type==LatticeType::SquareLatticeIsh)
                return ReductionTechnique::Decimation;
            if(lattice_type==LatticeType::CubicLatticeIsh)
                return ReductionTechnique::DecimationCubic;
        } else if(s=="Blockspin"){
            if(lattice_type==LatticeType::SquareLatticeIsh)
                return ReductionTechnique::Blockspin;
            if(lattice_type==LatticeType::CubicLatticeIsh)
                return ReductionTechnique::BlockspinCubic;
        } else if(s=="FerroBlockspin"){
            return ReductionTechnique::FerroBlockspin;
        } else if(s=="Blockspin4x4"){
            return ReductionTechnique::Blockspin4x4;
        } else if(s=="FerroBlockspin4x4"){
            return ReductionTechnique::FerroBlockspin4x4;
        } else if(s=="IsingMajority"){
            return ReductionTechnique::IsingMajority;
        } else if(s=="IsingTieBreaker"){
            return ReductionTechnique::IsingTieBreaker;
        } else {
            std::cerr << "Did not recognise the reduction process to use for MCRG, aborting..."<<std::endl;
            std::exit(21); 
        }
    }
    LatticeType init_lattice_type(std::string s) const {
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
    inline int n_interactions(){
        return interactions.size();
    }
    std::vector<std::vector<shift_t>> init_interactions(std::string s) const {
        std::vector<std::vector<shift_t>> i; 
        //Choose Interaction Set
        if(s=="small"){
            i=mcrg_utilities::small;
        } else if(s=="very small"){
            i=mcrg_utilities::very_small;
        } else if(s=="medium interaction"){
            i=mcrg_utilities::medium_interaction;
        } else if(s=="medium range" || s=="medium" ){
            i=mcrg_utilities::medium_range;
        } else if(s=="dXY" ){
            i=mcrg_utilities::dXY_handcrafted;
        } else if(s=="dIsing" ){
            i=mcrg_utilities::dIsing;
        //} else if(s=="dIsing2"){
        //    i=mcrg_utilities::dIsing2;
        } else if(s=="ising"){
            i=mcrg_utilities::small_Ising;
        } else if(s=="xy-3d-small"){
            i=mcrg_utilities::xy_3d_small;
        } else if(s=="xy-3d-very-small"){
            i=mcrg_utilities::xy_3d_very_small;
        } else if(s=="xy-3d-medium"){
            i=mcrg_utilities::xy_3d_medium;
        } else if(s=="xy-3d-massive"){
            i=mcrg_utilities::xy_3d_massive;
        } else if(s=="xy-3d-swendsen"){
            i=mcrg_utilities::xy_3d_swendsen;
        } else {
            std::cerr<<"Did not recognise the interaction set to use for MCRG, aborting..."<<std::endl;
            std::exit(20);
        }
        for(auto& in : i) in.shrink_to_fit();
        i.shrink_to_fit();
        return i;
    }
    int init_b(ReductionTechnique rt){
        switch (rt) {
            case ReductionTechnique::Decimation:
               return 2; 
            case ReductionTechnique::IsingTieBreaker:
               return 2;
            case ReductionTechnique::DecimationCubic:
               return 2; 
            case ReductionTechnique::Blockspin:
               return 2; 
            case ReductionTechnique::BlockspinCubic:
               return 2; 
            case ReductionTechnique::FerroBlockspin:
               return 2;
            case ReductionTechnique::IsingMajority:
               return 3;
            case ReductionTechnique::Blockspin4x4:
               return 4;
            case ReductionTechnique::FerroBlockspin4x4:
               return 4;
            default:
               std::cerr<< "Did not recognize the ReductionTechnique. Aborting...."<<std::endl;
               std::exit(4);
        }
    }
    bool inline is_last_iteration(){
        return max_iterations<=iteration;
    }
    bool inline is_first_iteration(){
        return !iteration; //if iteration==0 this is true
    }

    spin_t inline mod2Pi(double s){
        return mcrg_utilities::mod2Pi(s); 
    }


    inline int index_o_neighbour (const int& i, std::vector<int> const& dxdydz) const { 
        switch(lattice_type){
            case SquareLatticeIsh:
                return mcrg_utilities::index_o_neighbour_square(i,dxdydz[0],dxdydz[1],L);
            case CubicLatticeIsh:
                return mcrg_utilities::index_o_neighbour_cubic(i,dxdydz[0],dxdydz[1],dxdydz[2],L);
        }
    }

    //This generally calculates all the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    double S_alpha (const std::vector<spin_t>& spins, std::vector<shift_t> shifts) {
        double S_a=0;
        for(int i =0;i<N;++i){
            double tmp=1.;
            for(int j = 0; j <shifts.size();++j){ //Pairwise build up the scalar products
                shift_t s = shifts[j];
                double angle=spins[index_o_neighbour(i,std::get<1>(s))];
                int component=std::get<0>(s);
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

    //this function reduces the spin lattice to a quarter of the size by adding the vectors of 4 neighbouring sites together to one and then calculates the angle back and saves it in OUT
    std::vector<spin_t> reduce(const std::vector<spin_t>& spins, const int& entry_point){ //TODO formulate as switch
        if(reduction_type==ReductionTechnique::Decimation){
            return mcrg_utilities::decimate(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::DecimationCubic){
            return mcrg_utilities::decimate_cubic(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::BlockspinCubic){
            return mcrg_utilities::blockspin_cubic(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::Blockspin){
            return mcrg_utilities::blockspin(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::FerroBlockspin){
            if(is_first_iteration()){
                std::vector<spin_t> working_data=mcrg_utilities::ferromagnetic_transformation(spins,L,entry_point);
                return mcrg_utilities::blockspin(working_data,L,entry_point);
            }
            else
                return mcrg_utilities::blockspin(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::Blockspin4x4){
            return mcrg_utilities::blockspin4x4(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::IsingMajority){
            return mcrg_utilities::ising_majority(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::IsingTieBreaker){
            return mcrg_utilities::ising_tie_breaker(spins,L,entry_point);
        }
        if(reduction_type==ReductionTechnique::FerroBlockspin4x4){
            if(is_first_iteration()){
                std::vector<spin_t> working_data=mcrg_utilities::ferromagnetic_transformation(spins,L,0);
                return mcrg_utilities::blockspin4x4(working_data,L,entry_point);
            }
            else
                return mcrg_utilities::blockspin4x4(spins,L,entry_point);
        }
    }
};
#endif //MCPP_MCRG_H_
