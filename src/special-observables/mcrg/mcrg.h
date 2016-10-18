#ifndef MCPP_MCRG_H_
#define MCPP_MCRG_H_

#include <alps/lattice.h>
#include <cassert>
#include <valarray>
#include <alps/parameter.h>
#include <tuple>
#include "interactions.h"
#include "reductions.h"
class mcrg; //Needed as it is recursively called
class mcrg {
public:
    typedef double spin_t;
    typedef mcrg_utilities::shift_t shift_t; //(dx,dy, component)
    enum ReductionTechnique {Decimation,Blockspin, FerroBlockspin, Blockspin4x4, FerroBlockspin4x4};
    mcrg(const alps::Parameters& p, int Iteration_, int MCRG_It_) :
    iteration(Iteration_),
    max_iterations(MCRG_It_),
    L(p["L"]),
    N(L*L)
    {
        //Choose Interaction Set
        if(p["MCRG Interactions"]=="small"){
            interactions=mcrg_utilities::small;
        } else if(p["MCRG Interactions"]=="medium interactions"){
            interactions=mcrg_utilities::medium1;
        } else if(p["MCRG Interactions"]=="medium range" || p["MCRG Interactions"]=="medium" ){
            interactions=mcrg_utilities::medium2;
        } else if(p["MCRG Interactions"]=="massive"){
            std::cout << "Consider this to be a way to large interaction set. This is just me having fun implementing interactions, however if you wanna use it, feel free to do so, however it will take forever and a day to finish";
            interactions=mcrg_utilities::massive;
        } else {
            std::cerr<<"Did not recognise the interaction set to use for MCRG, aborting..."<<std::endl;
            std::exit(20);
        }
        for(auto& i : interactions) i.shrink_to_fit();
        interactions.shrink_to_fit();
        
        //Choose reduction technique
        if(p["MCRG Reduction Technique"]=="Decimation"){
            reduction_type=ReductionTechnique::Decimation;
        } else if(p["MCRG Reduction Technique"]=="Blockspin"){
            reduction_type=ReductionTechnique::Blockspin;
        } else if(p["MCRG Reduction Technique"]=="FerroBlockspin"){
            reduction_type=ReductionTechnique::FerroBlockspin;
        } else if(p["MCRG Reduction Technique"]=="Blockspin4x4"){
            reduction_type=ReductionTechnique::Blockspin4x4;
        } else if(p["MCRG Reduction Technique"]=="FerroBlockspin4x4"){
            reduction_type=ReductionTechnique::FerroBlockspin4x4;
        } else {
            std::cerr << "Did not recognise the reduction process to use for MCRG, aborting..."<<std::endl;
            std::exit(21); 
        }
        assert(lattice_name_=="square lattice"); //not sure if this works
        assert(!(L%2));
       
        //If this is not the last iteration, construct the MCRG class for the next smaller lattice 
        if(!is_last_iteration()) {
            alps::Parameters p_descendant=p;
            p_descendant["L"]=L/2;
            descendant=std::make_shared<mcrg>(p_descendant,iteration+1,MCRG_It_);
        }
    }
                
    //first index: which interaction
    //second index: listing the vectors of the shift
    //third index (shift_t): a pair of ints denoting the shift wrt to the first one (always (0,0))
    std::vector<std::vector<shift_t>> interactions;
    inline int n_interactions(){
        return interactions.size();
    }


    std::tuple<std::valarray<double>,std::valarray<double>> measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs){
        //std::valarray<double> OUT(n_interactions());
        std::vector<double> oute, outo;
        int counto=0;
        int counte=0;
        //measure the S_alpha in the not reduced lattice
        for(int i=0;i<interactions.size();++i) {
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
        save_outo/=N;
        std::valarray<double> save_oute=OUTe;
        save_oute/=N;
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
        outoute/=N;
        outouto/=N;
        obs["MCRGo S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outouto;
        obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outoute;
        
        if(!is_last_iteration()){//This is not the last instantation of the mcrg
            //reduce the lattice and call the descendant on the renormalized system
            std::vector<spin_t> reduced_spins = reduce(spins);
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
            outine/=N;
            obs["MCRGe S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outine;
            std::valarray<double> outino(INo.size()*INo.size());
            for(int i=0;i<INo.size();++i){
                for(int j=0;j<INo.size();++j){
                    outino[i*INo.size()+j]=(OUTo[i]*INo[j]);
                }
            }
            outino/=N;
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

private:
    std::shared_ptr<mcrg> descendant;
    int iteration;
    int max_iterations; 
    int L,N;

    ReductionTechnique reduction_type;

    bool inline is_last_iteration(){
        return max_iterations<iteration;
    }
    bool inline is_first_iteration(){
        return !iteration; //if iteration==0 this is true
    }

    spin_t inline mod2Pi(double s){
        return mcrg_utilities::mod2Pi(s); 
    }

    inline int index_o_neighbour (int i, int dx,int dy){ return mcrg_utilities::index_o_neighbour(i,dx,dy,L);}

    //This generally calculates all the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    //Example: field term along x: shift[j].size()==1, {{0,0,0}}
    //Example: NNx along xy: shift[j].size()==2, { {0,0,0} , {1,0,1} }
    //For more examples have a look at the interactions.h file
    double S_alpha (const std::vector<spin_t>& spins, std::vector<shift_t> shifts) {
        double S_a=0;
        for(int i =0;i<N;++i){
            double tmp=1.;
            for(int j = 0; j <shifts.size();++j){ //Pairwise build up the scalar products
                shift_t s = shifts[j];
                double angle=spins[index_o_neighbour(i,std::get<0>(s),std::get<1>(s))];
                int component=std::get<2>(s);
                if(component==0){
                    tmp*=std::cos(angle);
                } else if(component==1){
                    tmp*=std::sin(angle);
                } else {
                    std::cout << "In S_alpha with a not recognised component, The values is "<< component<<"... Aborting..."<< std::endl;
                    std::exit(19);
                }
            }
            S_a+=tmp;
        }
        return S_a;
    }

    //this function reduces the spin lattice to a quarter of the size by adding the vectors of 4 neighbouring sites together to one and then calculates the angle back and saves it in OUT
    std::vector<spin_t> reduce(const std::vector<spin_t>& spins){
        if(reduction_type==ReductionTechnique::Decimation){
            return mcrg_utilities::decimate(spins,L,0);//TODO make the entry point random
        }
        if(reduction_type==ReductionTechnique::Blockspin){
            return mcrg_utilities::blockspin(spins,L,0);//TODO make the entry point random
        }
        if(reduction_type==ReductionTechnique::FerroBlockspin){
            if(is_first_iteration()){
                std::vector<spin_t> working_data=mcrg_utilities::ferromagnetic_transformation(spins,L,0); //TODO make entry point random
                return mcrg_utilities::blockspin(working_data,L,0);
            }
            else
                return mcrg_utilities::blockspin(spins,L,0);//TODO make the entry point random
        }
        if(reduction_type==ReductionTechnique::Blockspin4x4){
            return mcrg_utilities::blockspin4x4(spins,L,0);//TODO make the entry point random
        }
        if(reduction_type==ReductionTechnique::FerroBlockspin4x4){
            if(is_first_iteration()){
                std::vector<spin_t> working_data=mcrg_utilities::ferromagnetic_transformation(spins,L,0); //TODO make entry point random
                return mcrg_utilities::blockspin4x4(working_data,L,0);
            }
            else
                return mcrg_utilities::blockspin4x4(spins,L,0);//TODO make the entry point random
        }
    }
};
#endif //MCPP_MCRG_H_
