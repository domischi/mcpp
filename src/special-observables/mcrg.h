#include <alps/lattice.h>
#include <cassert>
#include <valarray>
#include <alps/parameter.h>
#include <utility> //pair
#include <tuple>

class mcrg;
/*template<class G = alps::graph_helper<>::graph_type>*/
class mcrg /*: public alps::graph_helper<G>*/{
public:
    typedef double spin_t;
    typedef std::tuple<int,int,int> shift_t; //(dx,dy, component)
    mcrg(const alps::Parameters& p, int Iteration_, int MCRG_It_) :
    iteration(Iteration_),
    max_iterations(MCRG_It_),
    L(p["L"]),
    N(L*L)
    {
        interactions=interactions_e;
        interactions.insert(interactions.end(),interactions_o.begin(),interactions_o.end());
        for(auto& i : interactions) i.shrink_to_fit();
        interactions.shrink_to_fit();
        assert(lattice_name_=="square lattice"); //not sure if this works
        assert(!(L%2));
        alps::Parameters p_descendant=p;
        p_descendant["L"]=L/2;
        init_blocks();
        if(!is_last_iteration()) {
            descendant=std::make_shared<mcrg>(p_descendant,iteration+1,MCRG_It_);
        }
    }
                
    //first index: which interaction
    //second index: listing the vectors of the shift
    //third index (shift_t): a pair of ints denoting the shift wrt to the first one (always (0,0))
    /*static const*/ std::vector<std::vector<shift_t>> interactions;
    static const std::vector<std::vector<shift_t>> interactions_e;
    static const std::vector<std::vector<shift_t>> interactions_o;
    static inline int n_interactions_odd(){ 
        return interactions_o.size();
    }
    static inline int n_interactions_even(){ 
        return interactions_e.size();
    }
    static inline int n_interactions(){
        return n_interactions_odd()+n_interactions_even();
    }


    std::valarray<double> measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs){
        std::valarray<double> OUT(n_interactions());
        std::vector<spin_t> working_data;
        //if(is_first_iteration()){
        //    working_data=ferromagnetic_transformation(spins);
        //}
        //else
        working_data=spins;
        int count_o=0;
        int count_e=0;
        //std::vector<double> OUT_vec_even;
        //std::vector<double> OUT_vec_odd;
        //measure the S_alpha in the not reduced lattice
        for(int i=0;i<interactions.size();++i) {
            //tmp=S_alpha(working_data,interactions[i]);
            OUT[i]=S_alpha(working_data,interactions[i]);
        }
        obs["MCRG S_alpha"+ std::to_string(iteration)]<<OUT;
        //measure <S_alpha n S_beta n>
        std::valarray<double> outout(OUT.size()*OUT.size());
        
        for(int i=0;i<n_interactions();++i){
            for(int j=0;j<n_interactions();++j){
                outout[i*OUT.size()+j]=OUT[i]*OUT[j];
            }
        }
        obs["MCRG S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration)]<<outout;
        
        if(!is_last_iteration()){//This is not the last instantation of the mcrg
            //reduce the lattice and call the descendant on the renormalized system
            std::vector<spin_t> reduced_spins = reduce(working_data);
            std::valarray<double> IN(OUT.size());
            IN =  descendant->measure(reduced_spins,obs);
            //calculate <S_alpha n-1 S_beta n>
            //and save it into obs
            //even
            std::valarray<double> outin(IN.size()*IN.size());
            for(int i=0;i<IN.size();++i){
                for(int j=0;j<IN.size();++j){
                    outin[i*IN.size()+j]=(OUT[i]*IN[j]);
                }
            }
            obs["MCRG S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1)]<<outin;
        }
        return OUT;
    }

    void init_observables(alps::ObservableSet& obs){
        obs << alps::RealVectorObservable("MCRG S_alpha"+ std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration));
        obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration) +" S_beta"+std::to_string(iteration+1));
        if(!is_last_iteration()){
            descendant->init_observables(obs);
        }
    }

private:
    std::shared_ptr<mcrg> descendant;
    int iteration;
    int max_iterations; 
    int L,N;

    std::vector<std::vector<int>> blocks;

    bool inline is_last_iteration(){
        return max_iterations<iteration;
    }
    bool inline is_first_iteration(){
        return !iteration; //if iteration==0 this is true
    }

    //std::vector<spin_t> ferromagnetic_transformation(const std::vector<spin_t>& spins){
    //    std::vector<spin_t> ret=spins;//check if this function exists
    //    for(auto& b: blocks){
    //        for(int i=0;i<b.size();++i){
    //            int site=b[i];
    //            switch(i){
    //                case 0: //left lower 
    //                    break;
    //                case 1: //left upper
    //                    ret[i]=mod2Pi(M_PI-spins[i]); // mirror at y axis (x -> -x)
    //                    break;
    //                case 2: //right lower
    //                    ret[i]=mod2Pi(-spins[i]); // mirror at x axis (y -> -y) 
    //                    break;
    //                case 3: //right upper
    //                    ret[i]=mod2Pi(spins[i]+M_PI); // mirror at origin (x -> -x, y->-y)
    //                    break;
    //                default:
    //                    std::cout<<"Something went wrong in the ferromagn. transfo in MCRG, index encountered was "<<i<<"... Aborting..."<<std::endl;
    //                    std::exit(8);
    //            }
    //        }
    //    }
    //    return ret;
    //}

    spin_t inline mod2Pi(double s){
        return s-std::floor(s/(2*M_PI))*2*M_PI;
    }

    //return the index of a site j which is i+dx+dy in the periodic lattice
    inline int index_o_neighbour (int i, int dx,int dy){
         int x,y;
         i+=N; dx+=N; dy+=N;
         x=(i/L+dx)%L;
         y=(i%L+dy)%L;
         assert(x<L&&y<L);
         return x+L*y;
    }

    //This generally calculates all the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    //For example, just field term (one spin correlator) shifts.size()==0
    //For example, NN along x: shifts.size()==1, shifts[0]=(1,0)
    //For example, 4 spin interaction in a vortex type : shifts.size()==3 shifts[0]==(1,0), shifts[1]==(1,1), shifts[2]==(0,1), 
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
        //return S_a/N;
    }

    //Initialize the block list, such that I know which blocks to be reduced.
    //This generates the vectors of blocks
    //each block contains 4 spins, being in the order of left lower, left upper, right lower, right upper
    void init_blocks() {
        std::vector<int> in1,in2,in3,in4;
        for(int i=0;i<N;i+=2) {
            if((i/L)%2){//odd (L+1,L+3,L+5,...)
                in3.push_back(i);
                in4.push_back(i+1);
            }
            else {//even (1,3,5,...)
                in1.push_back(i);
                in2.push_back(i+1);
            }
        }
        for(int i=0;i<N/4;++i){
            std::vector<int> tmp(4);
            tmp[0]=in1[i];
            tmp[1]=in2[i];
            tmp[2]=in3[i];
            tmp[3]=in4[i];
            blocks.push_back(tmp);
        }
        blocks.shrink_to_fit();
    }

    //this function reduces the spin lattice to a quarter of the size by adding the vectors of 4 neighbouring sites together to one and then calculates the angle back and saves it in OUT
    std::vector<spin_t> reduce(const std::vector<spin_t>& spins){
        std::vector<spin_t> OUT(blocks.size());
        for(int b=0;b<blocks.size();++b) {
            OUT[b]=spins[(blocks[b])[0]]; // decimation
            //double c=0.;
            //double s=0.;
            //for(int j : blocks[b]) {
            //    c+=std::cos(spins[j]); 
            //    s+=std::sin(spins[j]); 
            //}
            //OUT[b]=std::atan2(s,c);
        }
        return OUT;
    }
};

const std::vector<std::vector<mcrg::shift_t>> mcrg::interactions_o = {
        //ODD
        {std::make_tuple(0,0,0)},// field term along x
        {std::make_tuple(0,0,1)} // field term along y
        //{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(0,1)},//3 part interaction for simple nn
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(0,1)},//3 part interaction for simple nn
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(1,0)},//3 part interaction for simple nn
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(-1,0)},//3 part interaction for simple nn
        ////{{std::make_pair(0,1),s>td::make_pair(1,0),std::make_pair(1,1)}},//3 part interaction for simple nn//this makes problems
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,2)},//3 part interaction for simple nn 
        //{std::make_pair(0,0),std::make_pair(2,0),std::make_pair(1,0)},//3 part interaction for simple nn
        //{std::make_pair(0,0),std::make_pair(0,2),std::make_pair(0,1)},//3 part interaction for simple nn
        //{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(2,1)}, 
        //{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(2,-1)},
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,1)},
        //{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(-2,1)},
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(1,2)},
        //{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(-1,2)},
        //{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(1,2)},
        //{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(-1,2)},
        //{std::make_pair(0,0),std::make_pair(2,1),std::make_pair(1,2)},
        //{std::make_pair(0,0),std::make_pair(-2,1),std::make_pair(-1,2)},
        //{std::make_pair(0,0),std::make_pair(-2,-1),std::make_pair(-1,-2)},
        //{std::make_pair(0,0),std::make_pair(2,-1),std::make_pair(1,-2)},
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(0,2)}, 
        //{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(0,2)},
        //{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,0)}
    };
const std::vector<std::vector<mcrg::shift_t>> mcrg::interactions_e = {
        //EVEN
        {std::make_tuple(0,0,0),std::make_tuple(1,0,0)}, //NNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,1)}, //NNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,0,1)}, //NNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,0)}, //NNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0)}, //NNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1)}, //NNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1)}, //NNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0)}, //NNy yx comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,0)}, //d11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,1)}, //d11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,1)}, //d11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,0)}, //d11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,0)},//d-11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,1)},//d-11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,1)},//d-11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,0)},//d-11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(2,0,0)}, //NNNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,1)}, //NNNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(2,0,1)}, //NNNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,0)}, //NNNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,0)}, //NNNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,1)}, //NNNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,1)}, //NNNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,0)}//, //NNNy yx comp
        //{std::make_pair(0,0),std::make_pair(2,0)},
        //{std::make_pair(0,0),std::make_pair(0,2)}, //NNN
        //{std::make_pair(0,0),std::make_pair(3,0)}, //NNNNx
        //{std::make_pair(0,0),std::make_pair(0,3)}, //NNNNy
        //{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(1,0),std::make_pair(1,1)}, // 4 spin
        //{std::make_pair(0,0),std::make_pair(1,2)},
        //{std::make_pair(0,0),std::make_pair(2,1)},
        //{std::make_pair(0,0),std::make_pair(-2,1)},
        //{std::make_pair(0,0),std::make_pair(-1,2)},
        //{std::make_pair(0,0),std::make_pair(-2,2)},
        //{std::make_pair(0,0),std::make_pair(4,0)},
        //{std::make_pair(0,0),std::make_pair(0,4)},
        //{std::make_pair(0,0),std::make_pair(-1,3)},
        //{std::make_pair(0,0),std::make_pair(1,3)},
        //{std::make_pair(0,0),std::make_pair(3,-1)},
        //{std::make_pair(0,0),std::make_pair(3,1)}
    };
