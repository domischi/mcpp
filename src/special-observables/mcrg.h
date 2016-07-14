#include <alps/lattice.h>
#include <cassert>
#include <valarray>
#include <alps/parameter.h>
#include <utility> //pair

class mcrg;

/*template<class G = alps::graph_helper<>::graph_type>*/
class mcrg /*: public alps::graph_helper<G>*/{
public:
    typedef double spin_type;
    mcrg(const alps::Parameters& p, int Iteration_, int MCRG_It_) :
    //graph_helper<G>(p),
    iteration(Iteration_),
    max_iterations(MCRG_It_),
    L(p["L"]),
    N(L*L)
    {
        interactions=interactions_e;
        interactions.insert(interactions.end(),interactions_o.begin(),interactions_o.end());
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
    //second index: symmetric parts (for example two diag terms, here)
    //third index: all pairs which are in a specific S_alpha (for example in a four spin interaction)
    /*static const*/ std::vector<std::vector<std::vector<std::pair<int,int>>>> interactions;
    static const std::vector<std::vector<std::vector<std::pair<int,int>>>> interactions_e;
    static const std::vector<std::vector<std::vector<std::pair<int,int>>>> interactions_o;
    static inline int n_interactions_odd(){ 
        return interactions_o.size();
    }
    static inline int n_interactions_even(){ 
        return interactions_e.size();
    }

    std::pair<std::valarray<double>,std::valarray<double>> measure(const std::vector<spin_type>& spins, alps::ObservableSet& obs){
        std::vector<double> OUT_vec_even;
        std::vector<double> OUT_vec_odd;
        //measure the S_alpha in the not reduced lattice
        for(int i=0;i<interactions.size();++i) {
            if(interactions[i][0].size()%2){//odd
                std::vector<double> tmp = {0.,0.};
                for(int j=0;j<interactions[i].size();++j){
                    auto tmp2=S_alpha_odd(spins,interactions[i][j]);
                    tmp[0]+=tmp2[0];
                    tmp[1]+=tmp2[1];
                }
                tmp[0]/=interactions[i].size();
                tmp[1]/=interactions[i].size();
                OUT_vec_odd.push_back(tmp[0]);
                OUT_vec_odd.push_back(tmp[1]);
            }
            else{
                double tmp=0; //add S_alpha, consists of more than one term in case of diag (as there 2 are the same and both should be counted)
                for(int j=0;j<interactions[i].size();++j){
                    tmp+=S_alpha_even(spins,interactions[i][j]);   
                }
                OUT_vec_even.push_back(tmp/interactions[i].size());
            }
        }
        std::valarray<double> OUT_even(OUT_vec_even.data(),OUT_vec_even.size());
        std::valarray<double> OUT_odd(OUT_vec_odd.data(),OUT_vec_odd.size());
        if(OUT_even.size()/n_interactions_even()*2!=OUT_odd.size()/n_interactions_odd()){
            std::cerr<<"SOMETHING WENT WRONG IN THE MEASUREMENT OF THE MCRG AS THE ODD OUT VECTOR HASN'T TWICE THE SIZE OF THE EVEN ONE, ABORTING...";
            throw 1;
        }
        obs["MCRG S_alpha"+ std::to_string(iteration-1)+" even"]<<OUT_even;
        obs["MCRG S_alpha"+ std::to_string(iteration-1)+" odd"]<<OUT_even;
        if(!is_last_iteration()){//This is not the last instantation of the mcrg
            //reduce the lattice and call the descendant on the renormalized system
            std::vector<spin_type> reduced_spins = reduce(spins);
            std::valarray<double> IN_even(OUT_even.size());
            std::valarray<double> IN_odd(OUT_odd.size());
            std::tie(IN_even, IN_odd) =  descendant->measure(reduced_spins,obs);
            //calculate <S_alpha S_beta> - <S_alpha> <S_beta>
            //and save it into obs
            //even
            std::valarray<double> outout_even(IN_even.size()*IN_even.size()),outin_even(IN_even.size()*IN_even.size());
            for(int i=0;i<IN_even.size();++i){
                for(int j=0;j<IN_even.size();++j){
                    outout_even[i*IN_even.size()+j]=OUT_even[i]*OUT_even[j];
                    outin_even[i*IN_even.size()+j]=OUT_even[i]*IN_even[j];
                }
            }
            obs["MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration-1)+" even"]<<outout_even;
            obs["MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration)+" even"  ]<<outout_even;
            //odd
            std::valarray<double> outout_odd(IN_odd.size()*IN_odd.size()/4),outin_odd(IN_odd.size()*IN_odd.size()/4);
            for(int i=0;i<IN_odd.size();i+=2){
                for(int j=0;j<IN_odd.size();j+=2){ 
                    outout_odd[(i*IN_odd.size()/2+j)/2]=OUT_odd[i]*OUT_odd[j]+OUT_odd[i+1]*OUT_odd[j+1];
                    outin_odd[(i*IN_odd.size()/2+j)/2]=OUT_odd[i]*IN_odd[j]+OUT_odd[i+1]*IN_odd[j+1];
                }
            }
            obs["MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration-1)+" odd"]<<outout_odd;
            obs["MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration)+" odd"]<<outin_odd;
        }
        return std::make_pair(OUT_even,OUT_odd); 
    }

    void init_observables(alps::ObservableSet& obs){
        obs << alps::RealVectorObservable("MCRG S_alpha"+ std::to_string(iteration-1)+" even");
        obs << alps::RealVectorObservable("MCRG S_alpha"+ std::to_string(iteration-1)+" odd");
        if(!is_last_iteration()){
            obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration-1)+" even");
            obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration)+" even");
            obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration-1)+" odd");
            obs << alps::RealVectorObservable("MCRG S_alpha"+std::to_string(iteration-1) +" S_beta"+std::to_string(iteration)+" odd");

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
        return max_iterations<=iteration;
    }

    //return the index of a site j which is i+dx+dy in the periodic lattice
    inline int index_o_neighbour (int i, int dx,int dy){
         int x,y;
         i+=N; dx+=N; dy+=N;
         x=(i%L+dx)%L;
         y=(i/L+dy)%L;
         assert(x<L&&y<L);
         return x+L*y;
    }

    //This generally calculates all the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    //For example, just field term (one spin correlator) shifts.size()==0
    //For example, NN along x: shifts.size()==1, shifts[0]=(1,0)
    //For example, 4 spin interaction in a vortex type : shifts.size()==3 shifts[0]==(1,0), shifts[1]==(1,1), shifts[2]==(0,1), 
    double S_alpha_even (const std::vector<spin_type>& spins, std::vector<std::pair<int,int>> shifts) {
        double S_a=0;
        for(int i =0;i<N;++i){
            double tmp=1;
            for(int j = 0; j <shifts.size();j+=2){
                auto& s = shifts[j];
                auto& t = shifts[j+1];
                tmp*=std::cos(spins[index_o_neighbour(i,s.first,s.second)])*std::cos(spins[index_o_neighbour(i,t.first,t.second)])+std::sin(spins[index_o_neighbour(i,s.first,s.second)])*std::sin(spins[index_o_neighbour(i,t.first,t.second)]);
            }
            S_a+=tmp;
        }
        return S_a/N;
    }
    std::vector<double> S_alpha_odd (const std::vector<spin_type>& spins, std::vector<std::pair<int,int>> shifts) {//TODO this might be wrong
        std::vector<double> S_a={0.,0.};
        for(int i =0;i<N;++i){
            double tmp=1; 
            for(int i=1;i<shifts.size();i+=2){
                auto& s = shifts[i];
                auto& t = shifts[i+1];
                tmp*=std::cos(spins[index_o_neighbour(i,s.first,s.second)])*std::cos(spins[index_o_neighbour(i,t.first,t.second)])+std::sin(spins[index_o_neighbour(i,s.first,s.second)])*std::sin(spins[index_o_neighbour(i,t.first,t.second)]);
            }
            S_a[0]+=tmp*std::cos(spins[i]);
            S_a[1]+=tmp*std::sin(spins[i]);
        }
        S_a[0]/=N;
        S_a[1]/=N;
        return S_a;
    }

    void init_blocks() {
        std::vector<int> in1,in2,in3,in4;
        for(int i=0;i<N;i+=2) { //should definitly check this....
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
            std::vector<int> tmp;
            tmp.push_back(in1[i]);
            tmp.push_back(in2[i]);
            tmp.push_back(in3[i]);
            tmp.push_back(in4[i]);
            blocks.push_back(tmp);
        }
    }

    //this function reduces the spin lattice to a quarter of the size by adding the vectors of 4 neighbouring sites together to one and then calculates the angle back and saves it in OUT
    std::vector<spin_type> reduce(const std::vector<spin_type>& spins){
        std::vector<spin_type> OUT;
        for(int i=0;i<blocks.size();++i) {
            double c=0;
            double s=0;
            for(int j : blocks[i]) {
                c+=std::cos(spins[j]); 
                s+=std::sin(spins[j]); 
            }
            OUT.push_back(std::atan2(s,c));
        }
        return OUT;
    }
};

const std::vector<std::vector<std::vector<std::pair<int,int>>>> mcrg::interactions_o = {
        //ODD
        {{std::make_pair(0,0)}},
        {{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(0,1)}},//3 part interaction for simple nn
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(0,1)}},//3 part interaction for simple nn
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(1,0)}},//3 part interaction for simple nn
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(-1,0)}},//3 part interaction for simple nn
        //{{std::make_pair(0,1),std::make_pair(1,0),std::make_pair(1,1)}},//3 part interaction for simple nn//this makes problems
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,2)}},//3 part interaction for simple nn 
        {{std::make_pair(0,0),std::make_pair(2,0),std::make_pair(1,0)}},//3 part interaction for simple nn
        {{std::make_pair(0,0),std::make_pair(0,2),std::make_pair(0,1)}},//3 part interaction for simple nn
        {{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(2,1)}}, 
        {{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(2,-1)}},
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,1)}},
        {{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(-2,1)}},
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(1,2)}},
        {{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(-1,2)}},
        {{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(1,2)}},
        {{std::make_pair(0,0),std::make_pair(1,0),std::make_pair(-1,2)}},
        {{std::make_pair(0,0),std::make_pair(2,1),std::make_pair(1,2)}},
        {{std::make_pair(0,0),std::make_pair(-2,1),std::make_pair(-1,2)}},
        {{std::make_pair(0,0),std::make_pair(-2,-1),std::make_pair(-1,-2)}},
        {{std::make_pair(0,0),std::make_pair(2,-1),std::make_pair(1,-2)}},
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(0,2)}}, 
        {{std::make_pair(0,0),std::make_pair(-1,1),std::make_pair(0,2)}},
        {{std::make_pair(0,0),std::make_pair(1,1),std::make_pair(2,0)}}
    };
const std::vector<std::vector<std::vector<std::pair<int,int>>>> mcrg::interactions_e = {
        //EVEN
        {{std::make_pair(0,0),std::make_pair(1,0)}},//NNx
        {{std::make_pair(0,0),std::make_pair(0,1)}},//NNy
        {{std::make_pair(0,0),std::make_pair(1,1)},{std::make_pair(0,0),std::make_pair(-1,1)}}, //diagonal
        {{std::make_pair(0,0),std::make_pair(2,0)},{std::make_pair(0,0),std::make_pair(0,2)}}, //NNN
        {{std::make_pair(0,0),std::make_pair(3,0)}}, //NNNNx
        {{std::make_pair(0,0),std::make_pair(0,3)}}, //NNNNy
        {{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(1,0),std::make_pair(1,1)}}, // 4 spin
        {{std::make_pair(0,0),std::make_pair(1,2)}},
        {{std::make_pair(0,0),std::make_pair(2,1)}},
        {{std::make_pair(0,0),std::make_pair(-2,1)}},
        {{std::make_pair(0,0),std::make_pair(-1,2)}},
        {{std::make_pair(0,0),std::make_pair(-2,2)}},
        {{std::make_pair(0,0),std::make_pair(4,0)}},
        {{std::make_pair(0,0),std::make_pair(0,4)}},
        {{std::make_pair(0,0),std::make_pair(-1,3)}},
        {{std::make_pair(0,0),std::make_pair(1,3)}},
        {{std::make_pair(0,0),std::make_pair(3,-1)}},
        {{std::make_pair(0,0),std::make_pair(3,1)}}
    }; 

//const std::vector<std::vector<std::vector<std::pair<int,int>>>> mcrg::interactions=std::copy(interactions_e.begin(),interactions_e.end(),interactions.end()).insert(interactions.end(),interactions_o.begin(),interactions_o.end());


