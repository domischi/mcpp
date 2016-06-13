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
    mcrg(const alps::Parameters& p, int Iteration_) :
    //graph_helper<G>(p),
    iteration(Iteration_),
    L(p["L"]),
    N(L*L)
    {
        assert(lattice_name_=="square lattice"); //not sure if this works
        assert(!(L%2));
        alps::Parameters p_descendant=p;
        p_descendant["L"]=L/2;
        init_blocks();
        if(iteration>0) {
            descendant=std::make_shared<mcrg>(p_descendant,Iteration_-1);
        }
        //init_interactions();
    }
                
    //first index: which interaction
    //second index: symmetric parts (for example two diag terms, here)
    //third index: all pairs which are in a specific S_alpha (for example in a four spin interaction)
    static const std::vector<std::vector<std::vector<std::pair<int,int>>>> interactions;
    static inline int n_interactions(){ //TODO check if this needs first to be initialized or the compiler does this at compile time
        return interactions.size();
    }

    std::valarray<double> measure(const std::vector<spin_type>& spins, alps::ObservableSet& obs){
        std::vector<double> OUT_vec;
        //measure the S_alpha in the not reduced lattice
        for(int i=0;i<interactions.size();++i) {
            double tmp=0; //add S_alpha, consists of more than one term in case of diag (as there 2 are the same and both should be counted)
            for(int j=0;j<interactions[i].size();++j){
                tmp+=S_alpha(spins,interactions[i][j]);   
            }
            OUT_vec.push_back(tmp/interactions[i].size());
        }
        std::valarray<double> OUT(OUT_vec.data(),OUT_vec.size());
        if(iteration>0){//This is not the last instantation of the mcrg
            //reduce the lattice and call the descendant on the renormalized system
            std::vector<spin_type> reduced_spins = reduce(spins);
            std::valarray<double> IN =  descendant->measure(reduced_spins,obs);
            //calculate <S_alpha S_beta> - <S_alpha> <S_beta>
            //and save it into obs
            assert(OUT.size()==IN.size());
            obs["S_alpha" + std::to_string(iteration)]<<OUT;
            std::valarray<double> outout(IN.size()*IN.size()),outin(IN.size()*IN.size());
            for(int i=0;i<IN.size();++i){
                for(int j=0;j<IN.size();++j){
                    outout[i*IN.size()+j]=OUT[i]*OUT[j];
                    outin[i*IN.size()+j]=OUT[i]*IN[j];
                }
            }
            obs["S_alpha S_beta same iteration"+ std::to_string(iteration)]<<outout;
            obs["S_alpha S_beta next iteration"+ std::to_string(iteration)]<<outin;
        }
        return OUT; 
    }

    void init_observables(alps::ObservableSet& obs){
        for(int i =1; i<=iteration;++i){
            obs << alps::RealVectorObservable("S_alpha" + std::to_string(i));
            obs << alps::RealVectorObservable("S_alpha S_beta same iteration" + std::to_string(i));
            obs << alps::RealVectorObservable("S_alpha S_beta next iteration" + std::to_string(i));
        }
    }

private:
    std::shared_ptr<mcrg> descendant;
    int iteration;
     
    int L,N;

    std::vector<std::vector<int>> blocks;

    //return the index of a site j which is i+dx+dy in the periodic lattice
    inline int index_o_neighbour (int i, int dx,int dy){
         int x,y;
         i+=N; dx+=N; dy+=N;
         x=(i%L+dx)%L;
         y=(i/L+dy)%L;
         assert(x<L&&y<L);
         return x+L*y;
    }

    //TODO generalize this to the case if 2 are the same... out of symmetry (i.e. diagonal terms) and do this in the innermost for loop, this should speed up the calculations of this observable dramatically as it reduces the runtime by a factor of the number of different terms 
    //This generally calculates all the S_alpha
    //for this it iterates over all the lattice sites i. For each of them all the shifts are calculated and taken into consideration with the correlation
    //For example, just field term (one spin correlator) shifts.size()==0
    //For example, NN along x: shifts.size()==1, shifts[0]=(1,0)
    //For example, 4 spin interaction in a vortex type : shifts.size()==3 shifts[0]==(1,0), shifts[1]==(1,1), shifts[2]==(0,1), 
    double S_alpha (const std::vector<spin_type>& spins, std::vector<std::pair<int,int>> shifts) {
        double S_a=0;
        for(int i =0;i<N;++i){
            double contribution_of_i_c=1;
            double contribution_of_i_s=1;
            for(auto& s : shifts){
                contribution_of_i_c*=std::cos(spins[index_o_neighbour(i,s.first,s.second)]);
                contribution_of_i_s*=std::sin(spins[index_o_neighbour(i,s.first,s.second)]);
            }
            S_a+=contribution_of_i_c+contribution_of_i_s;
        }
        return S_a/N;
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

const std::vector<std::vector<std::vector<std::pair<int,int>>>> mcrg::interactions = {
        {{std::make_pair(0,0)}},
        {{std::make_pair(0,0),std::make_pair(1,0)}},//NNx
        {{std::make_pair(0,0),std::make_pair(0,1)}},//NNy
        {{std::make_pair(0,0),std::make_pair(1,1)},{std::make_pair(0,0),std::make_pair(-1,1)}}, //diagonal
        {{std::make_pair(0,0),std::make_pair(2,0)}}, //NNNx
        {{std::make_pair(0,0),std::make_pair(0,2)}}, //NNNy
        {{std::make_pair(0,0),std::make_pair(0,1),std::make_pair(1,0),std::make_pair(1,1)}} // 4 spin
    }; 

