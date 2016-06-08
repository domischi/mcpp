#include <alps/graph_helper.h>
#include <cassert>

template<typename G = graph_helper<>::graph_type>
class mcrg : public graph_helper<G>{
public:
    typedef double spin_type;
    mcrg(const alps::Parameters& p, int Iteration_) :
    gh(p),
    iteration(Iteration_),
    L(p["L"]),
    N(L*L)
    {
        //TODO change parameters L in a copy...
        alps::Parameter p_descendant=p;
        p_descendant["L"]=L/2;

        if(iteration>0) {
            descendance=mcrg(p_descendant,Iteration_-1);    
        }
    }

    void measure(const std::vector<spin_type>& spins){
        
        if(iteration>0){
            
        } 
    }

private:
    mcrg descendant;
    int iteration;
    
    int L,N;

    std::vector<std::vector<int>> blocks;

    double NN(const std::vector<spin_type>& spins) {//S_alpha=(S_iS_i+y + S_iS_i+x)/2
        double S_alpha=0;
        int N=0;
        for(site_iterator siter=sites().first; siter!=sites().second;++siter){
            for (neighbour_iterator niter=neighbors(*siter).first; niter!=neighbors(*siter).second; ++niter){
                if(*siter<*niter) { //only check once
                    S_alpha+=std::cos(spins[*siter])*std::cos(spins[*niter])
                            +std::sin(spins[*siter])*std::sin(spins[*niter]);
                }
            }
        }
        return S_alpha/N;
    }
    double diagonal(const std::vector<spin_type>& spins) {//S_alpha=S_iS_i+x+y
        double S_alpha=0;
        int N=0;
        for(site_iterator siter=sites().first; siter!=sites().second;++siter){
            std::vector<int> candidates, diagnonal;
            for (neighbour_iterator niter=neighbors(*siter).first; niter!=neighbors(*siter).second; ++niter)
                for(neighbour_iterator niter2=neighbors(*niter).first; niter2!=neighbors(*niter).second; ++niter2) candidates.push_back(*niter);
                std::sort(candidates.first(),candidates.last());
                std::vector<double>::iterator it=candidates.begin();
                while(it!=candidates.last()){
                    it=std::adjacency_find(it,candidates.end());
                    diagonal.push_back(it);
                    ++(++it);
                }
                for(int i : diagonal)
                if(*siter<i) { //only check once
                    S_alpha+=std::cos(spins[*siter])*std::cos(spins[i])
                            +std::sin(spins[*siter])*std::sin(spins[i]);
                }
            
        }
        return S_alpha/N;
    }
    double NNN(const std::vector<spin_type>& spins) {//S_alpha=(S_iS_i+2x + S_iS_i+2y)/2 
        double S_alpha=0;
        int N=0;
        for(site_iterator siter=sites().first; siter!=sites().second;++siter){
            std::vector<int> candidates, nnn;
            for (neighbour_iterator niter=neighbors(*siter).first; niter!=neighbors(*siter).second; ++niter)
                for(neighbour_iterator niter2=neighbors(*niter).first; niter2!=neighbors(*niter).second; ++niter2) candidates.push_back(*niter);
                std::sort(candidates.first(),candidates.last());
                std::vector<double>::iterator it=candidates.begin();
                while(it!=candidates.last()){
                    it=std::adjacency_find(it,candidates.end());
                    candidates[it]=-1;
                    candidates[it+1]=-1;
                    ++(++it);
                }
                
                for(int i : candidates)
                if(i>=0&&*siter<i) { //only check once
                    S_alpha+=std::cos(spins[*siter])*std::cos(spins[i])
                            +std::sin(spins[*siter])*std::sin(spins[i]);
                }
            
        }
        return S_alpha/N;
    }

    void init_blocks() {
        assert(lattice_name_=="square lattice"); //not sure if this works
        assert(!(L%2));
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
        for(int i=0;i<L/4;++i){
            std::vector<int> tmp;
            tmp.push_back(in1[i]);
            tmp.push_back(in2[i]);
            tmp.push_back(in3[i]);
            tmp.push_back(in4[i]);
            blocks.push_back(tmp);
        }
        for(int i=0;i<blocks.size():++i){ //for debugging my code
            std::cout<<"Block "<< i << ":";
            for (int j : blocks[i]) std::cout<<" "<<j:
            std::cout<<std::endl;   
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
