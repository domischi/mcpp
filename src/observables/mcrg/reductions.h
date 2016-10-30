#ifndef MCPP_MCRG_REDUCTIONS_H_
#define MCPP_MCRG_REDUCTIONS_H_

#include <vector>
#include <cmath>
namespace mcrg_utilities{
    typedef double spin_t;

   	inline int index_o_neighbour (int i, int dx,int dy, int L){
        int x,y;
		int N=L*L;
        i+=N; dx+=N; dy+=N;
        x=(i/L+dx)%L;
        y=(i%L+dy)%L;
        assert(x<L&&y<L);
        return x+L*y;
    }

	std::vector<spin_t> decimate(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2)
                ret.push_back(spins[index_o_neighbour(EntryPoint,dx,dy,L)]);
        ret.shrink_to_fit();
        return ret; 
	}

	std::vector<spin_t> blockspin(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2){
                double s00=spins[index_o_neighbour(EntryPoint,  dx,  dy,L)]; 
                double s10=spins[index_o_neighbour(EntryPoint,1+dx,  dy,L)];
                double s01=spins[index_o_neighbour(EntryPoint,  dx,1+dy,L)];
                double s11=spins[index_o_neighbour(EntryPoint,1+dx,1+dy,L)];
                double c=std::cos(s00)+std::cos(s10)+std::cos(s01)+std::cos(s11);
                double s=std::sin(s00)+std::sin(s10)+std::sin(s01)+std::sin(s11);
                ret.push_back(std::atan2(s,c));
            }
        ret.shrink_to_fit();
        return ret; 
	}
	std::vector<spin_t> blockspin4x4(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=4)
            for(int dy=0; dy<L;dy+=4){
                double c=0.;
                double s=0.;
                for(int i =0; i<4;++i)
                for(int j =0; j<4;++j){
                    double spin=spins[index_o_neighbour(EntryPoint,i+dx,j+dy,L)];
                    c+=std::cos(spin);
                    s+=std::sin(spin);
                }
                ret.push_back(std::atan2(s,c));
            }
        ret.shrink_to_fit();
        return ret; 
	}
 
    spin_t inline mod2Pi(double s){
        return s-std::floor(s/(2*M_PI))*2*M_PI;
    }
    
    std::vector<spin_t> ferromagnetic_transformation(const std::vector<spin_t>& spins, int L, int EntryPoint){
        std::vector<spin_t> ret=spins;
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2){
                int index;
                //LowerLeft ok
                
                //UpperLeft
                index=index_o_neighbour(EntryPoint,  dx,1+dy,L);
                ret[index]=mod2Pi(M_PI-spins[index]);
                
                //LowerRight
                index=index_o_neighbour(EntryPoint,1+dx,dy,L);
                ret[index]=mod2Pi(-spins[index]);
                
                //UpperRight
                index=index_o_neighbour(EntryPoint,1+dx,1+dy,L);
                ret[index]=mod2Pi(M_PI+spins[index]);
            }
        return ret;
    }
}
#endif //MCPP_MCRG_REDUCTIONS_H_
