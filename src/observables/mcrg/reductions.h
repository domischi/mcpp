#ifndef MCPP_MCRG_REDUCTIONS_H_
#define MCPP_MCRG_REDUCTIONS_H_

#include <vector>
#include <cmath>

#include "utilities.h"

namespace mcrg_utilities{
    typedef double spin_t;

	std::vector<spin_t> decimate(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2)
                ret.push_back(spins[index_o_neighbour_square(EntryPoint,dx,dy,L)]);
        ret.shrink_to_fit();
        return ret; 
	}

    std::vector<spin_t> ising_majority(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint) {
        std::vector<spin_t> ret;
        for(int dx=0; dx<L;dx+=3)
            for(int dy=0; dy<L;dy+=3){
                double angle_sum=
                    mod2Pi(spins[index_o_neighbour_square(EntryPoint,  dx,  dy,L)]) 
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,1+dx,  dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,2+dx,  dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,  dx,1+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,1+dx,1+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,2+dx,1+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,  dx,2+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,1+dx,2+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,2+dx,2+dy,L)]);
                if(angle_sum>4.5*M_PI)
                    ret.push_back(M_PI);
                else
                    ret.push_back(0.);
            }
        ret.shrink_to_fit();
        return ret; 
    }

	std::vector<spin_t> blockspin_cubic(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2)
                for(int dz=0; dz<L;dz+=2){
                    double s000=spins[index_o_neighbour_cubic(EntryPoint,  dx,  dy,  dz,L)]; 
                    double s100=spins[index_o_neighbour_cubic(EntryPoint,1+dx,  dy,  dz,L)];
                    double s010=spins[index_o_neighbour_cubic(EntryPoint,  dx,1+dy,  dz,L)];
                    double s110=spins[index_o_neighbour_cubic(EntryPoint,1+dx,1+dy,  dz,L)];
                    double s001=spins[index_o_neighbour_cubic(EntryPoint,  dx,  dy,1+dz,L)]; 
                    double s101=spins[index_o_neighbour_cubic(EntryPoint,1+dx,  dy,1+dz,L)];
                    double s011=spins[index_o_neighbour_cubic(EntryPoint,  dx,1+dy,1+dz,L)];
                    double s111=spins[index_o_neighbour_cubic(EntryPoint,1+dx,1+dy,1+dz,L)];
                    double c=std::cos(s000)+std::cos(s100)+std::cos(s010)+std::cos(s110)+std::cos(s001)+std::cos(s101)+std::cos(s011)+std::cos(s111);
                    double s=std::sin(s000)+std::sin(s100)+std::sin(s010)+std::sin(s110)+std::sin(s001)+std::sin(s101)+std::sin(s011)+std::sin(s111);
                    ret.push_back(std::atan2(s,c));
                }
        ret.shrink_to_fit();
        return ret; 
	}
	std::vector<spin_t> blockspin(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2){
                double s00=spins[index_o_neighbour_square(EntryPoint,  dx,  dy,L)]; 
                double s10=spins[index_o_neighbour_square(EntryPoint,1+dx,  dy,L)];
                double s01=spins[index_o_neighbour_square(EntryPoint,  dx,1+dy,L)];
                double s11=spins[index_o_neighbour_square(EntryPoint,1+dx,1+dy,L)];
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
                    double spin=spins[index_o_neighbour_square(EntryPoint,i+dx,j+dy,L)];
                    c+=std::cos(spin);
                    s+=std::sin(spin);
                }
                ret.push_back(std::atan2(s,c));
            }
        ret.shrink_to_fit();
        return ret; 
	}
 
    std::vector<spin_t> ferromagnetic_transformation(const std::vector<spin_t>& spins, int L, int EntryPoint){
        std::vector<spin_t> ret=spins;
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2){
                int index;
                //LowerLeft ok
                
                //UpperLeft
                index=index_o_neighbour_square(EntryPoint,  dx,1+dy,L);
                ret[index]=mod2Pi(M_PI-spins[index]);
                
                //LowerRight
                index=index_o_neighbour_square(EntryPoint,1+dx,dy,L);
                ret[index]=mod2Pi(-spins[index]);
                
                //UpperRight
                index=index_o_neighbour_square(EntryPoint,1+dx,1+dy,L);
                ret[index]=mod2Pi(M_PI+spins[index]);
            }
        return ret;
    }
}
#endif //MCPP_MCRG_REDUCTIONS_H_
