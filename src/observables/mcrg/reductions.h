#ifndef MCPP_MCRG_REDUCTIONS_H_
#define MCPP_MCRG_REDUCTIONS_H_

#include <vector>
#include <cmath>

#include "utilities.h"
namespace mcpp {
namespace mcrg{
    typedef double spin_t;
    enum ReductionTechnique {Decimation,Blockspin, FerroBlockspin, Blockspin4x4, FerroBlockspin4x4, IsingMajority, IsingTieBreaker, BlockspinCubic, DecimationCubic};

    ReductionTechnique init_reduction_technique(std::string const& s, LatticeType const lattice_type) {
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
    ReductionTechnique init_reduction_technique(alps::Parameters const& p){
        return init_reduction_technique(static_cast<std::string>(p["MCRG Reduction Technique"]), init_lattice_type(p));
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
    int init_b(alps::Parameters const& p) {
        return init_b(init_reduction_technique(p));
    }
	std::vector<spin_t> decimate(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2)
                ret.push_back(spins[index_o_neighbour_square(EntryPoint,dx,dy,L)]);
        ret.shrink_to_fit();
        return ret;
	}
	std::vector<spin_t> decimate_cubic(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint){
        std::vector<spin_t> ret; //TODO do a bit of index magic, to not need to use the vector capabilities, such that it get's fastened up a bit
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2)
                for(int dz=0; dz<L;dz+=2)
                    ret.push_back(spins[index_o_neighbour_cubic(EntryPoint,dx,dy,dz,L)]);
        ret.shrink_to_fit();
        return ret;
	}

    std::vector<spin_t> ising_tie_breaker(const std::vector<spin_t>& spins, const int& L, const int& EntryPoint) {
        std::vector<spin_t> ret;
        for(int dx=0; dx<L;dx+=2)
            for(int dy=0; dy<L;dy+=2){
                double angle_sum=
                    mod2Pi(spins[index_o_neighbour_square(EntryPoint,  dx,  dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,1+dx,  dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,  dx,1+dy,L)])
                   +mod2Pi(spins[index_o_neighbour_square(EntryPoint,1+dx,1+dy,L)]);
                if(angle_sum>5./2.*M_PI)
                    ret.push_back(M_PI);
                else if(angle_sum<3./2.*M_PI)
                    ret.push_back(0.);
                else
                    ret.push_back(spins[index_o_neighbour_square(EntryPoint,  dx,  dy,L)]);
            }
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
        std::vector<spin_t> ret;
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

    //this function reduces the spin lattice to a quarter of the size by adding
    //the vectors of 4 neighbouring sites together to one and then calculates
    //the angle back and saves it in OUT
    std::vector<spin_t> reduce(const std::vector<spin_t>& spins, const int& entry_point, const int& L, ReductionTechnique const& reduction_type, bool const& is_first_iteration){
        switch(reduction_type){
            case ReductionTechnique::Decimation:
                return mcpp::mcrg::decimate(spins,L,entry_point);
            case ReductionTechnique::DecimationCubic:
                return mcpp::mcrg::decimate_cubic(spins,L,entry_point);
            case ReductionTechnique::BlockspinCubic:
                return mcpp::mcrg::blockspin_cubic(spins,L,entry_point);
            case ReductionTechnique::Blockspin:
                    return mcpp::mcrg::blockspin(spins,L,entry_point);
            case ReductionTechnique::FerroBlockspin:
                if(is_first_iteration){
                    std::vector<spin_t> working_data=mcpp::mcrg::ferromagnetic_transformation(spins,L,entry_point);
                    return mcpp::mcrg::blockspin(working_data,L,entry_point);
                }
                else
                    return mcpp::mcrg::blockspin(spins,L,entry_point);
            case ReductionTechnique::Blockspin4x4:
                return mcpp::mcrg::blockspin4x4(spins,L,entry_point);
            case ReductionTechnique::IsingMajority:
                return mcpp::mcrg::ising_majority(spins,L,entry_point);
            case ReductionTechnique::IsingTieBreaker:
                return mcpp::mcrg::ising_tie_breaker(spins,L,entry_point);
            case ReductionTechnique::FerroBlockspin4x4:
                if(is_first_iteration){
                    std::vector<spin_t> working_data=mcpp::mcrg::ferromagnetic_transformation(spins,L,0);
                    return mcpp::mcrg::blockspin4x4(working_data,L,entry_point);
                }
                else
                    return mcpp::mcrg::blockspin4x4(spins,L,entry_point);
        }
    }
}//mcrg
}//mcpp
#endif //MCPP_MCRG_REDUCTIONS_H_
