#ifndef MCPP_MCRG_UTILITIES_H_
#define MCPP_MCRG_UTILITIES_H_

#include <vector>
#include <cmath>
namespace mcrg_utilities{
    typedef double spin_t;

    spin_t inline mod2Pi(double s){
        return s-std::floor(s/(2*M_PI))*2*M_PI;
    }

   	inline int index_o_neighbour_square (int i, int dx,int dy, int L){
        int x,y;
		int N=L*L;
        i+=N; dx+=N; dy+=N;
        x=(i/L+dx)%L;
        y=(i%L+dy)%L;
        assert(x<L&&y<L);
        return x+L*y;
    }
   	inline int index_o_neighbour_cubic (int i, int dx,int dy, int dz, int L){
        int x,y,z;
		int N=L*L*L;
        i+=N; dx+=N; dy+=N; dz+=N;
        x=(i/(L*L)+dx)%L;
        y=(i/L+dy)%L;
        z=(i%L+dz)%L;
        assert(x<L&&y<L&&z<L);
        return x+L*(y+L*z);
    }

}
#endif //MCPP_MCRG_UTILITIES_H_
