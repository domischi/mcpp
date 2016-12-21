#ifndef MCPP_UTILITIES_H_
#define MCPP_UTILITIES_H_

#include <alps/lattice.h>

namespace mcpp{
    int init_N (const alps::Parameters& p){
        return alps::graph_helper<>(p).num_sites();
    }

}
#endif //MCPP_UTILITIES_H_
