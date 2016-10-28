#ifndef MCPP_MCRG_INTERACTIONS_H_
#define MCPP_MCRG_INTERACTIONS_H_

#include <vector>
#include <tuple>

//TODO now bigger sets copy the content of the smaller ones, this is not desirable as you should never do copy paste...

namespace mcrg_utilities {
typedef std::tuple<int,int,int> shift_t;

/*
    This is only a small set, including NNx, NNy, d11, d-11
    In the components x, y, xy, yx
    Odd:  only field term
    Even: only 2 particle interactions
*/
std::vector<std::vector<shift_t>> small = {
		//Odd
		{std::make_tuple(0,0,0)}, // field term along x
        {std::make_tuple(0,0,1)}, // field term along y
		//Even
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
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,0)} //d-11 yx comp
    };
/*
    This is a small set designed for Ising
    Only in the component x -> as the simulation is running for Ising
    Odd:  only field term
    Even: NNx, NNy, d11, d-11 and 4 Spin
*/
//TODO check if One restricts itself to only non-equivalent sites (for Ising e.g. NNx and NNy are equivalent) the strange degeneracy gets lifted, this might also be the case for the dipolar interaction as essentially the sites are equivalent
std::vector<std::vector<shift_t>> small_Ising = {
		//Odd
		{std::make_tuple(0,0,0)}, // field term along x
		//Even
		{std::make_tuple(0,0,0),std::make_tuple(1,0,0)}, //NNx 
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0)}, //NNy 
        {std::make_tuple(0,0,0),std::make_tuple(1,1,0)}, //d11 
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,0)},//d-11 
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0),std::make_tuple(1,0,0),std::make_tuple(1,1,0)}, // 4 Spin
    };

/*
    This is medium set, including NNx, NNy, d11, d-11
    In the components x, y, xy, yx
    Odd:  field term and some 3 particle interactions
    Even: 4 particle interaction for square with all components, all 2 particle interactions
	Remark: only very near neighbour interactions are implemented, but more complicated ones
*/
//TODO this is broken.....      there is at least one interaction that is more than once in, as simulations break with this...
std::vector<std::vector<shift_t>> medium1 = {
		//Odd
		{std::make_tuple(0,0,0)}, // field term along x
        {std::make_tuple(0,0,1)}, // field term along y
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,0), std::make_tuple(0,-1,0)}, // ‾| xxx 
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,0), std::make_tuple(0,-1,1)}, // ‾| xxy
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,1), std::make_tuple(0,-1,0)}, // ‾| xyx
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,1), std::make_tuple(0,-1,1)}, // ‾| xyy
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,0), std::make_tuple(0,-1,0)}, // ‾| yxx
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,0), std::make_tuple(0,-1,1)}, // ‾| yxy
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,1), std::make_tuple(0,-1,0)}, // ‾| yyx
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,1), std::make_tuple(0,-1,1)}, // ‾| yyy
        {std::make_tuple(0,0,0), std::make_tuple(1,0,0),  std::make_tuple(0,-1,0)}, // |‾ xxx 
        {std::make_tuple(0,0,0), std::make_tuple(1,0,0),  std::make_tuple(0,-1,1)}, // |‾ xxy
        {std::make_tuple(0,0,0), std::make_tuple(1,0,1),  std::make_tuple(0,-1,0)}, // |‾ xyx
        {std::make_tuple(0,0,0), std::make_tuple(1,0,1),  std::make_tuple(0,-1,1)}, // |‾ xyy
        {std::make_tuple(0,0,1), std::make_tuple(1,0,0),  std::make_tuple(0,-1,0)}, // |‾ yxx
        {std::make_tuple(0,0,1), std::make_tuple(1,0,0),  std::make_tuple(0,-1,1)}, // |‾ yxy
        {std::make_tuple(0,0,1), std::make_tuple(1,0,1),  std::make_tuple(0,-1,0)}, // |‾ yyx
        {std::make_tuple(0,0,1), std::make_tuple(1,0,1),  std::make_tuple(0,-1,1)}, // |‾ yyy
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,0), std::make_tuple(0,1,0)} , // _| xxx 
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,0), std::make_tuple(0,1,1)} , // _| xxy
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,1), std::make_tuple(0,1,0)} , // _| xyx
        {std::make_tuple(0,0,0), std::make_tuple(-1,0,1), std::make_tuple(0,1,1)} , // _| xyy
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,0), std::make_tuple(0,1,0)} , // _| yxx
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,0), std::make_tuple(0,1,1)} , // _| yxy
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,1), std::make_tuple(0,1,0)} , // _| yyx
        {std::make_tuple(0,0,1), std::make_tuple(-1,0,1), std::make_tuple(0,1,1)} , // _| yyy
        {std::make_tuple(0,0,0), std::make_tuple(1,0,0),  std::make_tuple(0,1,0)} , // |_ xxx 
        {std::make_tuple(0,0,0), std::make_tuple(1,0,0),  std::make_tuple(0,1,1)} , // |_ xxy
        {std::make_tuple(0,0,0), std::make_tuple(1,0,1),  std::make_tuple(0,1,0)} , // |_ xyx
        {std::make_tuple(0,0,0), std::make_tuple(1,0,1),  std::make_tuple(0,1,1)} , // |_ xyy
        {std::make_tuple(0,0,1), std::make_tuple(1,0,0),  std::make_tuple(0,1,0)} , // |_ yxx
        {std::make_tuple(0,0,1), std::make_tuple(1,0,0),  std::make_tuple(0,1,1)} , // |_ yxy
        {std::make_tuple(0,0,1), std::make_tuple(1,0,1),  std::make_tuple(0,1,0)} , // |_ yyx
        {std::make_tuple(0,0,1), std::make_tuple(1,0,1),  std::make_tuple(0,1,1)} , // |_ yyy
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,0),  std::make_tuple(0,1,0)} , // 2 ^ 1 xxx
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,0),  std::make_tuple(0,1,1)} , // 2 ^ 1 xxy
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,1),  std::make_tuple(0,1,0)} , // 2 ^ 1 xyx=yxx
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,1),  std::make_tuple(0,1,1)} , // 2 ^ 1 xyy=yxy
        //{std::make_tuple(0,0,1), std::make_tuple(0,0,1),  std::make_tuple(0,1,0)} , // 2 ^ 1 yyx
        //{std::make_tuple(0,0,1), std::make_tuple(0,0,1),  std::make_tuple(0,1,1)} , // 2 ^ 1 yyy
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,0),  std::make_tuple(1,0,0)} , // 2-1 xxx 
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,0),  std::make_tuple(1,0,1)} , // 2-1 xxy
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,1),  std::make_tuple(1,0,0)} , // 2-1 xyx=yxx
        //{std::make_tuple(0,0,0), std::make_tuple(0,0,1),  std::make_tuple(1,0,1)} , // 2-1 xyy=yxy
        //{std::make_tuple(0,0,1), std::make_tuple(0,0,1),  std::make_tuple(1,0,0)} , // 2-1 yyx
        //{std::make_tuple(0,0,1), std::make_tuple(0,0,1),  std::make_tuple(1,0,1)} , // 2-1 yyy
		//Even
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
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0),std::make_tuple(1,0,0),std::make_tuple(1,1,0)}, // 4 spin xxxx 
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0),std::make_tuple(1,0,0),std::make_tuple(1,1,1)}, // 4 spin xxxy
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0),std::make_tuple(1,0,1),std::make_tuple(1,1,0)}, // 4 spin xxyx
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0),std::make_tuple(1,0,1),std::make_tuple(1,1,1)}, // 4 spin xxyy
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1),std::make_tuple(1,0,0),std::make_tuple(1,1,0)}, // 4 spin xyxx
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1),std::make_tuple(1,0,0),std::make_tuple(1,1,1)}, // 4 spin xyxy
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1),std::make_tuple(1,0,1),std::make_tuple(1,1,0)}, // 4 spin xyyx
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1),std::make_tuple(1,0,1),std::make_tuple(1,1,1)}, // 4 spin xyyy
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0),std::make_tuple(1,0,0),std::make_tuple(1,1,0)}, // 4 spin yxxx
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0),std::make_tuple(1,0,0),std::make_tuple(1,1,1)}, // 4 spin yxxy
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0),std::make_tuple(1,0,1),std::make_tuple(1,1,0)}, // 4 spin yxyx
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0),std::make_tuple(1,0,1),std::make_tuple(1,1,1)}, // 4 spin yxyy
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1),std::make_tuple(1,0,0),std::make_tuple(1,1,0)}, // 4 spin yyxx
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1),std::make_tuple(1,0,0),std::make_tuple(1,1,1)}, // 4 spin yyxy
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1),std::make_tuple(1,0,1),std::make_tuple(1,1,0)}, // 4 spin yyyx
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1),std::make_tuple(1,0,1),std::make_tuple(1,1,1)}  // 4 spin yyyy
    
    };
/*
    This is medium set, including NNx, NNy, d11, d-11, NNNx, NNNy, d-21, d-12, d21, d12
    In the components x, y, xy, yx
    Odd:  field term 
    Even: only 2 particle interactions
	Remark: here a medium set is introduced via a longer range interaction, instead of more complicated interactions
*/
std::vector<std::vector<shift_t>> medium2 = {
		//Odd
		{std::make_tuple(0,0,0)}, // field term along x
        {std::make_tuple(0,0,1)}, // field term along y
		//Even
		{std::make_tuple(0,0,0),std::make_tuple(1,0,0)}, //NNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,1)}, //NNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,0,1)}, //NNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,0)}, //NNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0)}, //NNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1)}, //NNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1)}, //NNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0)}, //NNy yx comp
		{std::make_tuple(0,0,0),std::make_tuple(2,0,0)}, //NNNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,1)}, //NNNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(2,0,1)}, //NNNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,0)}, //NNNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,0)}, //NNNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,1)}, //NNNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,1)}, //NNNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,0)}, //NNNy yx comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,0)}, //d11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,1)}, //d11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,1)}, //d11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,0)}, //d11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,0)},//d-11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,1)},//d-11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,1)},//d-11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,0)},//d-11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(2,1,0)}, //d21 x comp
        {std::make_tuple(0,0,1),std::make_tuple(2,1,1)}, //d21 y comp
        {std::make_tuple(0,0,0),std::make_tuple(2,1,1)}, //d21 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(2,1,0)}, //d21 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-2,1,0)},//d-21 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-2,1,1)},//d-21 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-2,1,1)},//d-21 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-2,1,0)},//d-21 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(1,2,0)}, //d12 x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,2,1)}, //d12 y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,2,1)}, //d12 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,2,0)}, //d12 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,2,0)},//d-12 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,2,1)},//d-12 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,2,1)},//d-12 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,2,0)} //d-12 yx comp
    };

/*
    This is massive, gigantic, totally useless, most exaggerated set, including NNx, NNy, d11, d-11, NNNx, NNNy, d-21, d-12, d21, d12
    In the components x, y, xy, yx
    Odd:  field term 
    Even: only 2 particle interactions
	Remark: here a medium set is introduced via a longer range interaction, instead of more complicated interactions
*/
std::vector<std::vector<shift_t>> massive = {
		//Odd
		{std::make_tuple(0,0,0)}, // field term along x
        {std::make_tuple(0,0,1)}, // field term along y
		//Even
		{std::make_tuple(0,0,0),std::make_tuple(1,0,0)}, //NNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,1)}, //NNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,0,1)}, //NNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,0,0)}, //NNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,0)}, //NNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,1)}, //NNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,1,1)}, //NNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,1,0)}, //NNy yx comp
		{std::make_tuple(0,0,0),std::make_tuple(2,0,0)}, //NNNx x comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,1)}, //NNNx y comp
        {std::make_tuple(0,0,0),std::make_tuple(2,0,1)}, //NNNx xy comp
        {std::make_tuple(0,0,1),std::make_tuple(2,0,0)}, //NNNx yx comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,0)}, //NNNy x comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,1)}, //NNNy y comp
        {std::make_tuple(0,0,0),std::make_tuple(0,2,1)}, //NNNy xy comp
        {std::make_tuple(0,0,1),std::make_tuple(0,2,0)}, //NNNy yx comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,0)}, //d11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,1)}, //d11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,1,1)}, //d11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,1,0)}, //d11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,0)},//d-11 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,1)},//d-11 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,1,1)},//d-11 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,1,0)},//d-11 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(2,1,0)}, //d21 x comp
        {std::make_tuple(0,0,1),std::make_tuple(2,1,1)}, //d21 y comp
        {std::make_tuple(0,0,0),std::make_tuple(2,1,1)}, //d21 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(2,1,0)}, //d21 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-2,1,0)},//d-21 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-2,1,1)},//d-21 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-2,1,1)},//d-21 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-2,1,0)},//d-21 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(1,2,0)}, //d12 x comp
        {std::make_tuple(0,0,1),std::make_tuple(1,2,1)}, //d12 y comp
        {std::make_tuple(0,0,0),std::make_tuple(1,2,1)}, //d12 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(1,2,0)}, //d12 yx comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,2,0)},//d-12 x comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,2,1)},//d-12 y comp
        {std::make_tuple(0,0,0),std::make_tuple(-1,2,1)},//d-12 xy comp
        {std::make_tuple(0,0,1),std::make_tuple(-1,2,0)} //d-12 yx comp
    };
}
#endif //MCPP_MCRG_INTERACTIONS_H_
