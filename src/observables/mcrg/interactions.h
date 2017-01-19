#ifndef MCPP_MCRG_INTERACTIONS_H_
#define MCPP_MCRG_INTERACTIONS_H_

#include <vector>
#include <tuple>

//TODO now bigger sets copy the content of the smaller ones, this is not desirable as you should never do copy paste...

namespace mcrg_utilities {
typedef std::tuple<int,std::vector<int>> shift_t; //component, (dx,dy,dz,...)

/*
    This is only a small set, including NNx, NNy, d11, d-11
    In the components x, y, xy, yx
    Odd:  only field term
    Even: only 2 particle interactions
*/
std::vector<std::vector<shift_t>> small = {
    //Odd
    {std::make_tuple(0,std::vector<int> {0,0} )}, // field term along x
    {std::make_tuple(1,std::vector<int> {0,0} )}, // field term along y
    //Even
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 0} )}, //NNx x comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 0} )}, //NNx y comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 0} )}, //NNx xy comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, //NNx yx comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, //NNy x comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 0, 1} )}, //NNy y comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 0, 1} )}, //NNy xy comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, //NNy yx comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, //d11 x comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}, //d11 y comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}, //d11 xy comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, //d11 yx comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 1} )}, //d-11 x comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {-1, 1} )}, //d-11 y comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {-1, 1} )}, //d-11 xy comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 1} )}  //d-11 yx comp
};


std::vector<std::vector<shift_t>> dIsing = {
    //Odd
        //field term
        {std::make_tuple(0,std::vector<int> {0,0} )},
        //3 body interaction
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,1} ), std::make_tuple(0,std::vector<int> {0,2} )}, // (00) (01) (02) 
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,1} ), std::make_tuple(0,std::vector<int> {1,0} )}, // (00) (01) (10) 
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,1} ), std::make_tuple(0,std::vector<int> {2,0} )}, // (00) (01) (20) 
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,2} ), std::make_tuple(0,std::vector<int> {2,0} )}, // (00) (02) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {1,1} ), std::make_tuple(0,std::vector<int> {2,2} )}, // (00) (11) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {1,1} ), std::make_tuple(0,std::vector<int> {2,1} )}, // (00) (11) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {1,1} ), std::make_tuple(0,std::vector<int> {2,0} )}, // (00) (11) (20)
    //Even
        //same component
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 0} )}, // (00)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 3} )}, // (33)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (50)
        //quartic terms
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // plaquette
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // tetris tile
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}  // tetris L tile
};
// Do a hand crafted, good set for dXY //TODO understand a good handcrafted interaction set
std::vector<std::vector<shift_t>> dXY_handcrafted = {
    //Odd
        //field term
        {std::make_tuple(0,std::vector<int> {0,0} )},
        //3 body interaction
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} )}, // cubic site term
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {1,0} )}, // (00) (00) (10) 
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {1,1} )}, // (00) (00) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {1,1} )}, // (00) (00) (11) cross
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {2,0} )}, // (00) (00) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {2,1} )}, // (00) (00) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {2,1} )}, // (00) (00) (21) cross
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {2,2} )}, // (00) (00) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {2,2} )}, // (00) (00) (22) cross
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {3,0} )}, // (00) (00) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {3,1} )}, // (00) (00) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {3,1} )}, // (00) (00) (31) cross
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {3,2} )}, // (00) (00) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {3,2} )}, // (00) (00) (32) cross
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {3,3} )}, // (00) (00) (33)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> {3,3} )}, // (00) (00) (33) cross
    //Even
        //same component
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 0} )}, // (00)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 3} )}, // (33)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (50)
        //cross component
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}, // (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 2, 1} )}, // (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 2, 2} )}, // (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 3, 1} )}, // (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 3, 2} )}, // (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 3, 3} )}, // (33)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 4, 1} )}, // (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 4, 2} )}, // (42)
        //Biquadratic terms
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 0, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (11) same component
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 0, 0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}  // (11) cross component
};
std::vector<std::vector<shift_t>> very_small = { 
    //Odd
    {std::make_tuple(0,std::vector<int> {0,0} )}, // field term along x
    //Even
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 0} )}, //NNx x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 0} )}, //NNx xy comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, //NNy x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, //d11 x comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}, //d11 y comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(1,std::vector<int> { 1, 1} )}, //d11 xy comp
    {std::make_tuple(1,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}  //d11 yx comp
};
/*
    This is a small set designed for Ising
    Only in the component x -> as the simulation is running for Ising
    Odd:  only field term
    Even: NNx, NNy, d11, d-11 and 4 Spin
*/
std::vector<std::vector<shift_t>> small_Ising = {
    //Odd
    {std::make_tuple(0,std::vector<int> {0,0} )}, // field term along x
    //Even
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, //NNx x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, //NNy x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, //d11 x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 1} )}, //d-11 x comp
    {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )} // 4 spin
};
//ATTENTION: This interaction was autogenerated with the following parameters for the function
//(3,[1, 3, 1], RESTRICTED_4=True)
//It consists of 137 different interactions.
std::vector<std::vector<shift_t>> medium_range = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 3, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 3}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 3, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 3}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-3}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-2, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-2, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-2, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-2,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-2,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-3, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 3, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 2, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 2, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 2, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 2,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 2,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1,-2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 3}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 2}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) }
};
//ATTENTION: This interaction was autogenerated with the following parameters for the function
//(4,[1, 1, 1, 1, 1], RESTRICTED_4=False)
//It consists of 105 different interactions.
std::vector<std::vector<shift_t>> medium_interaction = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 0,-1}), std::make_tuple(1,std::vector<int>{-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0}), std::make_tuple(1,std::vector<int>{ 0, 1}), std::make_tuple(1,std::vector<int>{ 1, 1}) }
};

//ATTENTION: This interaction was autogenerated with the following parameters for the function
//(2,[1, 2]) and num_components=3
//It consists of 17 different interactions.
std::vector<std::vector<shift_t>> xy_3d_small = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(1,std::vector<int>{ 0, 0, 0}) }
};
std::vector<std::vector<shift_t>> xy_3d_very_small = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(1,std::vector<int>{ 0, 0, 0}), std::make_tuple(1,std::vector<int>{ 1, 0, 0}) },
};
std::vector<std::vector<shift_t>> xy_3d_medium = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 2}) }
};

//ATTENTION: This interaction was autogenerated with the following parameters for the function
//(4,[1, 3, 1, 1]) and num_components=3 , restricted to x direction
//It consists of 86 different interactions.
std::vector<std::vector<shift_t>> xy_3d_massive = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 3, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2,-2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 2,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1,-2,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 3, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 2,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1,-2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 3}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 1, 0}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 1}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0,-1}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0,-1, 0}), std::make_tuple(0,std::vector<int>{-1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0}), std::make_tuple(0,std::vector<int>{ 0, 1}), std::make_tuple(0,std::vector<int>{ 1, 1}) }
};
// This is the interaction set also Swendsen used in his publication 1983
std::vector<std::vector<shift_t>> xy_3d_swendsen = 
{
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 1, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2, 1}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 2, 2, 2}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 0, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 0}) },
	{ std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 0, 0, 0}), std::make_tuple(0,std::vector<int>{ 1, 1, 1}) }
};
}
#endif //MCPP_MCRG_INTERACTIONS_H_
