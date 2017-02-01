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
        // // (10)
        // // (20)
        // // (30)
        // // (40)
        // // (50)
        // // (60)
        // // (01)
        // // (11)
        // // (21)
        // // (31)
        // // (41)
        // // (51)
        // // (02)
        // // (12)
        // // (22)
        // // (32)
        // // (42)
        // // (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (01) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (01) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (01) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (01) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (01) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (01) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (01) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (01) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (01) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (01) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (01) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (01) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (01) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (01) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (01) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (01) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (01) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (02) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (02) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (02) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (02) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (02) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (02) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (02) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (02) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (02) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (02) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (02) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (02) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (02) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (02) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (02) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (02) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (02) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (03) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (03) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (03) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (03) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (03) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (03) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (03) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (03) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (03) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (03) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (03) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (03) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (03) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (03) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (03) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (03) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (03) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (10) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (10) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (10) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (10) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (10) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (10) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (10) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (10) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (10) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (10) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (10) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (10) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (10) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (10) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (10) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (10) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (10) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (10) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (11) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (11) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (11) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (11) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (11) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (11) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (11) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (11) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (11) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (11) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (11) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (11) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (11) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (11) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (11) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (11) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (11) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (12) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (12) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (12) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (12) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (12) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (12) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (12) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (12) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (12) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (12) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (12) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (12) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (12) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (12) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (12) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (12) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (12) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (20) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (20) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (20) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (20) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (20) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (20) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (20) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (20) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (20) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (20) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (20) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (20) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (20) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (20) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (20) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (20) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (20) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (21) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (21) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (21) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (21) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (21) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (21) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (21) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (21) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (21) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (21) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (21) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (21) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (21) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (21) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (21) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (21) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (21) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (22) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (22) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (22) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (22) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (22) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (22) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (22) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (22) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (22) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (22) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (22) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (22) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (22) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (22) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (22) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (22) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (22) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (30) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (30) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (30) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (30) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (30) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (30) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (30) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (30) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (30) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (30) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (30) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (30) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (30) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (30) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (30) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (30) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (30) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (31) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (31) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (31) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (31) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (31) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (31) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (31) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (31) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (31) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (31) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (31) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (31) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (31) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (31) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (31) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (31) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (31) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (32) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (32) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (32) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (32) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (32) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (32) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (32) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (32) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (32) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (32) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (32) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (32) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (32) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (32) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (32) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (32) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (32) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (40) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (40) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (40) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (40) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (40) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (40) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (40) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (40) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (40) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (40) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (40) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (40) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (40) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (40) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (40) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (40) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (40) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (41) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (41) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (41) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (41) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (41) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (41) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (41) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (41) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (41) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (41) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (41) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (41) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (41) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (41) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (41) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (41) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (41) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (42) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (42) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (42) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (42) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (42) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (42) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (42) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (42) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (42) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (42) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (42) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (42) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (42) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (42) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (42) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (42) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (42) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (50) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (50) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (50) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (50) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (50) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (50) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (50) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (50) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (50) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (50) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (50) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (50) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (50) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (50) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (50) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (50) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (50) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (51) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (51) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (51) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (51) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (51) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (00) (51) (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (51) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (51) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (51) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (51) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (51) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (51) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (51) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (51) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (51) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (51) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (51) (03)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (00) (60) (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (00) (60) (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (00) (60) (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (00) (60) (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (00) (60) (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (00) (60) (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (00) (60) (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (00) (60) (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (00) (60) (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (00) (60) (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (00) (60) (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (00) (60) (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (00) (60) (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (00) (60) (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (00) (60) (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (00) (60) (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (00) (60) (03)
    //Even
        //same component
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // (10)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}, // (20)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 0} )}, // (30)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 0} )}, // (40)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 0} )}, // (50)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 6, 0} )}, // (60)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // (01)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // (11)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 1} )}, // (21)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 1} )}, // (31)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 1} )}, // (41)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 5, 1} )}, // (51)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 2} )}, // (02)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 1, 2} )}, // (12)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 2, 2} )}, // (22)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 3, 2} )}, // (32)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 4, 2} )}, // (42)
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 3} )}, // (03)
        //quartic terms
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // plaquette
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // tetris tile
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0,-1} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} )}, // tetris tile rot
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> {-1, 0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // tetris lightning tile
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0,-1} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 1, 1} )}, // tetris lightning tile rot
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 0, 2} ), std::make_tuple(0,std::vector<int> { 1, 0} )}, // tetris L tile 
        {std::make_tuple(0,std::vector<int> {0,0} ), std::make_tuple(0,std::vector<int> { 0, 1} ), std::make_tuple(0,std::vector<int> { 1, 0} ), std::make_tuple(0,std::vector<int> { 2, 0} )}  // tetris L tile rot
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
