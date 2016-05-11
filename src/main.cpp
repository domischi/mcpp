//#include <alps/osiris/comm.h>

#include "single-spin-flip/ssf.h"

#include <iostream>

int main(int argc,char** argv){
    //if (argc!=2){
    //    std::cerr<< "Usage: mc++ input_file.xml"<<std::endl;
    //    return 1;
    //}
    //else{
    // TODO this needs to be generalized for any potential factory.... 
        return alps::scheduler::start(argc,argv, ssf_factory());
    //}
}
