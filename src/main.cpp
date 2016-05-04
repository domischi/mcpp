//#include <alps/osiris/comm.h>

#include "single-spin-flip/ssf.h"

#include <iostream>

int main(int argc,char** argv){
    //if (argc!=2){
    //    std::cerr<< "Usage: mc++ input_file.xml"<<std::endl;
    //    return 1;
    //}
    //else{
        return alps::scheduler::start(argc,argv, ssf_factory());
    //}
}
