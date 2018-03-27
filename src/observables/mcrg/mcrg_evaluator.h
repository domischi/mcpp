#ifndef MCPP_MCRG_EVALUATOR_H_
#define MCPP_MCRG_EVALUATOR_H_

// HACK: has to be used to provide set_zero for std::valarray
// This is in the mrcg_evaluator class, to extract the jackknife information
// The following lines have to be added into alps/src/alps/utility/set_zero.hpp
// Otherwise the following wont compile.
//#include <iterator>
//namespace alps{
//  template <class X>
//  inline typename boost::enable_if<is_sequence<std::valarray<X> >,void>::type
//  set_zero(std::valarray<X> & x)
//  {
//    for(int i=0;i<x.size();++i) x[i]=X();
//	//std::fill(std::begin(x),std::end(x),typename element_type<X>::type());
//  }
//}

#include <valarray>
#include <alps/parameter.h>
#include <iostream>
#include "utilities.h"
#ifndef NARMADILLO
#include <armadillo>
#endif //NARMADILLO

class mcrg_evaluator {
public:
    mcrg_evaluator(alps::Parameters const& params)
    #ifndef NARMADILLO
        :
        mcrg_it_depth(params["mcrg_iteration_depth"]),
        error_msg_already_printed(false),
        d(alps::graph_helper<>(params).dimension()),
        b(mcpp::mcrg::init_b(mcpp::mcrg::init_reduction_technique(params)))
    #else //no Armadillo
    #endif// NARMADILLO
    {
        #ifndef NARMADILLO
        #else //no Armadillo
            std::cerr<<"mc++ was built without Armadillo, analysis not possible"<<std::endl;
        #endif// NARMADILLO
    }
    void evaluate(alps::ObservableSet& obs) const {
        #ifndef NARMADILLO
        for(int iteration=1;iteration<=mcrg_it_depth;++iteration){
            for(auto sector : {std::string("e"), std::string("o")}){
                std::string start="MCRG"+sector+" ";
                if(  obs.has(start+"S_alpha"+std::to_string(iteration))
                   &&obs.has(start+"S_alpha"+std::to_string(iteration-1))
                   &&obs.has(start+"S_alpha"+std::to_string(iteration)+" S_beta"+std::to_string(iteration))
                   &&obs.has(start+"S_alpha"+std::to_string(iteration-1)+" S_beta"+std::to_string(iteration))){
                    const alps::RealVectorObsevaluator ObsEval_Sa   =obs[start+"S_alpha"+std::to_string(iteration)];
                    const alps::RealVectorObsevaluator ObsEval_Sap  =obs[start+"S_alpha"+std::to_string(iteration-1)];
                    const alps::RealVectorObsevaluator ObsEval_SaSb =obs[start+"S_alpha"+std::to_string(iteration)+" S_beta"+std::to_string(iteration)];
                    const alps::RealVectorObsevaluator ObsEval_SaSbp=obs[start+"S_alpha"+std::to_string(iteration-1)+" S_beta"+std::to_string(iteration)];
                    #ifdef ALPS_HAS_VALARRAY_SET_ZERO
                        evaluate_w_errors(obs , "lambda"+sector+std::to_string(iteration), ObsEval_Sa   ,ObsEval_Sap  ,ObsEval_SaSb ,ObsEval_SaSbp);
                    #else  // ALPS_HAS_VALARRAY _SET_ZERO
                        evaluate_wo_errors(obs, "lambda"+sector+std::to_string(iteration), ObsEval_Sa   ,ObsEval_Sap  ,ObsEval_SaSb ,ObsEval_SaSbp);
                    #endif // ALPS_HAS_VALARRAY_SET_ZERO
                }
            }
            alps::RealObsevaluator l_e = obs["lambdae"+std::to_string(iteration)];
            alps::RealObsevaluator l_o = obs["lambdao"+std::to_string(iteration)];
            if(l_e.count()>0 && l_o.bin_number()==l_e.bin_number()){
                alps::RealObsevaluator alpha("MCRG alpha"+std::to_string(iteration));
                alps::RealObsevaluator beta ("MCRG beta" +std::to_string(iteration));
                alps::RealObsevaluator gamma("MCRG gamma"+std::to_string(iteration));
                alps::RealObsevaluator delta("MCRG delta"+std::to_string(iteration));
                alps::RealObsevaluator eta  ("MCRG eta"  +std::to_string(iteration));
                alps::RealObsevaluator nu   ("MCRG nu"   +std::to_string(iteration));
                alpha=-d*log(b)/log(l_e)+2;
                beta =(d*log(b)-log(l_o))/log(l_e);
                gamma=(-d*log(b)+2*log(l_o))/log(l_e);
                delta=log(l_o)/(d*log(b)-log(l_o));
                eta  =d+2.-2*log(l_o)/log(b);
                nu   =log(b)/log(l_e);
                obs.addObservable(alpha);
                obs.addObservable(beta );
                obs.addObservable(gamma);
                obs.addObservable(delta);
                obs.addObservable(eta  );
                obs.addObservable(nu   );
            }
        }
        #else //no Armadillo
        #endif// NARMADILLO
    }
private:
    #ifndef NARMADILLO
    const int mcrg_it_depth;
    const int b;
    const int d;
    mutable bool error_msg_already_printed;
    std::complex<double> get_EV(std::valarray<double> const& Data_Sa, std::valarray<double> const& Data_Sap, std::valarray<double> const& Data_SaSb, std::valarray<double> const& Data_SaSbp) const {
        arma::mat T=mcpp::mcrg::get_T(Data_Sa, Data_Sap, Data_SaSb, Data_SaSbp);
        std::complex<double> ev=arma::max(arma::eig_gen(T));
        if(std::abs(ev.real())<10*std::abs(ev.imag())){
            std::cerr << "Strange EV for MCRG: "<<ev<<std::endl;
        }
        return ev;
    }
    #ifdef ALPS_HAS_VALARRAY_SET_ZERO
    //Note: The error gets estimated around a factor 100 times smaller
    //This is probably due as the jackknife data uses most of the information
    //whereas a real observable would just use one bin.
    void evaluate_w_errors( alps::ObservableSet& obs, const std::string& obs_name,
                    const alps::RealVectorObsevaluator& ObsEval_Sa   ,
                    const alps::RealVectorObsevaluator& ObsEval_Sap  ,
                    const alps::RealVectorObsevaluator& ObsEval_SaSb ,
                    const alps::RealVectorObsevaluator& ObsEval_SaSbp) const {
        alps::RealObservable l_r(obs_name);
        const std::vector<std::valarray<double>> J_Sa   =alps::alea::mcdata<std::valarray<double>>(ObsEval_Sa   ).jackknife();
        const std::vector<std::valarray<double>> J_Sap  =alps::alea::mcdata<std::valarray<double>>(ObsEval_Sap  ).jackknife();
        const std::vector<std::valarray<double>> J_SaSb =alps::alea::mcdata<std::valarray<double>>(ObsEval_SaSb ).jackknife();
        const std::vector<std::valarray<double>> J_SaSbp=alps::alea::mcdata<std::valarray<double>>(ObsEval_SaSbp).jackknife();
        const int n_j=J_Sa.size();
        assert(n_bins=SaSb.bin_number() && n_bins=Sap.bin_number()&&n_bins=SaSbp.bin_number());
        for(int del=1;del<n_j;++del) { //delete 1 analysis, starting from 1 as 0th is meanvalue
            std::valarray<double> Data_Sa   =J_Sa   [del];
            std::valarray<double> Data_Sap  =J_Sap  [del];
            std::valarray<double> Data_SaSb =J_SaSb [del];
            std::valarray<double> Data_SaSbp=J_SaSbp[del];
            std::complex<double> ev=get_EV(Data_Sa,Data_Sap, Data_SaSb, Data_SaSbp);
            if(ev.imag()>.1*ev.real()){
                std::cerr<<"Strange EV encountered in MCRG: "<<ev<<std::endl;
            }
            else{
                l_r<<ev.real();
            }
        }
        if(l_r.count()>0)
            obs.addObservable(l_r);
    }
    #else  // ALPS_HAS_VALARRAY_SET_ZERO
    #endif // ALPS_HAS_VALARRAY_SET_ZERO
    void evaluate_wo_errors( alps::ObservableSet& obs, const std::string& obs_name,
                    const alps::RealVectorObsevaluator& ObsEval_Sa   ,
                    const alps::RealVectorObsevaluator& ObsEval_Sap  ,
                    const alps::RealVectorObsevaluator& ObsEval_SaSb ,
                    const alps::RealVectorObsevaluator& ObsEval_SaSbp) const {
        alps::SimpleRealObservable l_r(obs_name);
        if(!error_msg_already_printed){
            std::cerr<<"WARNING: MCRG EVs are only calculated with the mean, not jackknifed\n"
                <<"to support jackknifed data, set ALPS_HAS_VALARRAY_SET_ZERO (and make sure it is that way), then recompile"<<std::endl;
            error_msg_already_printed=true;
        }
        const std::valarray<double> Data_Sa   =ObsEval_Sa   .mean();
        const std::valarray<double> Data_Sap  =ObsEval_Sap  .mean();
        const std::valarray<double> Data_SaSb =ObsEval_SaSb .mean();
        const std::valarray<double> Data_SaSbp=ObsEval_SaSbp.mean();
        auto ev=get_EV(Data_Sa,Data_Sap, Data_SaSb, Data_SaSbp);
        if(ev.imag()>.1*ev.real()){
            std::cerr<<"Strange EV encountered in MCRG: "<<ev<<std::endl;
        }
        else{
            l_r<<ev.real();
            obs.addObservable(l_r);
        }
    }
    #else //no Armadillo
    #endif// NARMADILLO
};
#endif//MCPP_MCRG_EVALUATOR_H_
