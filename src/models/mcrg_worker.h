#ifndef MCPP_MODELS_MCRG_H
#define MCPP_MODELS_MCRG_H
#include <alps/parapack/worker.h>
#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris.h>
#include <alps/osiris/dump.h>
#include <alps/expression.h>
#include <alps/lattice.h>

#include "xy.h"
#include "../observables/mcrg/utilities.h"
template < typename WORKER>
class mcrg_worker : public alps::parapack::mc_worker{
public :
    typedef typename WORKER::spin_t spin_t;
    mcrg_worker(alps::Parameters const& params) :
        alps::parapack::mc_worker(params),
        T(mcpp::init_T(params)),
        L(params["L"]),
        n_even(init_n_even(params)),
        Each_Measurement(params.value_or_default("Each_Measurement",10)),
        mcrg_averages(params.value_or_default("MCRG Temperature Averages",2<<11)),
        Step_Number(0),
        lattice_type(mcpp::mcrg::init_lattice_type(params)),
        reduction_technique(mcpp::mcrg::init_reduction_technique(params)),
        scale_factor_b(mcpp::mcrg::init_b(mcpp::mcrg::init_reduction_technique(params)))
        {
            std::cout
                << alps::logger::header()
                << "Initialize workers for MCRG"<<std::endl;
            alps::Parameters p1=params;
            alps::Parameters p2=params;
            if(!params.defined("MCRG Interactions")){
                p1["MCRG Interactions"]="minimal";
                p2["MCRG Interactions"]="minimal";
            }
            if(static_cast<int>(p2["L"])%scale_factor_b){
                std::cerr<<"L not divisible by the renormalization factor b, aborting..."<<std::endl;
                std::exit(5);
            }
            p2["L"]=static_cast<int>(p2["L"])/scale_factor_b;
            workers.push_back(std::make_shared<WORKER>(p1));
            workers.push_back(std::make_shared<WORKER>(p2));
            mcrg_obs.push_back(std::make_shared<mcrg>(p1, 1, nullptr, nullptr));
            mcrg_obs.push_back(std::make_shared<mcrg>(p2, 1, nullptr, nullptr));
        }
    void init_observables(alps::Parameters const& p, alps::ObservableSet& obs){
        obs<<alps::RealObservable("Tc");
        obs<<alps::RealObservable("Tcp");
        obs<<alps::RealObservable("LHS0");
        obs<<alps::RealObservable("delta_K0");
        workers[0]->init_observables(p,obsset_large);
        workers[1]->init_observables(p,obsset_small);
        mcrg_obs[0]->init_observables(obsset_large);
        mcrg_obs[1]->init_observables(obsset_small);
    }
    void save(alps::ODump &dump) const{
        for(std::shared_ptr<WORKER> const& w:workers)
            w->save(dump);
    }
    void load(alps::IDump &dump){
        for(std::shared_ptr<WORKER> const& w:workers)
            w->load(dump);
    }
    void run(alps::ObservableSet& obs){
        workers[0]->run(obsset_large);
        workers[1]->run(obsset_small);
        if(is_thermalized()){
            if(!(Step_Number%Each_Measurement)){
                //const std::vector<spin_t> spins_big_system   = mcpp::mcrg::reduce(workers[0]->get_spins(),0, L  , mcpp::mcrg::ReductionTechnique::IsingTieBreaker, true);
                //const std::vector<spin_t> spins_small_system = mcpp::mcrg::reduce(workers[1]->get_spins(),0, L/2, mcpp::mcrg::ReductionTechnique::IsingTieBreaker, true);
                const std::vector<spin_t> spins_big_system   = workers[0]->get_spins();
                const std::vector<spin_t> spins_small_system = workers[1]->get_spins();
                mcrg_obs[0]->measure(spins_big_system, obsset_large);
                mcrg_obs[1]->measure(spins_small_system, obsset_small);
                if(mcrg_averages-1==(Step_Number/Each_Measurement)%mcrg_averages){
                    const std::string start="MCRGe ";
                    std::valarray<double> large_Sa1   (alps::RealVectorObsevaluator(obsset_large[start+"S_alpha1"]).mean());
                    std::valarray<double> large_Sa0   (alps::RealVectorObsevaluator(obsset_large[start+"S_alpha0"]).mean());
                    std::valarray<double> large_SaSb00(alps::RealVectorObsevaluator(obsset_large[start+"S_alpha0 S_beta0"]).mean());
                    std::valarray<double> large_SaSb01(alps::RealVectorObsevaluator(obsset_large[start+"S_alpha0 S_beta1"]).mean());
                    std::valarray<double> small_Sa0   (alps::RealVectorObsevaluator(obsset_small[start+"S_alpha0"]).mean());
                    std::valarray<double> small_SaSb00(alps::RealVectorObsevaluator(obsset_small[start+"S_alpha0 S_beta0"]).mean());//they are computed correctly (cross check with mcrg code)
                    large_Sa1   /=std::pow(L/scale_factor_b,2  );
                    large_Sa0   /=std::pow(L               ,2  );
                    large_SaSb00/=std::pow(L               ,2*2);
                    large_SaSb01/=std::pow(L               ,2  )/scale_factor_b;
                    small_Sa0   /=std::pow(L/scale_factor_b,2  );
                    small_SaSb00/=std::pow(L/scale_factor_b,2*2);
                    large_Sa1   /=std::pow(n_even,1);
                    large_Sa0   /=std::pow(n_even,1);
                    large_SaSb00/=std::pow(n_even,2);
                    large_SaSb01/=std::pow(n_even,2);
                    small_Sa0   /=std::pow(n_even,1);
                    small_SaSb00/=std::pow(n_even,2);
                    arma::mat delta_dSdK= mcpp::mcrg::dSdK(large_Sa1, large_Sa0, large_SaSb01)
                                         -mcpp::mcrg::dSdK(small_Sa0, small_Sa0, small_SaSb00);
                    //delta_dSdK/=std::pow(L/scale_factor_b,2*n_even);
                    //delta_dSdK/=std::pow(L/scale_factor_b,2);
                    arma::vec LHS= arma::vec(&large_Sa1[0],large_Sa1.size())
                                  -arma::vec(&small_Sa0[0],small_Sa0.size());
                    arma::vec delta_K=arma::solve(delta_dSdK,LHS);
                    delta_K*=2./scale_factor_b;
                    delta_K/=std::pow(-std::pow(L/scale_factor_b,2), n_even-2);
                    //delta_K*=-1./scale_factor_b;
                    obs["delta_K0"]<<delta_K[0];
                    obs["Tc"] <<1./(1./T-delta_K[0]);
                    obsset_small.reset();
                    obsset_large.reset();
                }
            }
            ++Step_Number;
        }
    }

    bool is_thermalized() const {
        bool s=true;
        for(std::shared_ptr<WORKER> const& w:workers)
            s&=w->is_thermalized();
        return s;
    }
    double progress() const {
        double s=0.;
        for(std::shared_ptr<WORKER> const& w:workers)
            s+=w->progress();
        return s/workers.size();
    }
private:
    double T;
    const int L;
    //int entry_point;
    int Step_Number;
    const int Each_Measurement, mcrg_averages;
    const mcpp::mcrg::LatticeType lattice_type;
    const int scale_factor_b;
    const int n_even;
    //alps::RealObservable delta_K0;
    const mcpp::mcrg::ReductionTechnique reduction_technique;
    std::vector<std::shared_ptr<WORKER>> workers;
    std::vector<std::shared_ptr<mcrg>> mcrg_obs;
    alps::ObservableSet obsset_small, obsset_large;

    int init_n_even(const alps::Parameters p) const{
        const auto interactions=mcpp::mcrg::init_interactions(p);  
        int r=0;
        for (const auto& q : interactions){
            if (q.size()%2==0){
                ++r;
            }
        }
        return r;
    }
    void set_beta(double beta_) {
        T=1/beta_;
        for(auto& w : workers) w->set_beta(beta_);
    }
    //void update_entry_point(){
    //    entry_point=mcpp::mcrg::new_entry_point(entry_point, L,scale_factor_b);
    //}
};
#endif //MCPP_MODELS_MCRG_H
