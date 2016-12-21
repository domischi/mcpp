#ifndef MCPP_STRUCTURE_FACTOR_H_
#define MCPP_STRUCTURE_FACTOR_H_

#include <alps/lattice.h>
#include <cassert>
#include <algorithm>
#include <valarray>
#include <alps/parameter.h>
#if !NFFTW
    #include <fftw3.h>
#endif
#include "../../utilities.h"
class structure_factor {
public:
    typedef double spin_t;

    structure_factor(const alps::Parameters& p) :
    #if !NFFTW 
    fftw_timelimit_(p.value_or_default("FFTW timelimit for measurement", 10.)),
    #endif //NFFTW
    L(p["L"]),
    N(mcpp::init_N(p))
    {
        alps::graph_helper<> gh(p);
        alps::graph_helper<>::basis_vector_iterator v,v_end;
        for(std::tie(v,v_end)=gh.reciprocal_basis_vectors() ; v!=v_end ; ++v){
            vector_type vec=*v;
            for(auto& e: vec) e/=L; //calculate the Fouriertransform of the super cell, not the elementary one...
            reciprocal_vectors.push_back(vec);
        }
        for(std::tie(v,v_end)=gh.basis_vectors() ; v!=v_end ; ++v){
            vector_type vec=*v;
            basis_vectors.push_back(vec);
        } 
        assert(reciprocal_vectors.size()==basis_vectors.size());
        assert(reciprocal_vectors.size()==2);
        assert(p["LATTICE"]=="square lattice");
        #if !NFFTW
        fftw_in=fftw_alloc_complex(N);
        fftw_out=fftw_alloc_complex(N);
        fftw_set_timelimit(fftw_timelimit_);
        plan=fftw_plan_dft_2d(L,L,fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
        #endif //NFFTW
    }
    // Rule of 3 already spams, but I dont need move semantics, therefore I can leave the move constructor and the move assign operator default 
    structure_factor(const structure_factor &original) :
        L(original.L),
        N(original.N),
        #if !NFFTW 
        fftw_timelimit_(original.fftw_timelimit_),
        #endif //NFFTW
        reciprocal_vectors(original.reciprocal_vectors),
        basis_vectors(original.basis_vectors)
    {
        #if !NFFTW
        fftw_in=fftw_alloc_complex(N);
        fftw_out=fftw_alloc_complex(N);
        fftw_set_timelimit(fftw_timelimit_);
        plan=fftw_plan_dft_2d(L,L,fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE | FFTW_DESTROY_INPUT);
        memcpy(fftw_in , original.fftw_in , N*sizeof(fftw_complex));
        memcpy(fftw_out, original.fftw_out, N*sizeof(fftw_complex));
        #endif //NFFTW
    }

    structure_factor& operator=(const structure_factor& original) = delete; //as there are const members which could be very dangerous, just forbid the use of an assignment

    ~structure_factor(){
        #if !NFFTW
        fftw_free(fftw_in ); 
        fftw_free(fftw_out); 
        fftw_destroy_plan(plan);
        #endif //NFFTW
    }

    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) {
        std::valarray<std::complex<double>> Sx(spins.size());
        std::valarray<std::complex<double>> Sy(spins.size());
        std::transform(spins.begin(),spins.end(),begin(Sx),[](double d) {return std::cos(d);});
        std::transform(spins.begin(),spins.end(),begin(Sy),[](double d) {return std::sin(d);});
        Sx=fourier_transform(Sx);
        Sy=fourier_transform(Sy);
        std::valarray<double> S2(spins.size());
        for(int i =0;i<spins.size();++i)
            S2[i]=std::pow(std::abs(Sx[i]),2)+std::pow(std::abs(Sy[i]),2);
        obs["|Structure Factor|^2"]<< S2;
    }

    void init_observables(alps::ObservableSet& obs) const {
        obs << alps::RealVectorObservable("|Structure Factor|^2");
    }
    //Intentionally left empty
    void save(alps::ODump &dump) const{ }
private:
    typedef typename alps::graph_helper<>::vector_type vector_type;
    std::vector<vector_type> reciprocal_vectors; 
    std::vector<vector_type> basis_vectors; 
    
    #if !NFFTW
    fftw_complex *fftw_in;
    fftw_complex *fftw_out;
    const double fftw_timelimit_; 
    fftw_plan plan; 
    #endif //NFFTW

    static constexpr std::complex<double> I=std::complex<double>(0.,1.);
    const int L,N;
    
    #if !NFFTW
    // ATTENTION: THIS IS HACKED, THERE DOES NOT EXIST A NICE VERSION OF THAT PART, 
    // BUT REMIND THAT THERE IS NO C++ GUARANTEE THAT THIS WILL WORK, 
    // IT JUST HAPPENS TO BE THE CASE, HOWEVER DO NOT ASSUME THIS TO PERSIST
    void copy_in(std::valarray<std::complex<double>> IN) {
        memcpy(fftw_in, std::begin(IN),N*sizeof(fftw_complex));
    }
    std::valarray<std::complex<double>> copy_out() {
        std::valarray<std::complex<double>> ret(N);
        memcpy(std::begin(ret), fftw_out, N*sizeof(std::complex<double>));
        return ret;
    }
    #endif //NFFTW

    #if NFFTW
    std::valarray<std::complex<double>> fourier_transform(const std::valarray<std::complex<double>>& real, int prefactor=1) const {
        std::valarray<std::complex<double>> reciprocal(N);
        for(int n1=0;n1<L;++n1)
        for(int n2=0;n2<L;++n2)
        for(int m1=0;m1<L;++m1)
        for(int m2=0;m2<L;++m2) {
             double k1=n1*(reciprocal_vectors[0])[0]+n2*(reciprocal_vectors[1])[0];
             double k2=n1*(reciprocal_vectors[0])[1]+n2*(reciprocal_vectors[1])[1];
             double r1=m1*(basis_vectors[0])[0]+m2*(basis_vectors[1])[0];
             double r2=m1*(basis_vectors[0])[1]+m2*(basis_vectors[1])[1];
             std::complex<double> pf(1,0);
             reciprocal[n1+L*n2]+=real[m1+L*m2]*std::exp(-I*pf*(k1*r1+k2*r2));
        }
        return reciprocal;
    }
    #else // !NFFTW (FFTW around)
    std::valarray<std::complex<double>> fourier_transform(const std::valarray<std::complex<double>>& real, int prefactor=1) {
        copy_in(real);
        fftw_execute(plan);
        return copy_out();
    }
    #endif //NFFTW
};
#endif //MCPP_STRUCTURE_FACTOR_H_
