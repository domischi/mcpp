#ifndef MCPP_STRUCTURE_FACTOR_H_
#define MCPP_STRUCTURE_FACTOR_H_

#include <alps/lattice.h>
#include <cassert>
#include <algorithm>
#include <valarray>
#include <alps/parameter.h>

class structure_factor {
public:
    typedef double spin_t;

    structure_factor(const alps::Parameters& p) :
    L(p["L"]),
    N(L*L)
    {
        alps::graph_helper<> gh(p);
        alps::graph_helper<>::basis_vector_iterator v,v_end;
        for(std::tie(v,v_end)=gh.reciprocal_basis_vectors() ; v!=v_end ; ++v){
            vector_type vec=*v;
            for(auto& e: vec) e/=L; //calculate the Fouriertransform of the super cell, not the elementary one...
            std::cout <<vec[0]<<std::endl;
            reciprocal_vectors.push_back(vec);
        }
        for(std::tie(v,v_end)=gh.basis_vectors() ; v!=v_end ; ++v){
            vector_type vec=*v;
            basis_vectors.push_back(vec);
        } 
        assert(reciprocal_vectors.size()==basis_vectors.size());
        assert(reciprocal_vectors.size()==2);
        assert(p["LATTICE"]=="square lattice");
    }
                
    void measure(const std::vector<spin_t>& spins, alps::ObservableSet& obs) const{
        std::valarray<std::complex<double>> Sx(spins.size());
        std::valarray<std::complex<double>> Sy(spins.size());
        std::transform(spins.begin(),spins.end(),begin(Sx),[](double d) {return std::cos(d);});
        std::transform(spins.begin(),spins.end(),begin(Sy),[](double d) {return std::sin(d);});
        Sx=fourier_transform(Sx);
        Sy=fourier_transform(Sy);
        //obs["Structure Factor X"]<<Sx; 
        //obs["Structure Factor Y"]<<Sy; 
        std::valarray<double> S2(spins.size());
        for(int i =0;i<spins.size();++i)
            S2[i]=std::pow(std::abs(Sx[i]),2)+std::pow(std::abs(Sy[i]),2);
        obs["|Structure Factor|^2"]<< S2;
    }

    void init_observables(alps::ObservableSet& obs) const {
        //obs << alps::RealVectorObservable("Structure Factor X");
        //obs << alps::RealVectorObservable("Structure Factor Y");
        obs << alps::RealVectorObservable("|Structure Factor|^2");
    }

private:
    typedef typename alps::graph_helper<>::vector_type vector_type;
    std::vector<vector_type> reciprocal_vectors; 
    std::vector<vector_type> basis_vectors; 
    static constexpr std::complex<double> I=std::complex<double>(0.,1.);
    int L,N;
    
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

};
#endif //MCPP_STRUCTURE_FACTOR_H_
