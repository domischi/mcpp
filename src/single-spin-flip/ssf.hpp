#ifndef MCPP_SSF_H
#define MCPP_SSF_H

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris.h>
#include <alps/osiris/dump.h>
#include <alps/expression.h>
#include <alps/lattice.h>

class ssf : public alps::scheduler::LatticeMCRun{
public :
MyMonteCarlo(const alps::ProcessList &,const alps::Parameters &,int);
static void print_copyright(std::ostream &);
void save(alps::ODump &) const;
void load(alps::IDump &);
void dostep();
bool is_thermalized() const;
double work_done() const;
private :
// your own internal data here ...
};

typedef alps::scheduler::SimpleMCFactory<ssf> ssf_factory();
#endif //MCPP_SSF_H
