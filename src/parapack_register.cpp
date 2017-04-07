#include "version/version.h"
#include "models/xy.h"
#include <alps/parapack/exchange.h>
#include <alps/parapack/exchange_multi.h>
#include <alps/parapack/temperature_scan.h>

PARAPACK_SET_VERSION(GIT_COMMIT_HASH);
PARAPACK_SET_COPYRIGHT("Using mc++\n  copyright (c) by Dominik Schildknecht <dominik.schildknecht@psi.ch>\n  if you reuse this project, please mention the ALPS project and me\n");
//PARAPACK_REGISTER_WORKER(xy_worker, "xy");
PARAPACK_REGISTER_ALGORITHM(xy_worker, "xy");
PARAPACK_REGISTER_EVALUATOR(xy_evaluator, "xy");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<xy_worker>,"exmc");
PARAPACK_REGISTER_EVALUATOR(xy_evaluator, "exmc");
#ifdef USE_MPI
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::parallel_exchange_worker<xy_worker>,"exmc");
#endif
