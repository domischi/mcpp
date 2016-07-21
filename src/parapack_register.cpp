#include "single-spin-flip/ssf.h"

PARAPACK_SET_VERSION("NO VERSION NUMBER AS I ONLY WORK WITH GIT");
PARAPACK_SET_COPYRIGHT("Using mc++\n  copyright (c) by Dominik Schildknecht <dominik.schildknecht@psi.ch>\n  if you reuse this project, please mention the ALPS project and me as a fair user");
PARAPACK_REGISTER_WORKER(ssf_worker, "ssf");
PARAPACK_REGISTER_EVALUATOR(ssf_evaluator, "ssf");

