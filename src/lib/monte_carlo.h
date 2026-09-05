/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "ssd.h"
#include "incrementable.h"

// Which search this binary is. One tree builds several: the targets differ only
// by this define on the compile line, the way vina and vina_split differ only by
// their main. `vina` is stock AutoDock Vina; the QuickVina targets switch on the
// history buffer described in search_database.h.
//
//     -DQVINA_SEARCH=SEARCH_QVINA_2
//
#define SEARCH_VINA    0
#define SEARCH_QVINA_2 1
#define SEARCH_QVINA_W 2

#ifndef QVINA_SEARCH
#define QVINA_SEARCH SEARCH_VINA
#endif

struct monte_carlo {
	unsigned num_steps;
	bool use_database; // consult the history buffer before a local search
	fl temperature;
	vec hunt_cap;
	fl min_rmsd;
	sz num_saved_mins;
	fl mutation_amplitude;
	ssd ssd_par;
	monte_carlo() : num_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) , use_database(QVINA_SEARCH != SEARCH_VINA) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

	output_type operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	output_type many_runs(model& m, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;

	void single_run(model& m, output_type& out, const precalculate& p, const igrid& ig, rng& generator) const;
	// out is sorted
	void operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const;
	void many_runs(model& m, output_container& out, const precalculate& p, const igrid& ig, const vec& corner1, const vec& corner2, sz num_runs, rng& generator) const;

};

#endif

