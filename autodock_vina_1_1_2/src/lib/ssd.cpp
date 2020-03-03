/*

   Copyright (c) 2006-2010, The Scripps Research Institute
   Copyright (c) 2015, The University of Georgia

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
           
   Modifications for Vina-Carb 1.0 By: Anita K. Nivedha <nivedha@uga.edu>
                                       The Woods' Lab
                                       Complex Carbohydrate Research Center
                                       The University of Georgia

*/
#include "ssd.h"

// clean up
void ssd::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v, const fl chi_coeff, const fl chi_cutoff) const { // g must have correct size
	out.e = m.eval_deriv(p, ig, v, out.c, g,chi_coeff,chi_cutoff);
	fl factor = initial_factor;
	VINA_U_FOR(i, evals) {
		if(factor < min_factor) break;
		output_type candidate(out);
		candidate.c.increment(g, -factor);
		change candidate_g(g); 
		candidate.e = m.eval_deriv(p, ig, v, candidate.c, candidate_g, chi_coeff, chi_cutoff);
		if(candidate.e <= out.e) {
			out = candidate;
			g = candidate_g;
			factor *= up;
		}
		else { 
			factor *= down;
		}
	}
	out.coords = m.get_heavy_atom_movable_coords();
}
