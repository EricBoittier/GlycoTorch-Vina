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

#include "quasi_newton.h"
#include "bfgs.h"
//#include <boost/thread/mutex.hpp>

//boost::mutex cout_mutex;

struct quasi_newton_aux {
	model* m;
	const precalculate* p;
	const igrid* ig;
	const vec v;
	fl tmpq;
	const fl chi_coeff;
	const fl chi_cutoff;
	//tmp=(fl*)calloc(2,sizeof(fl));
	//std::vector< std::vector<int> > glyco_info;
	quasi_newton_aux(model* m_, const precalculate* p_, const igrid* ig_, const vec& v_/*, std::vector< std::vector<int> > glyco_info_*/, fl tmpq, fl chi_coeff_, fl chi_cutoff_) : m(m_), p(p_), ig(ig_), v(v_), chi_coeff(chi_coeff_), chi_cutoff(chi_cutoff_)/*, glyco_info(glyco_info_)*/ {}
	fl operator()(const conf& c, change& g) {
//		boost::mutex::scoped_lock lock(cout_mutex);
		//fl* tmp;
		//tmpq=(fl*)calloc(2,sizeof(fl));
		tmpq= m->eval_deriv(*p, *ig, v, c, g/*, glyco_info*/,chi_coeff,chi_cutoff);
//		std::cout<<"quasinewton tmpq =  "<<tmpq<<"\n";
		//const fl tmp = m->eval_deriv(*p, *ig, v, c, g/*, glyco_info*/)[0];
		return tmpq;
	}
};

void quasi_newton::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v/*, std::vector< std::vector<int> > glyco_info*/, fl tmp, const fl chi_coeff, const fl chi_cutoff) const { // g must have correct size
	quasi_newton_aux aux(&m, &p, &ig, v/*, glyco_info*/,tmp,chi_coeff,chi_cutoff);
//	{boost::mutex::scoped_lock lock(cout_mutex);
//	std::cout<<"before bfgs\n";
	fl res = bfgs(aux, out.c, g, max_steps, average_required_improvement, 10);
//	std::cout<<"after bfgs\n";
	out.e = res;
//	fl chi_energy = 0.0;
//        out.vce=m.eval_chi(chi_coeff,chi_cutoff);
//	std::cout<<"quasinewton out.e = "<<out.e<<" and chi energy = "<<out.vce<<"\n";
//	out.e+=out.vce;
	//out.vce=aux.tmpq[1];
//	out.vce=0;
	//free(aux.tmpq);
//	}//end of mutex
}

