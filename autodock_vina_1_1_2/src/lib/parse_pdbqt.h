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

#ifndef VINA_PARSE_PDBQT_H
#define VINA_PARSE_PDBQT_H

#include "model.h"


struct parsed_atom : public atom_vc {
	unsigned number; 
	std::string resname;
	std::string resnum; 
	std::string atomname;
	parsed_atom(sz ad_ , fl charge_, const std::string& resname, const std::string& resnum, const std::string& atomname, const vec& coords_, unsigned number_) :  number(number_) , resname(resname) , resnum(resnum) , atomname(atomname) {
		ad = ad_;
		charge = charge_;
		coords = coords_;
	}
};

model parse_receptor_pdbqt(const path& rigid, const path& flex); // can throw parse_error
model parse_receptor_pdbqt(const path& rigid); // can throw parse_error
model parse_ligand_pdbqt  (const path& name); // can throw parse_error
std::vector<parsed_atom> liginfo_return(); //AKN
std::vector< std::vector<size_t*> > glycan_info_func();
#endif

