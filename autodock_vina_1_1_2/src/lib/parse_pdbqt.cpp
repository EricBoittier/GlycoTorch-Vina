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

#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "parse_pdbqt.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "parse_error.h"
#include "glylib.h"

coord_3D get_puckering_parm(coord_3D **r);

std::vector<int> branch_atom1;
std::vector<int> branch_atom2;

// A new int vector to hold the size of branches
std::vector<int> branch_sizes;

std::vector<size_t*> glyco_info;
std::vector< std::vector<size_t*> > ligand_glyco_info;
int hetatm=0;
std::ofstream VC_log;

struct stream_parse_error {
	unsigned line;
	std::string reason;
	stream_parse_error(unsigned line_, const std::string& reason_) : line(line_), reason(reason_) {}
	parse_error to_parse_error(const path& name) const {
		return parse_error(name, line, reason);
	}
};


std::vector<parsed_atom> ligand_info;

void add_context(context& c, std::string& str) {
	c.push_back(parsed_line(str, boost::optional<sz>()));
}

std::string omit_whitespace(const std::string& str, sz i, sz j) {
	if(i < 1) i = 1;
	if(j < i-1) j = i-1; // i >= 1
	if(j < str.size()) j = str.size();

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	// omit trailing whitespace
	while(i <= j && std::isspace(str[j-1]))
		--j;

	VINA_CHECK(i-1 < str.size());
	VINA_CHECK(j-i+1 < str.size());

	return str.substr(i-1, j-i+1);
}

struct atom_syntax_error {
	std::string nature;
	atom_syntax_error(const std::string& nature_) : nature(nature_) {}
};

template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
	VINA_CHECK(i >= 1);
	VINA_CHECK(i <= j+1);
	if(j > str.size()) throw atom_syntax_error("The line is too short");

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	const std::string substr = str.substr(i-1, j-i+1);
	try {
		return boost::lexical_cast<T>(substr);
	}
	catch(...) {
		throw atom_syntax_error(std::string("\"") + substr + "\" is not a valid " + dest_nature);
	}
}

parsed_atom parse_pdbqt_atom_string(const std::string& str) { 
	unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "atom number");
	vec coords(checked_convert_substring<fl>(str, 31, 38, "coordinate"),
			   checked_convert_substring<fl>(str, 39, 46, "coordinate"),
			   checked_convert_substring<fl>(str, 47, 54, "coordinate"));
	fl charge = 0;
	if(!substring_is_blank(str, 69, 76))
		charge = checked_convert_substring<fl>(str, 69, 76, "charge");
	std::string name = omit_whitespace(str, 78, 79);
	std::string resname;
	std::string atomname;
	std::string resnum;
	sz ad = string_to_ad_type(name);
	resname=str.substr(17,3); 
	resnum=str.substr(23,4); 
	atomname=str.substr(13,3); 
	parsed_atom tmp(ad, charge, resname, resnum, atomname, coords, number);
	if(is_non_ad_metal_name(name))
		tmp.xs = XS_TYPE_Met_D;
	if(tmp.acceptable_type()) 
		return tmp;
	else 
		throw atom_syntax_error(std::string("\"") + name + "\" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
}

struct atom_reference {
	sz index;
	bool inflex;
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct movable_atom : public atom_vc {
	vec relative_coords;
	movable_atom(const atom_vc& a, const vec& relative_coords_) : atom_vc(a) {
		relative_coords = relative_coords_;
	}
};

struct rigid {
	atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
	vector_mutable<ligand> ligands;
	vector_mutable<residue_vc> flex;

	mav atoms;
	atomv inflex;

	distance_type_matrix atoms_atoms_bonds;
	matrix<distance_type> atoms_inflex_bonds;
	distance_type_matrix inflex_inflex_bonds;

	distance_type_matrix mobility_matrix() const {
		distance_type_matrix tmp(atoms_atoms_bonds);
		tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
		return tmp;
	}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
	struct node_t {
		sz context_index;
		parsed_atom a;
		std::vector<T> ps;
		node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

		// inflex atom insertion
		void insert_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
			nr.inflex.push_back(a);
		}
		void insert_immobiles_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile_inflex(nr);
		}

		// insertion into non_rigid_parsed
		void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
			vec relative_coords; relative_coords = a.coords - frame_origin;
			c[context_index].second = nr.atoms.size();
			nr.atoms.push_back(movable_atom(a, relative_coords));
		}
		void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile(nr, c, frame_origin);
		}
	};

	typedef node_t<parsing_struct> node; 
	boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
	boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;  

	void add(const parsed_atom& a, const context& c) { 
		VINA_CHECK(c.size() > 0);
		atoms.push_back(node(a, c.size()-1)); 
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr);
		}
	}

	// insertion into non_rigid_parsed
	void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.atoms.size(), false);
			atoms[immobile_atom.get()].insert(nr, c, frame_origin);
		}
	}

	bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
		VINA_FOR_IN(i, atoms) {
			if(immobile_atom && immobile_atom.get() != i)
				return false;
			const node& nd = atoms[i];
			if(!nd.ps.empty())
				return false; // FIXME : iffy
		}
		return true;
	}
};


int classify_ring_using_cremer_pople(coord_3D **ring)
// ring_atoms[0]=C1;
// ring_atoms[1]=C2;
// ring_atoms[2]=C3;
// ring_atoms[3]=C4;
// ring_atoms[4]=C5;
// ring_atoms[5]=O5;
{



coord_3D p ;
p = get_puckering_parm(ring);

double phi;
double theta;
double q;

phi = p.k * 180/PI ;
theta = p.j  * 180/PI;
q = p.i;

std::cout  << "Phi: " << phi << " Theta: " << theta << " Q: " << q << " Ring: ";


if (theta < 40.00) {
	std::cout  << " 4c1 \n";
	return 1;
}

if (theta > 140.00) {
	std::cout  << " 1c4 \n";
	return 2;
}

if ((112.5 > theta) && (theta > 67.5) && (130 < phi) && (phi < 170)) {
	std::cout  << " 2sO \n";
	return 3;
}

std::cout  << " ???\n ";

return 0;
}


coord_3D get_puckering_parm(coord_3D **r){
int n = 6;
int jval[n];
int l;
double z[n];
double H1[n];
double H2[n];
double q[n];
double m1=0.0;
double n1=0.0;
double q1=0.0;
double M=0.0;
double N=0.0;
plane pval;//plane for returning the values of plane 
vectormag_3D Rj,Rcos,Rsin,R1,R2,R1xR2,avg_coord;
R1.i=R1.j=R1.k=R2.i=R2.j=R2.k=0;
avg_coord.i=avg_coord.j=avg_coord.k=0;
for(l=0;l<n;l++){
	jval[l]=l+1;//getting j vals for further equations
	}
for(l=0;l<n;l++){//for loop for calculating the Rj vals
	Rj.i= r[l]->i;
	Rj.j= r[l]->j;
	Rj.k= r[l]->k;
	Rsin.i = (Rj.i*sin(2*PI*(jval[l]-1)/n));
	Rsin.j = (Rj.j*sin(2*PI*(jval[l]-1)/n));
	Rsin.k = (Rj.k*sin(2*PI*(jval[l]-1)/n));
	Rcos.i = (Rj.i*cos(2*PI*(jval[l]-1)/n));
	Rcos.j = (Rj.j*cos(2*PI*(jval[l]-1)/n));
	Rcos.k = (Rj.k*cos(2*PI*(jval[l]-1)/n));
	R1.i=R1.i+Rsin.i;
	R1.j=R1.j+Rsin.j;
	R1.k=R1.k+Rsin.k;
	R2.i=R2.i+Rcos.i;
	R2.j=R2.j+Rcos.j;
	R2.k=R2.k+Rcos.k;
	}
R1xR2 = get_crossprod(R1,R2);
pval.A = (R1xR2.i/R1xR2.d);
pval.B = (R1xR2.j/R1xR2.d);
pval.C = (R1xR2.k/R1xR2.d);
for(l=0;l<n;l++){//for loop for calculating the avg x,y,z coordinates
	avg_coord.i=avg_coord.i+r[l]->i;
	avg_coord.j=avg_coord.j+r[l]->j;
	avg_coord.k=avg_coord.k+r[l]->k;
	}
avg_coord.i=(avg_coord.i/n);
avg_coord.j=(avg_coord.j/n);
avg_coord.k=(avg_coord.k/n);

pval.D= -(pval.A*avg_coord.i+pval.B*avg_coord.j+pval.C*avg_coord.k);

for(l=0;l<n;l++){
	z[l]=(pval.A*r[l]->i+pval.B*r[l]->j+pval.C*r[l]->k+pval.D)/1;
	H1[l]=(z[l]*cos(4*PI*(jval[l]-1)/n));
	H2[l]=(z[l]*sin(4*PI*(jval[l]-1)/n));
	q[l]=z[l]*(pow((-1),jval[l]-1));
	m1=m1+H1[l];
	n1=n1+H2[l];
	q1=q1+q[l];
	}
double x=(1.0/3.0);
double y=(1.0/6.0);
coord_3D *parm; 
parm=(coord_3D*)calloc(2,sizeof(coord_3D));
M=m1*(sqrt(x));
N=-(n1*(sqrt(x)));

/*
        parm[0].i = q_2
        parm[0].j = q_3
        parm[0].k = phi_2
        parm[1].i = Q
        parm[1].j = theta
        parm[1].k = phi
*/
parm[0].k=atan(N/M);  /* phi_2 */
parm[0].i = N/sin(parm[0].k); /* q_2 */
if(parm[0].i<0)
    {/* if the q_2<0 */
    parm[0].k += PI; /* rotate phi_2 by PI */
    parm[0].i = -parm[0].i; /* make q_2 positive */
    } 
if(parm[0].k<0)
    {/* if phi_2 is negative */
    parm[0].k += 2*PI; /* enforce acceptable range of phi_2 */
    }
parm[0].j = q1*(sqrt(y)); /* q_3 */ 
/* */
parm[1].i = sqrt ( M*M + N*N + parm[0].j*parm[0].j ); /* Q */
parm[1].k = parm[0].k;  /* phi = phi_2 */
parm[1].j = acos(parm[0].j/parm[1].i); /* theta = acos ( q_3 / Q ) */

return parm[1];

}

unsigned parse_one_unsigned(const std::string& str, const std::string& start, unsigned count) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp;
	in_str >> tmp;
	if(!in_str || tmp < 0) 
		throw stream_parse_error(count, "Syntax error");
	return unsigned(tmp);
}

void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned count, unsigned& first, unsigned& second) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp1, tmp2;
	in_str >> tmp1;
	in_str >> tmp2;
	if(!in_str || tmp1 < 0 || tmp2 < 0) 
		throw stream_parse_error(count, "Syntax error");
	first = unsigned(tmp1);
	second = unsigned(tmp2);
}

void parse_pdbqt_rigid(const path& name, rigid& r) {
	ifile in(name);
	unsigned count = 0;
	std::string str;
	while(std::getline(in, str)) {
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "TER")) {} // ignore 
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				r.atoms.push_back(parse_pdbqt_atom_string(str));
			}
			catch(atom_syntax_error& e) {
				throw parse_error(name, count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw parse_error(name, count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
	unsigned line_count = 0;
	// a loop to add branch sizes
	while(std::getline(in, str)) {
		++line_count;
	}





}

double get_torsion_coords_vec_list(vec A, vec B, vec C, vec D)
{ //returns the torsion angle formed by 4 x,y,z co-ordinates.
double angle=0.0;
vec AB, BC, CD, ABCcross, BCDcross;
double ABC_BCD_dot, AB_BCD_dot, BCscalar_AB_BCD_dot;
AB=B-A;
BC=C-B;
CD=D-C;
ABCcross=cross_product(AB,BC);
BCDcross=cross_product(BC,CD);
ABC_BCD_dot=dot_product(ABCcross,BCDcross);
AB_BCD_dot=dot_product(AB,BCDcross);
BCscalar_AB_BCD_dot=magnitude(BC)*AB_BCD_dot;                        
angle=atan2(BCscalar_AB_BCD_dot,ABC_BCD_dot)*57.2957795; //angle in degrees
return angle;           
}                       

void parse_pdbqt_root_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				if(hetatm==1){
				parsed_atom LIGpa=parse_pdbqt_atom_string(str);//AKN
				ligand_info.push_back(LIGpa); 
				}
				p.add(parse_pdbqt_atom_string(str), c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "ENDROOT")) return;
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_root(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ROOT")) {
			parse_pdbqt_root_aux(in, count, p, c);
			break;
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux(std::istream& in, unsigned& count, const std::string& str, parsing_struct& p, context& c) {
	unsigned first, second;
	parse_two_unsigneds(str, "BRANCH", count, first, second); 
	sz i = 0;
	if(hetatm==1){
	branch_atom1.push_back(first);
	branch_atom2.push_back(second);
	}
	for(; i < p.atoms.size(); ++i)
		{
		if(p.atoms[i].a.number == first) {
			p.atoms[i].ps.push_back(parsing_struct());
			parse_pdbqt_branch(in, count, p.atoms[i].ps.back(), c, first, second);
			break;
		}
		}
	if(i == p.atoms.size())
		throw stream_parse_error(count, "No atom number " + boost::lexical_cast<std::string>(first) + " in this branch");
}

              

void parse_pdbqt_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue_vc) {
	
	parse_pdbqt_root(in, count, p, c);

	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(!residue_vc && starts_with(str, "TORSDOF")) {
			if(torsdof) throw stream_parse_error(count, "TORSDOF can occur only once");
			torsdof = parse_one_unsigned(str, "TORSDOF", count);
		}
		else if(residue_vc && starts_with(str, "END_RES")) return; 
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}

	int branch=0, i=0;
	VINA_FOR(i,ligand_info.size())	
	{//This loop helps to remove blank spaces from atomnames in case they have them so the actual names can be used for comparison instead of with trailing blank spaces.
		if(ligand_info[i].atomname.find(" ")<3)
		{
		ligand_info[i].atomname.replace(ligand_info[i].atomname.find(" "),1,"\0");
		}
	}

	std::cout << "\n";

	VINA_FOR(h,ligand_info.size())
	{
		std::cout << h << " " << ligand_info[h].atomname << "\n";
	}

	VINA_FOR(h,branch_atom1.size())
	{
	int ring1_conf=0, ring2_conf=0;
	size_t *S1_C1, *S1_O5, *S1_O1, *S1_C5, *S1_C2, *S1_C3, *S1_C4, *S2_Cxp1, *S2_C2, *S2_O5, *S2_C5, *S2_C1, *S2_O1, *S2_Ox, *S2_Cx, *S2_Cxm1, *S2_O4, *S2_C4, *S2_C3, *sizet_ring1_conf;
	size_t *S1_AB, *S2_AE, *S2_Link, *S2_6AE, *sizet_ring2_conf;
	coord_3D *C1, *C2, *C3, *C4, *C5, *O5, **ring_atoms;
	coord_3D *C1_, *C2_, *C3_, *C4_, *C5_, *O5_, **cp_ring_atoms;
    double S1_AB_angle=0.0, S2_AE_angle=0.0, S2_omega_angle=0.0;
	ring_atoms=(coord_3D**)calloc(6,sizeof(coord_3D*));
	cp_ring_atoms=(coord_3D**)calloc(6,sizeof(coord_3D*));
	C1=(coord_3D*)calloc(1,sizeof(coord_3D));
	C2=(coord_3D*)calloc(1,sizeof(coord_3D));
	C3=(coord_3D*)calloc(1,sizeof(coord_3D));
	C4=(coord_3D*)calloc(1,sizeof(coord_3D));
	C5=(coord_3D*)calloc(1,sizeof(coord_3D));
	O5=(coord_3D*)calloc(1,sizeof(coord_3D));
	S1_C1=(size_t*)calloc(1,sizeof(size_t));
	S1_O5=(size_t*)calloc(1,sizeof(size_t));
	S1_O1=(size_t*)calloc(1,sizeof(size_t));
	S1_C5=(size_t*)calloc(1,sizeof(size_t));
	S1_C2=(size_t*)calloc(1,sizeof(size_t));
	S1_C3=(size_t*)calloc(1,sizeof(size_t));
	S1_C4=(size_t*)calloc(1,sizeof(size_t));
	sizet_ring1_conf=(size_t*)calloc(1,sizeof(size_t));
	S2_Cxp1=(size_t*)calloc(1,sizeof(size_t));
	S2_C2=(size_t*)calloc(1,sizeof(size_t));
	S2_O5=(size_t*)calloc(1,sizeof(size_t));
	S2_O4=(size_t*)calloc(1,sizeof(size_t));
	S2_C4=(size_t*)calloc(1,sizeof(size_t));
	S2_C3=(size_t*)calloc(1,sizeof(size_t));
	S2_C5=(size_t*)calloc(1,sizeof(size_t));
	S2_C1=(size_t*)calloc(1,sizeof(size_t));
	S2_O1=(size_t*)calloc(1,sizeof(size_t));
	S2_Ox=(size_t*)calloc(1,sizeof(size_t));
	S2_Cx=(size_t*)calloc(1,sizeof(size_t));
	S2_Cxm1=(size_t*)calloc(1,sizeof(size_t));
	S1_AB=(size_t*)calloc(1,sizeof(size_t));
	S2_AE=(size_t*)calloc(1,sizeof(size_t));
	S2_Link=(size_t*)calloc(1,sizeof(size_t));
	S2_6AE=(size_t*)calloc(1,sizeof(size_t));
	sizet_ring2_conf=(size_t*)calloc(1,sizeof(size_t));
	S1_AB[0]=-1;
	S2_6AE[0]=-1;

	//  Must find these atoms before adding to glycoinfo
	bool S1_C1_, S1_C2_, S1_C3_, S1_C4_, S1_C5_, S1_O5_, S2_Ox_, S2_C1_, S2_C2_, S2_C3_, S2_C4_, S2_C5_, S2_O5_ ;
	S1_C1_ = S1_C2_ = S1_C3_ = S1_C4_ = S1_C5_ = S1_O5_ = S2_Ox_ = S2_C1_ = S2_C2_ = S2_C3_ = S2_C4_ = S2_C5_ = S2_O5_ = false;

	std::cout << h << " " << branch_atom1[h] - 1 << " " << ligand_info[branch_atom1[h]-1].atomname << " " << ligand_info[branch_atom1[h]-1].resnum 
	<< " " << branch_atom2[h] - 1 << " " << ligand_info[branch_atom2[h]-1].atomname << " " << ligand_info[branch_atom2[h]-1].resnum << "\n";

    if(
	// if the two atoms from the BRANCH line are from different residues 
	(ligand_info[branch_atom1[h]-1].resnum!=ligand_info[branch_atom2[h]-1].resnum) && 

	(ligand_info[branch_atom1[h]-1].resname.compare("OME")!=0 && 
	ligand_info[branch_atom1[h]-1].resname.compare("ROH")!=0 ) && 
	(ligand_info[branch_atom2[h]-1].resname.compare("OME")!=0 && 
	ligand_info[branch_atom2[h]-1].resname.compare("ROH")!=0 )  )

	{//found glycosidic

	std::cout << "\n";

	VINA_FOR(i,2){if(i==0){branch=branch_atom1[h];} else if(i==1){branch=branch_atom2[h];} 
	//using 2 values for i, one for branch_atom1 and another for branch_atom2

		if(ligand_info[branch-1].atomname.compare("C1")==0){S1_C1[0]=branch-1; S1_C1_=true; 
		// loop through the entire atom list to find atoms in sugar 1 (S1)
			VINA_FOR(j,ligand_info.size())
			{if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
				{if(ligand_info[j].atomname.compare("C2")==0){S1_C2[0]=j; S1_C2_=true;}
				if(ligand_info[j].atomname.compare("O5")==0){S1_O5[0]=j; S1_O5_=true;}
				if(ligand_info[j].atomname.compare("C5")==0){S1_C5[0]=j; S1_C5_=true;}
				if(ligand_info[j].atomname.compare("C4")==0){S1_C4[0]=j; S1_C4_=true;}
				if(ligand_info[j].atomname.compare("C3")==0){S1_C3[0]=j; S1_C3_=true;}
				}}}

			else if(ligand_info[branch-1].atomname[0]=='O'){S2_Ox[0]=branch-1; S2_Ox_=true;
			// loop through the entire atom list to find atoms in sugar 2 (S2)
			VINA_FOR(j,ligand_info.size()){
			if(ligand_info[j].resnum==ligand_info[branch-1].resnum){
				if(ligand_info[j].atomname[0]=='C')
					{if(ligand_info[j].atomname[1]==ligand_info[branch-1].atomname[1]) 
						{S2_Cx[0]=j;
						if(ligand_info[j].atomname[1]=='2'){S2_Link[0]=2;}
						if(ligand_info[j].atomname[1]=='3'){S2_Link[0]=3;}
						if(ligand_info[j].atomname[1]=='4'){S2_Link[0]=4;}
						if(ligand_info[j].atomname[1]=='6'){S2_Link[0]=6;}}

					if(ligand_info[j].atomname.compare("C2")==0){S2_C2[0]=j; S2_C2_=true;}
					if(ligand_info[j].atomname.compare("C5")==0){S2_C5[0]=j; S2_C5_=true;}
					if(ligand_info[j].atomname.compare("C4")==0){S2_C4[0]=j; S2_C4_=true;}
					if(ligand_info[j].atomname.compare("C3")==0){S2_C3[0]=j; S2_C3_=true;}
					if(ligand_info[j].atomname.compare("C1")==0){S2_C1[0]=j; S2_C1_=true;}}

				if(ligand_info[j].atomname.compare("O5")==0){S2_O5[0]=j; S2_O5_=true;}}
				} 

		// 	float distance=0.0, distance_1up=100.0, distance_1down=100.0;
		// 		VINA_FOR(j,ligand_info.size()){
		// 		if(ligand_info[j].resnum==ligand_info[branch-1].resnum){
		// 			if(ligand_info[j].atomname[0]=='C'){
		// 				distance=sqrt(pow(ligand_info[j].coords.data[0]-ligand_info[S2_Cx[0]].coords.data[0],2)+pow(ligand_info[j].coords.data[1]-ligand_info[S2_Cx[0]].coords.data[1],2)+pow(ligand_info[j].coords.data[2]-ligand_info[S2_Cx[0]].coords.data[2],2));
		// 					if(distance<distance_1up && ligand_info[j].atomname.compare(ligand_info[S2_Cx[0]].atomname)>0)
		// 					{
		// 					distance_1up=distance;
		// 					S2_Cxp1[0]=j;
		// 					}
		// 					if(distance<distance_1down && ligand_info[j].atomname.compare(ligand_info[S2_Cx[0]].atomname)<0)
		// 					{
		// 					distance_1down=distance;
        //                                                 S2_Cxm1[0]=j;
		// 					}
		// 				}
		// 			}
		// 		}
		// 		float distance_O1=100.0;
		// 		VINA_FOR(j,ligand_info.size())
		// 		{
		// 			if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
        //                                 {
		// 				if(ligand_info[j].atomname.compare("O5")==0)
		// 				{
		// 				S2_O5[0]=j;
		// 				}
		// 				if(ligand_info[j].atomname.compare("O1")==0)
        //                                         {//To find S2_O1 -- if present inside same residue as O1
        //                                         S2_O1[0]=j;
        //                                         }
		// 				if(ligand_info[j].atomname.compare("O4")==0)
        //                                         {//To find S2_O1 -- if present inside same residue as O1
        //                                         S2_O4[0]=j;
        //                                         }
        //                                         if(ligand_info[j].atomname.compare("C4")==0)
        //                                         {//To find S2_O1 -- if present inside same residue as O1
        //                                         S2_C4[0]=j;
        //                                         }
        //                                         if(ligand_info[j].atomname.compare("C3")==0)
        //                                         {//To find S2_O1 -- if present inside same residue as O1
        //                                         S2_C3[0]=j;
        //                                         }
		// 				else
		// 				{//To find S2_O1 -- when having to use Ox of neighbouring residue in it's place
		// 					VINA_FOR(k,ligand_info.size())
		// 					{//Going through all atoms im structure so that nearest linking Oxygen atom can be found
		// 						if(ligand_info[k].atomname[0]=='O' && ligand_info[k].atomname[1]!='5')
		// 						{
		// 						distance=sqrt(pow(ligand_info[k].coords.data[0]-ligand_info[S2_C1[0]].coords.data[0],2)+pow(ligand_info[k].coords.data[1]-ligand_info[S2_C1[0]].coords.data[1],2)+pow(ligand_info[k].coords.data[2]-ligand_info[S2_C1[0]].coords.data[2],2));
		// 							if(distance<distance_O1)
		// 							{
		// 							distance_O1=distance;
		// 							S2_O1[0]=k;
		// 							}
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}
		// 	} //end of if where Ox was found
		}//end of if i is either 0 or 1!
	}
		if ((S1_C1_) && (S1_C2_) && (S1_C3_) && (S1_C4_) && (S1_C5_) && (S1_O5_) && (S2_Ox_) 
		&& (S2_C1_) && (S2_C2_) && (S2_C3_) && (S2_C4_) && (S2_C5_) && (S2_O5_)) {

				std::cout << "here";
				C1[0].i=ligand_info[S1_C1[0]].coords[0];
				C1[0].j=ligand_info[S1_C1[0]].coords[1];
				C1[0].k=ligand_info[S1_C1[0]].coords[2];

				C2[0].i=ligand_info[S1_C2[0]].coords[0];
				C2[0].j=ligand_info[S1_C2[0]].coords[1];
				C2[0].k=ligand_info[S1_C2[0]].coords[2];

				C3[0].i=ligand_info[S1_C3[0]].coords[0];
				C3[0].j=ligand_info[S1_C3[0]].coords[1];
				C3[0].k=ligand_info[S1_C3[0]].coords[2];

				C4[0].i=ligand_info[S1_C4[0]].coords[0];
				C4[0].j=ligand_info[S1_C4[0]].coords[1];
				C4[0].k=ligand_info[S1_C4[0]].coords[2];

				C5[0].i=ligand_info[S1_C5[0]].coords[0];
				C5[0].j=ligand_info[S1_C5[0]].coords[1];
				C5[0].k=ligand_info[S1_C5[0]].coords[2];

				O5[0].i=ligand_info[S1_O5[0]].coords[0];
				O5[0].j=ligand_info[S1_O5[0]].coords[1];
				O5[0].k=ligand_info[S1_O5[0]].coords[2];
				
				S1_AB_angle=get_torsion_coords_vec_list(ligand_info[S1_C5[0]].coords,ligand_info[S1_O5[0]].coords,ligand_info[S1_C1[0]].coords,ligand_info[S2_Ox[0]].coords);
				S2_omega_angle=get_torsion_coords_vec_list(ligand_info[S2_O4[0]].coords,ligand_info[S2_C4[0]].coords,ligand_info[S2_C3[0]].coords,ligand_info[S2_C5[0]].coords);
				S2_AE_angle=get_angle_ABC(ligand_info[S2_Ox[0]].coords,ligand_info[S2_Cx[0]].coords,ligand_info[S2_O5[0]].coords);


				ring_atoms[0]=C1;
				ring_atoms[1]=C2;
				ring_atoms[2]=C3;
				ring_atoms[3]=C4;
				ring_atoms[4]=C5;
				ring_atoms[5]=O5;

				cp_ring_atoms[0]=O5;
				cp_ring_atoms[1]=C1;
				cp_ring_atoms[2]=C2;
				cp_ring_atoms[3]=C3;
				cp_ring_atoms[4]=C4;
				cp_ring_atoms[5]=C5;

				// Printing
				std::cout << ligand_info[S1_C1[0]].resname << " " << ligand_info[S1_C1[0]].resnum << " ";

				// std::cout << "\n " << branch_atom1.size() << "\n";

				ring1_conf=classify_ring_using_cremer_pople(cp_ring_atoms);

				sizet_ring1_conf[0]=ring1_conf;
				C1[0].i=ligand_info[S2_C1[0]].coords[0];
				C1[0].j=ligand_info[S2_C1[0]].coords[1];
				C1[0].k=ligand_info[S2_C1[0]].coords[2];

				C2[0].i=ligand_info[S2_C2[0]].coords[0];
				C2[0].j=ligand_info[S2_C2[0]].coords[1];
				C2[0].k=ligand_info[S2_C2[0]].coords[2];

				C3[0].i=ligand_info[S2_C3[0]].coords[0];
				C3[0].j=ligand_info[S2_C3[0]].coords[1];
				C3[0].k=ligand_info[S2_C3[0]].coords[2];

				C4[0].i=ligand_info[S2_C4[0]].coords[0];
				C4[0].j=ligand_info[S2_C4[0]].coords[1];
				C4[0].k=ligand_info[S2_C4[0]].coords[2];

				C5[0].i=ligand_info[S2_C5[0]].coords[0];
				C5[0].j=ligand_info[S2_C5[0]].coords[1];
				C5[0].k=ligand_info[S2_C5[0]].coords[2];

				O5[0].i=ligand_info[S2_O5[0]].coords[0];
				O5[0].j=ligand_info[S2_O5[0]].coords[1];
				O5[0].k=ligand_info[S2_O5[0]].coords[2];

				ring_atoms[0]=C1;
				ring_atoms[1]=C2;
				ring_atoms[2]=C3;
				ring_atoms[3]=C4;
				ring_atoms[4]=C5;
				ring_atoms[5]=O5;

				cp_ring_atoms[0]=O5;
				cp_ring_atoms[1]=C1;
				cp_ring_atoms[2]=C2;
				cp_ring_atoms[3]=C3;
				cp_ring_atoms[4]=C4;
				cp_ring_atoms[5]=C5;

				std::cout << ligand_info[S2_C1[0]].resname << " " << ligand_info[S2_C1[0]].resnum << " ";

				// std::cout << "\n " << branch_atom2.size() << "\n";

				ring2_conf=classify_ring_using_cremer_pople(cp_ring_atoms);

				std::cout << "\n";

				sizet_ring2_conf[0]=ring2_conf;
				if(sizet_ring1_conf[0]==0)
				{
				VC_log<<"CHI energy penalties NOT applied to phi torsion in linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring2_conf[0]==0)
				{
				VC_log<<"CHI energy penalties NOT applied to psi torsion in linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring1_conf[0]==0 && sizet_ring2_conf[0]==0)
				{
				VC_log<<"CHI energy penalities NOT applied to glycosidic linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring1_conf[0]!=0 && sizet_ring2_conf[0]!=0)
				{
				VC_log<<"CHI energy penalties applied to "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<" linkage.\n";
				}
				if( ((S1_AB_angle<-40) && (S1_AB_angle>-80)) || ((S1_AB_angle>40) && (S1_AB_angle<80)) )
				{
					//Phi_Alpha_L
					S1_AB[0]=0;
				}
				else if( ((S1_AB_angle<-160) && (S1_AB_angle>-200)) || ((S1_AB_angle>160) && (S1_AB_angle<200)) )
				{
                                        //Phi_Beta_D
					S1_AB[0]=1;
				}
				if( (S2_AE_angle<120 && S2_AE_angle>80) )
				{
				//Axial attachment
				S2_AE[0]=0;
				}
				else if(S2_AE_angle>130 && S2_AE_angle<170)
				{
				//Equatorial attachment
				S2_AE[0]=1;
				}
				if(S2_omega_angle>0){S2_6AE[0]=1;}
				else if (S2_omega_angle<0){S2_6AE[0]=0;}
		//Co-ordinates START
		glyco_info.push_back(S1_O5); //index 0
		glyco_info.push_back(S1_C1); //index 1
		glyco_info.push_back(S2_Ox); //index 2
		glyco_info.push_back(S2_Cx); //index 3
		glyco_info.push_back(S2_Cxm1); //index 4
		//Co-ordinates END
		glyco_info.push_back(S1_AB); //index 5
		glyco_info.push_back(S2_AE); //index 6
		glyco_info.push_back(S2_Link); //index 7 
        glyco_info.push_back(S2_6AE); //index 8  //0 -> positive; 1-> negative
        //Co-ordinates START
    	glyco_info.push_back(S2_O5); //index 9 
		glyco_info.push_back(sizet_ring1_conf); //index 10
		glyco_info.push_back(sizet_ring2_conf); //index 11
        //Co-ordinates ED
		ligand_glyco_info.push_back(glyco_info);
		glyco_info.clear();
		}
	}//end of finding glycosidic

	}//end of for for branch atom size
}


std::vector<parsed_atom> liginfo_return()
{
return ligand_info;
}

std::vector< std::vector<size_t*> > glycan_info_func()
{
return ligand_glyco_info;
}

void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
	if(atm)
		VINA_RANGE(i, r.begin, r.end) {
			atom_reference& ar = atm.get();
			if(ar.inflex) 
				nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
			else
				nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
		}
}

void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
	if(axis_begin && axis_end) {
		atom_reference& r1 = axis_begin.get();
		atom_reference& r2 = axis_end  .get();
		if(r2.inflex) {
			VINA_CHECK(r1.inflex); // no atom-inflex rotors
			nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
		}
		else
			if(r1.inflex)
				nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
			else
				nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
	}
}

typedef std::pair<sz, sz> axis_numbers;
typedef boost::optional<axis_numbers> axis_numbers_option;

void nr_update_matrixes(non_rigid_parsed& nr) {
	// atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)

	nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
	nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
	nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

template<typename B> // B == branch or main_branch or flexible_body 
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
	b.node.begin = nr.atoms.size();
	VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
		parsing_struct::node& p_node = p.atoms[i];
		if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
		else p_node.insert(nr, c, b.node.get_origin());
		p_node.insert_immobiles(nr, c, b.node.get_origin());
	}
	b.node.end = nr.atoms.size();
	nr_update_matrixes(nr);
	add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
	add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
	set_rotor(nr, p.axis_begin, p.axis_end);

	VINA_RANGE(i, b.node.begin, b.node.end)
		VINA_RANGE(j, i+1, b.node.end)
			nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME


	VINA_FOR_IN(i, p.atoms) { 	// postprocess children
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
				b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, b.children.back());
			}
		}
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
	VINA_CHECK(!p.atoms.empty());
	nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
	postprocess_branch(nr, p, c, nr.ligands.back());
	nr_update_matrixes(nr); // FIXME ?
}

void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		p_node.insert_inflex(nr);
		p_node.insert_immobiles_inflex(nr);
	}
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
				nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_//branch will assign begin and end
				postprocess_branch(nr, ps, c, nr.flex.back());
			}
		}
	}
	nr_update_matrixes(nr); // FIXME ?
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	VC_log.open("VC_log.txt");
	unsigned count = 0;
	parsing_struct p;
	boost::optional<unsigned> torsdof;
	try {
		parse_pdbqt_aux(in, count, p, c, torsdof, false);
		if(p.atoms.empty()) 
			throw parse_error(name, count, "No atoms in the ligand");
		if(!torsdof)
			throw parse_error(name, count, "Missing TORSDOF");
		postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
	}
	catch(stream_parse_error& e) {
		throw e.to_parse_error(name);
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VC_log.close();
}

void parse_pdbqt_residue(std::istream& in, unsigned& count, parsing_struct& p, context& c) { 
	boost::optional<unsigned> dummy;
	parse_pdbqt_aux(in, count, p, c, dummy, true);
}

void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	unsigned count = 0;
	std::string str;

	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BEGIN_RES")) {
			try {
				parsing_struct p;
				parse_pdbqt_residue(in, count, p, c);
				postprocess_residue(nr, p, c);
			}
			catch(stream_parse_error& e) {
				throw e.to_parse_error(name);
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} //ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(starts_with(str, "ENDBRANCH")) {
			unsigned first, second;
			parse_two_unsigneds(str, "ENDBRANCH", count, first, second);
			if(first != from || second != to) 
				throw stream_parse_error(count, "Inconsistent branch numbers");
			if(!p.immobile_atom) 
				throw stream_parse_error(count, "Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch");
			return;
		}
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				parsed_atom a = parse_pdbqt_atom_string(str);
				
				if(hetatm==1){
				ligand_info.push_back(a);
				}
				if(a.number == to)
					p.immobile_atom = p.atoms.size();
				p.add(a, c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}


//////////// new stuff //////////////////


struct pdbqt_initializer {
	model m;
	void initialize_from_rigid(const rigid& r) { // static really
		VINA_CHECK(m.grid_atoms.empty());
		m.grid_atoms = r.atoms;
	}
	void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
		VINA_CHECK(m.ligands.empty());
		VINA_CHECK(m.flex   .empty());

		m.ligands = nrp.ligands;
		m.flex    = nrp.flex;

		VINA_CHECK(m.atoms.empty());

		sz n = nrp.atoms.size() + nrp.inflex.size();
		m.atoms.reserve(n);
		m.coords.reserve(n);

		VINA_FOR_IN(i, nrp.atoms) {
			const movable_atom& a = nrp.atoms[i];
			atom_vc b = static_cast<atom_vc>(a);
			b.coords = a.relative_coords;
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_FOR_IN(i, nrp.inflex) {
			const atom_vc& a = nrp.inflex[i];
			atom_vc b = a;
			b.coords = zero_vec_vc; // to avoid any confusion; presumably these will never be looked at
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_CHECK(m.coords.size() == n);

		m.internal_coords.resize(m.coords.size(), zero_vec_vc); // FIXME

		m.minus_forces = m.coords;
		m.m_num_movable_atoms = nrp.atoms.size();

		if(is_ligand) {
			VINA_CHECK(m.ligands.size() == 1);
			m.ligands.front().cont = c;
		}
		else
			m.flex_context = c;

	}
	void initialize(const distance_type_matrix& mobility) {
		m.initialize(mobility);
	}
};

model parse_ligand_pdbqt  (const path& name) { // can throw parse_error
	non_rigid_parsed nrp;
	context c;
	hetatm=1;
	parse_pdbqt_ligand(name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_nrp(nrp, c, true);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name, const path& flex_name) { // can throw parse_error
	hetatm=0;
	rigid r;
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_rigid(rigid_name, r);
	parse_pdbqt_flex(flex_name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	tmp.initialize_from_nrp(nrp, c, false);
	tmp.initialize(nrp.mobility_matrix());
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name) { // can throw parse_error
	hetatm=0;
	rigid r;
	parse_pdbqt_rigid(rigid_name, r);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	distance_type_matrix mobility_matrix;
	tmp.initialize(mobility_matrix);
	return tmp.m;

}


