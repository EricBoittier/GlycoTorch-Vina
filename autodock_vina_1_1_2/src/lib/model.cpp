/*

   Copyright (c) 2006+6010, The Scripps Research Institute
   Copyright (c) 2015, The University of Georgia

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE+6.0

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

#include "model.h"
#include "file.h"
#include "curl.h"
#include "parse_pdbqt.h"
#include <boost/thread/mutex.hpp>
namespace
	{
	boost::mutex cout_mutex;
	}



template<typename T>
atom_range get_atom_range(const T& t) {
	atom_range tmp = t.node;
	VINA_FOR_IN(i, t.children) {
		atom_range r = get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}

struct branch_metrics {
	sz length;
	sz corner2corner;
	branch_metrics() : length(0), corner2corner(0) {}
};

template<typename T>
branch_metrics get_branch_metrics(const T& t) {
	branch_metrics tmp;
	if(!t.children.empty()) {
		sz corner2corner_max = 0;
		szv lengths;
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);
			if(corner2corner_max < res.corner2corner)
				corner2corner_max = res.corner2corner;
			lengths.push_back(res.length + 1); // FIXME? weird compiler warning (sz -> unsigned)
		}
		std::sort(lengths.begin(), lengths.end());

		tmp.length = lengths.back();

		tmp.corner2corner = tmp.length;
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 1];

		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

sz model::ligand_longest_branch(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).length;
}

sz model::ligand_length(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

void ligand::set_range() {
	atom_range tmp = get_atom_range(*this);
	begin = tmp.begin;
	end   = tmp.end;
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

struct appender_info {
	sz grid_atoms_size;
	sz m_num_movable_atoms;
	sz atoms_size;

	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

class appender {
	appender_info a_info;
	appender_info b_info;
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
	bool is_a;

	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
		if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);
		else            tmp.i = operator()(tmp.i);
			return tmp;
	}

	// type-directed old -> new transformations
	void update(interacting_pair& ip) const {
		ip.a = operator()(ip.a);
		ip.b = operator()(ip.b);
	}
	void update(vec& v) const { // coordinates & forces - do nothing
	}
	void update(ligand& lig) const {
		lig.transform(*this); // ligand as an atom_range subclass
		transform_ranges(lig, *this);
		VINA_FOR_IN(i, lig.pairs)
			this->update(lig.pairs[i]);
		VINA_FOR_IN(i, lig.cont)
			this->update(lig.cont[i]); // parsed_line update, below
	}
	void update(residue_vc& r) const {
		transform_ranges(r, *this);
	}
	void update(parsed_line& p) const {
		if(p.second)
			p.second = operator()(p.second.get());
	}
	void update(atom_vc& a) const {
		VINA_FOR_IN(i, a.bonds) {
			bond_vc& b = a.bonds[i];
			b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
		}
	}

	// ligands, flex, flex_context, atoms; also used for other_pairs
	template<typename T>
	void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
		sz a_sz = a.size();
		vector_append(a, b);

		is_a = true;
		VINA_FOR(i, a_sz)
			update(a[i]);

		is_a = false;
		VINA_RANGE(i, a_sz, a.size())
			update(a[i]);
	}

	// internal_coords, coords, minus_forces, atoms
	template<typename T>
	void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
		std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

		is_a = true;
		VINA_FOR_IN(i, a)
			update(a[i]);

		is_a = false;
		VINA_FOR_IN(i, b_copy)
			update(b_copy[i]);

		// interleave 
		typedef typename std::vector<T>::const_iterator const cci;
		cci b1 = b_copy.begin();
		cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
		cci b3 = b_copy.end();

		a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
		a.insert(a.end()                                , b2 , b3);
	}
};

void model::append(const model& m) {
	VINA_CHECK(atom_typing_used() == m.atom_typing_used());

	appender t(*this, m);

	t.append(other_pairs, m.other_pairs);

	VINA_FOR_IN(i, atoms)
		VINA_FOR_IN(j, m.atoms) {
			if(i >= m_num_movable_atoms && j >= m.m_num_movable_atoms) continue; // no need for inflex-inflex interactions

			const atom_vc& a =   atoms[i];
			const atom_vc& b = m.atoms[j];

			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

			if(t1 < n && t2 < n) {
				t.is_a =  true;
				sz new_i = t(i);
				t.is_a = false;
				sz new_j = t(j);
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
				other_pairs.push_back(interacting_pair(type_pair_index, new_i, new_j));
			}
		}

	VINA_CHECK(  minus_forces.size() ==   coords.size());
	VINA_CHECK(m.minus_forces.size() == m.coords.size());

	t.coords_append(internal_coords, m.internal_coords);
	t.coords_append(         coords, m         .coords);
	t.coords_append(   minus_forces, m   .minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

	t.append(ligands,         m.ligands);
	t.append(flex,            m.flex);
	t.append(flex_context,    m.flex_context);

	t       .append(grid_atoms, m.grid_atoms);
	t.coords_append(     atoms, m     .atoms);

	m_num_movable_atoms += m.m_num_movable_atoms;
}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////

atom_index model::sz_to_atom_index(sz i) const {
	if(i < grid_atoms.size()) return atom_index(i                    ,  true);
	else                      return atom_index(i - grid_atoms.size(), false);
}

distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
	if(i.in_grid && j.in_grid) return DISTANCE_FIXED;
	if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	assert(!i.in_grid);
	assert(!j.in_grid);
	assert(i.i < atoms.size());
	assert(j.i < atoms.size());
	sz a = i.i;
	sz b = j.i;
	if(a == b) return DISTANCE_FIXED;
	return (a < b) ? mobility(a, b) : mobility(b, a);
}

const vec& model::atom_coords(const atom_index& i) const {
	return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
	return vec_distance_sqr(atom_coords(a), atom_coords(b));
}

struct bond_less { // FIXME rm!?
	bool operator()(const bond_vc& a, const bond_vc& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};


bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
	fl r2 = distance_sqr_between(a, b);
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
		if(a == c || b == c) continue;
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
		if(ac != DISTANCE_VARIABLE &&
		   bc != DISTANCE_VARIABLE &&
		   distance_sqr_between(a, c) < r2 &&
		   distance_sqr_between(b, c) < r2)
			return true;
	}
	return false;
}

struct beads {
	fl radius_sqr;
	std::vector<std::pair<vec, szv> > data;
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }
	void add(sz index, const vec& coords) {
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
				data[i].second.push_back(index);
				return;
			}
		}
		// not found
		std::pair<vec, szv> tmp;
		tmp.first = coords;
		tmp.second.push_back(index);
		data.push_back(tmp);
	}
};

void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
	const fl bond_length_allowance_factor = 1.1;
	sz n = grid_atoms.size() + atoms.size();

	// construct beads
	const fl bead_radius = 15;
	beads beads_instance(n, sqr(bead_radius));
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		beads_instance.add(i, atom_coords(i_atom_index));
	}
	// assign bonds
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom_vc& i_atom = get_atom(i_atom_index);

		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;

		//find relevant atoms
		szv relevant_atoms;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
		VINA_FOR_IN(b, beads_instance.data) {
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom_vc& j_atom = get_atom(j_atom_index);
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
					//if(r2 < sqr(bond_length_allowance_factor * bond_length))
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
		// find bonded atoms
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
			if(j <= i) continue; // already considered
			atom_index j_atom_index = sz_to_atom_index(j);
			atom_vc& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
				bool rotatable = (dt == DISTANCE_ROTOR);
				fl length = std::sqrt(r2);
				i_atom.bonds.push_back(bond_vc(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond_vc(i_atom_index, length, rotatable));
			}

		}
	}
}

bool model::bonded_to_HD(const atom_vc& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond_vc& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}

bool model::bonded_to_heteroatom(const atom_vc& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond_vc& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom_vc& a = get_atom(ai);
		a.assign_el();
		sz& x = a.xs;

		bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score forumaltion apparently ignores SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

		switch(a.el) {
			case EL_TYPE_H    : break;
			case EL_TYPE_C    : x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H; break;
			case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
			case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
			case EL_TYPE_S    : x = XS_TYPE_S_P; break;
			case EL_TYPE_P    : x = XS_TYPE_P_P; break;
			case EL_TYPE_F    : x = XS_TYPE_F_H; break;
			case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;
			case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;
			case EL_TYPE_I    : x = XS_TYPE_I_H; break;
			case EL_TYPE_Met  : x = XS_TYPE_Met_D; break;
			case EL_TYPE_SIZE : break;
			default: VINA_CHECK(false);
		}
	}
}

sz model::find_ligand(sz a) const {
	VINA_FOR_IN(i, ligands)
		if(a >= ligands[i].begin && a < ligands[i].end)
			return i;
	return ligands.size();
}

void model::bonded_to(sz a, sz n, szv& out) const {
	if(!has(out, a)) { // not found
		out.push_back(a);
		if(n > 0) 
			VINA_FOR_IN(i, atoms[a].bonds) {
				const bond_vc& b = atoms[a].bonds[i];
				if(!b.connected_atom_index.in_grid)
					bonded_to(b.connected_atom_index.i, n-1, out);
			}
	}
}

szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}


void model::initialize_pairs(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, atoms) {
		sz i_lig = find_ligand(i);
		szv bonded_atoms = bonded_to(i, 3);
		VINA_RANGE(j, i+1, atoms.size()) {
			if(i >= m_num_movable_atoms && j >= m_num_movable_atoms) continue; // exclude inflex-inflex
			if(mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
				sz t1 = atoms[i].get  (atom_typing_used());
				sz t2 = atoms[j].get  (atom_typing_used());
				sz n  = num_atom_types(atom_typing_used());
				if(t1 < n && t2 < n) { // exclude, say, Hydrogens
					sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
					interacting_pair ip(type_pair_index, i, j);
					if(i_lig < ligands.size() && find_ligand(j) == i_lig)
						ligands[i_lig].pairs.push_back(ip);
					else
						other_pairs.push_back(ip);
				}
			}
		}
	}
}

void model::initialize(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, ligands)
		ligands[i].set_range();
	assign_bonds(mobility);
	assign_types();
	initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////


sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();
	return tmp;
}

szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
	sz n = num_atom_types(atom_typing_used_);
	VINA_FOR(i, m_num_movable_atoms) {
		const atom_vc& a = atoms[i];
		sz t = a.get(atom_typing_used_);
		if(t < n && !has(tmp, t))
			tmp.push_back(t);
	}
	return tmp;
}

conf_size model::get_size() const {
	conf_size tmp;
	tmp.ligands = ligands.count_torsions();
	tmp.flex    = flex   .count_torsions();
	return tmp;
}

conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
	return tmp;
}

grid_dims model::movable_atoms_box(fl add_to_each_dimension, fl granularity) const {
	vec corner1(0, 0, 0), corner2(0, 0, 0);
	VINA_FOR(i, num_movable_atoms()) {
		const vec& v = movable_coords(i);
		VINA_FOR_IN(j, v) {
			if(i == 0 || v[j] < corner1[j]) corner1[j] = v[j];
			if(i == 0 || v[j] > corner2[j]) corner2[j] = v[j];
		}
	}
	corner1 -= add_to_each_dimension / 2;
	corner2 += add_to_each_dimension / 2;

	grid_dims gd;
	{ // always doing this now FIXME ?
		vec center; center = 0.5 * (corner2 + corner1);
		VINA_FOR_IN(i, gd) {
			gd[i].n = sz(std::ceil((corner2[i] - corner1[i]) / granularity));
			fl real_span = granularity * gd[i].n;
			gd[i].begin = center[i] - real_span/2;
			gd[i].end = gd[i].begin + real_span;
		}
	}
	return gd;
}

void string_write_coord(sz i, fl x, std::string& str) {
	VINA_CHECK(i > 0);
	--i;
	std::ostringstream out;
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.setf(std::ios::showpoint);
	out << std::setw(8) << std::setprecision(3) << x;
	VINA_CHECK(out.str().size() == 8); 
	VINA_CHECK(str.size() > i + 8);
	VINA_FOR(j, 8)
		str[i+j] = out.str()[j];
}
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) { 
	std::string tmp(str);
	string_write_coord(31, coords[0], tmp);
	string_write_coord(39, coords[1], tmp);
	string_write_coord(47, coords[2], tmp);
	return tmp;
}

void model::write_context(const context& c, ofile& out) const {
	verify_bond_lengths();
	VINA_FOR_IN(i, c) {
		const std::string& str = c[i].first;
		if(c[i].second) {
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		}
		else
			out << str << '\n';
	}
}

void model::seti(const conf& c) {
	ligands.set_conf(atoms, internal_coords, c.ligands);
}

void model::sete(const conf& c) {
	VINA_FOR_IN(i, ligands)
		c.ligands[i].rigid.apply(internal_coords, coords, ligands[i].begin, ligands[i].end);
	flex.set_conf(atoms, coords, c.flex);
}

void model::set         (const conf& c) {
	ligands.set_conf(atoms, coords, c.ligands);
	flex   .set_conf(atoms, coords, c.flex);
}

fl model::gyration_radius(sz ligand_number) const {
	VINA_CHECK(ligand_number < ligands.size());
	const ligand& lig = ligands[ligand_number];
	fl acc = 0;
	unsigned counter = 0;
	VINA_RANGE(i, lig.begin, lig.end) {
		if(atoms[i].el != EL_TYPE_H) { // only heavy atoms are used
			acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
			++counter;
		}
	}
	return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


fl eval_interacting_pairs(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords) { // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
		if(r2 < cutoff_sqr) {
			fl tmp = p.eval_fast(ip.type_pair_index, r2);
			curl(tmp, v);
			e += tmp;
		}
	}
	return e;
}

fl eval_interacting_pairs_deriv(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces) { // adds to forces  // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	int count=0; 
	VINA_FOR_IN(i, pairs) {
		count++;
		const interacting_pair& ip = pairs[i];  
		vec r; r = coords[ip.b] - coords[ip.a]; // a -> b
		fl r2 = sqr(r);
		if(r2 < cutoff_sqr) {
			pr tmp = p.eval_deriv(ip.type_pair_index, r2);
			vec force; force = tmp.second * r;
			curl(tmp.first, force, v);
			e += tmp.first;
			// FIXME inefficient, if using hard curl
			forces[ip.a] -= force; // we could omit forces on inflex here
			forces[ip.b] += force;
		}
	}
	return e;
}

fl model::evali(const precalculate& p,                                  const vec& v                          ) const { // clean up
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, internal_coords); // probably might was well use coords here
	return e;
}

fl model::evale(const precalculate& p, const igrid& ig, const vec& v                          ) const { // clean up
	fl e = ig.eval(*this, v[1]);
	e += eval_interacting_pairs(p, v[2], other_pairs, coords);
	return e;
}

///// GlycoTorch Energy Equations (START)

double model::energy_2SO_1_e4_D_PHI(double angle)
{
double LH = 4.359708247316277, Lc = 29.27117168076318, LW = 3873.8975726702674, 
RH = -2.948491807983969, Rc = 126.32904338325498, RW = 826.0987631007995, 
aH = 5.848873167333682, ac = 172.7314013953158, aW = 3670.2219001991334, 
bH = 1.475893448498687, bc = 273.8161781189412, bW = 417.0692425804288, 
cH = 3.9598919873550558, cc = 310.15686340297225, cW = 897.1386567386608, 
dH = 3.7623279400842375, dc = 320.70035097807676, dW = 14299.427761495232, 
c = -0.8924033032135367,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_2SO_1_e4_D_PSI(double angle)
{
double LH = 2.5604611056949214, Lc = 303.50448722899597, LW = 516.7508796989758, 
RH = -3.2614993606521407, Rc = 303.28860490012204, RW = 931.3397339999491, 
aH = 2.820817462362836, ac = 190.3502725765731, aW = 846.8272888455898, 
bH = 1.3116652120885153, bc = 123.89981105264911, bW = 508.09132124471665, 
cH = 3.7652862325101504, cc = 61.62846391150891, cW = 952.9935834609684, 
dH = -97.80047812375938, dc = 1347.6760596440222, dW = 9553.573822154925, 
c = 0.9735579307057209,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_D_a1_4_2SO_PHI(double angle)
{
double LH = 6.9, Lc = 68.0, LW = 0.0, 
RH = 11.659152274463143, Rc = 119.79297559684207, RW = 6992.909966507498, 
aH = 1.3603650164374643, ac = 307.6251497928516, aW = 646.7882221895775, 
bH = -3.615161542405505, bc = 230.02532951385652, bW = 888.8198197109309, 
cH = -4.853804441506952, cc = 109.3070343908871, cW = 1634.8747888051992, 
dH = 2.4916833606078845, dc = 367.63640930223863, dW = 1040.4679111392395, 
c = 1.6670371052490376,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_D_a1_4_2SO_PSI(double angle)
{
double LH = 6.94871708950802, Lc = 54.578177790579005, LW = 1535.2698918725655, 
RH = 5.554362074751906, Rc = 194.3413033329821, RW = 3887.2784380306507, 
aH = 0.5602230089917873, ac = 299.9126055955525, aW = 211.4877981092501, 
bH = -4.369144200221367, bc = 226.70745099673786, bW = 1417.0122348174993, 
cH = -2.34661024226204, cc = 43.773110495046964, cW = 481.8675745806132, 
dH = 1.7914817393081224, dc = 363.9784183253017, dW = 182.84455394321645, 
c = 0.4840111578002676,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_L_a1_e4_D_PSI(double angle)
{
double LH = 12.962429302892675, Lc = 85.46229703521264, LW = 1458.7892164485636, 
RH = 3.1603145455476094, Rc = 191.46296594939247, RW = 1279.7830283754727, 
aH = -0.6246596560505729, ac = 305.87064526089017, aW = 2306.553297433773, 
bH = -0.40618765516122574, bc = 221.27228560794836, bW = 273.88862409591474, 
cH = -10.853485895524086, cc = 89.70184526684163, cW = 914.4743175542117, 
dH = 1.0, dc = 0.0, dW = 0.0, 
c = 0.6436456386363167,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_L_a1_e4_D_PHI(double angle)
{
double LH = 7.614592811265287, Lc = 289.67204203945215, LW = 3978.885414576388, 
RH = 5.388137729917904, Rc = 208.0520242319089, RW = 1898.5208167527662, 
aH = 4.279392016058717, ac = 179.69816073860844, aW = 971.5517628865543, 
bH = 0.5692306281653353, bc = 1988.5455350419015, bW = 1.121416106131602, 
cH = 2.156369999432613, cc = 61.77140760350418, cW = 985.5875515527207, 
dH = 2.9182925365463226, dc = 4.20233633251493, dW = 2252.5727175010406, 
c = -0.22526059365739395,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_D_1a_4a_L_PHI(double angle)
{
double LH = 8.215355663608781, Lc = 81.58142683598575, LW = 5242.4277869765265, 
RH = 6.0686645297686646, Rc = 157.71400487058838, RW = 1068.5912001366532, 
aH = 3.50984970930512, ac = 187.1650289881355, aW = 531.7794783561717, 
bH = 1.3468040262930794, bc = 279.4932995964358, bW = 256.18420542249663, 
cH = 2.1342418228040256, cc = 308.1583831414212, cW = 480.3928637353839, 
dH = 3.279983501237371, dc = 348.73877215579176, dW = 705.6382095328541, 
c = -0.02742111038999487,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_D_1a_4a_L_PSI(double angle)
{
double LH = 4.038572351124325, Lc = 288.1223228606298, LW = 1090.2320324585273, 
RH = 5.255277821355811, Rc = 237.0330034683259, RW = 746.534192812652, 
aH = 5.854973867113257, ac = 178.7272547545124, aW = 1291.908751246364, 
bH = 1.0899910315417805, bc = 104.95199123398154, bW = 110.7492576716796, 
cH = 0.5846991735544921, cc = 56.389700673389164, cW = 262.17912222577365, 
dH = 0.9112707352083366, dc = 4.035111186433551, dW = 28.681341540120915, 
c = 0.026484054464878493,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_4c1_a1_a4_4c1_PHI(double angle)
{
double LH = 8.004543186726876, Lc = 77.85233712724319, LW = 5455.960123111793, 
RH = 6.27840620609398, Rc = 155.5752083277422, RW = 1246.6074492935256, 
aH = 4.339147785143784, ac = 185.55908530905037, aW = 702.7872247711337, 
bH = 3.086522647595491, bc = 329.66839919445357, bW = 4730.334624920652, 
cH = 1.2326055148176038, cc = 356.4232838489615, cW = 264.75091371997235, 
dH = 6.647856428704409, dc = 385.98945416549094, dW = 30.341388740661834, 
c = -0.6689872272627778,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_4c1_a1_a4_4c1_PSI(double angle)
{
double LH = 3.409413378028371, Lc = 41.047256725062795, LW = 438.4188036833962, 
RH = 5.2074157902203435, Rc = 93.7610143227761, RW = 3478.619983211818, 
aH = 5.645865478301781, ac = 172.78854477370072, aW = 1099.3424881797048, 
bH = 0.5392939612690947, bc = 288.6282067422295, bW = 309.9954681996438, 
cH = 0.19777, cc = 3080.199, cW = 224.193, 
dH = 38.656778469013524, dc = 1104.964659112571, dW = 148373.97759200912, 
c = -0.4038347511817272,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_4c1_b1_e3_4c1_PHI(double angle)
{
double LH = 6.9, Lc = 68.0, LW = 0.0, 
RH = 4.32711277425824, Rc = 181.46914889527443, RW = 1235.9818531573799, 
aH = 2.8886672952308095, ac = 298.7912187800737, aW = 1280.316237492101, 
bH = -1.298180460482828, bc = 250.30335666232517, bW = 620.9712745436813, 
cH = -2.576379105563004, cc = 109.11095696281704, cW = 3043.200209080894, 
dH = -2.7057312256679777, dc = 354.16793044412776, dW = 38.45515297600506, 
c = 2.591561729544288,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_4c1_b1_e3_4c1_PSI(double angle)
{
double LH = -16.500350321423088, Lc = 81.04613249850134, LW = 3324.3093061958643, 
RH = 0.19332886656597498, Rc = -16.42785995686048, RW = 21.561100631640812, 
aH = 4.143131612469768, ac = 188.99388352199406, aW = 1285.2853933906251, 
bH = 2.6041941513620297, bc = 257.095292985663, bW = 501.9166686393258, 
cH = 4.134668202627348, cc = 311.2074873861674, cW = 1238.1800052376595, 
dH = 16.940905465703946, dc = 80.10601851046633, dW = 3642.862864087344, 
c = -0.2636164402052412,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

//

double model::energy_4c1_b1_4e_4c1_PHI(double angle)
{
double LH = 9.072264959383777, Lc = 13.903369029602617, LW = 19585.685526697027, 
RH = 1.0260158377303492, Rc = 126.97681191469557, RW = 350.9840047498267, 
aH = 4.552720962681688, ac = 157.74530574591387, aW = 824.8319296094436, 
bH = 7.738982922572851, bc = 194.97612776782356, bW = 1650.6175708315393, 
cH = 1.8460097576649908, cc = 290.44189995703096, cW = 974.8782280277637, 
dH = 10.407613696290978, dc = 315.6148749908362, dW = 8213.1086340685, 
c = -6.561879366631349,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

double model::energy_4c1_b1_4e_4c1_PSI(double angle)
{
double LH = 3.826247303166274, Lc = 62.65126784061971, LW = 974.9193945663436, 
RH = 1.8846986460682171, Rc = 130.5623697724467, RW = 447.57435633353515, 
aH = 3.289123567470172, ac = 189.4577579387754, aW = 1179.3346979215887, 
bH = 0.6635688943301651, bc = 268.13876288043014, bW = 416.7488587408713, 
cH = -2.489011981034135, cc = 620.2069640808243, cW = 633922.8496762958, 
dH = 1.338345311946964, dc = 362.5544913580919, dW = 151.28553574836548, 
c = 2.270729155804515,
Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + c;
return Totx;
}

///// GlycoTorch Energy Equations (END)


double model::phi_alpha_energy(double phi_angle)
{
double LH = 2.97696467271672, Lc = -199.494365163839, LW = 677.808323900125, RH = 102.253303636096, Rc = 170.599580473404, RW = 1696.78443699429, aH = 10.7448005875571, ac = -105.313553566706, aW = 4724.58364072706, bH = 3.67344580413578, bc = 6.20116671874232, bW = 1347.72056251564, cH = 2.06094652655659, cc = 91.6553021324274, cW = 1500.02002601097, Off = 1.00501e-30, dH = 6.19388683252667, dc = -22.9786969888816, dW = 2122.27783139301, eH = -2.11153017593601, ec = 83.6019123356148, eW = 1254.13371108961, fH = -98.0013005657107, fc = 170.012289132741, fW = 1598.73272567307, Leftx, Rightx, ax, bx, cx, dx, ex, fx, x, Totx;
x=phi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + fx;
return Totx;
}

double model::phi_beta_energy(double phi_angle)
{
double Lc = -330.769995527134, aH = 5.93533323829663, ac = -152.080139620062, aW = 6049.77220005964, bH = 22.467372096061, bc = -23.5159916173247, bW = 606.89715970453, cH = 10.0360057033439, cc = 120.962836525241, cW = 4037.89330459304, LH = 450.540038600828, LW = 4449.7622241787, RH = 23.7118506901333, Rc = 304.624980492529, RW = 8375.1929028027, /*Off = -2.27829251796721,*/ Off = -2.1283, dH = -18.1406478565247, dc = -24.2677756921736, dW = 543.050986049266, eH = 5.88226333077368, ec = 19.6321032903376, eW = 897.92664572344, Leftx, Rightx, ax, bx, cx, dx, ex, x, Totx;
x=phi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
Totx = Rightx + Leftx + ax + bx + cx + dx + ex + Off;
return Totx;
}

double model::psi_2A3E_energy(double psi_angle)
{
double LH = 4.62366277694845, Lc = 5.045583934308, LW = 5005.75236060956, RH = 4.61387450239844, Rc = 362.487847702007, RW = 2090.63190217702, aH = 4.94191727813274, ac = 121.202321824468, aW = 2093.75214491931, bH = 0.402901504292045, bc = 241.428583877882, bW = 456.828754790442, cH = 0.798876573705798, cc = 68.425080241155, cW = 678.807178379645, Off = -0.125645118474882, dc = 192.925748017071, dW = 347.244734136509, dH = 0.222992242737354, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
if(psi_angle<0)
{
psi_angle=360+psi_angle;}
x=psi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
return Totx;
}

double model::psi_2E3A_energy(double psi_angle)
{
double LH = 4.46811874171788, Lc = 1e-30, LW = 1279.58772056199, RH = 4.38204018882823, Rc = 357.770654336205, RW = 6050.14162479438, aH = 284.944948778136, ac = 146.644068129462, aW = 1551.75673776163, bH = 4.76134025362478, bc = 220.683118921686, bW = 5892.94143218231, cH = -169.197666368856, cc = 147.370828680332, cW = 1742.47541063603, Off = 1.0219924486158, dc = 146.05660843428, dW = 1359.82873591396, dH = -118.440552792375, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
if(psi_angle<0){
psi_angle=360+psi_angle;}
x=psi_angle;
Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
return Totx;
}

double model::psi_6A_energy(double psi_angle)
{
double aH = 67.9431348410598, ac = -59.5393395706705, aW = 993.323581145538, bH = 6.13421142432396, bc = 10.4786088782815, bW = 945.770771330812, cH = 3.27628727235978, cc = 54.2960678151208, cW = 851.528141801851, dH = 0.727486729062442, dc = 131.067737803489, dW = 1037.41211378392, eH = 2.57362265878937, ec = 245.102891425541, eW = 2012.99451568206, fH = 5.75995973448166, fc = 359.999988549478, fW = 1153.3974275673, gH = 3.47492643928157, gc = 321.677942414686, gW = 2080.97053159226, hH = -0.741000462200939, hc = 199.106903524814, hW = 522.180434119001, ax, bx, cx, dx, ex, fx, gx, hx, x, Totx;

if(psi_angle<0)
{
psi_angle=360+psi_angle;
}
x=psi_angle;
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
gx = gH * exp(-pow((x-(gc)),2.0)/gW);
hx = hH * exp(-pow((x-(hc)),2.0)/hW);
Totx = ax + bx + cx + dx + ex + fx + gx + hx;
return Totx;
}

double model::psi_6E_energy(double psi_angle)
{
double aH = 7.24858655753829, ac = 3.60600554520403, aW = 2459.23916629141, bH = 1.9, bc = 96.5930821702371, bW = 2683.88656516991, cH = 0.741022592342903, cc = 141.663521919709, cW = 1150.04756181103, dH = 0.2, dc = 162, dW = 400, eH = 0.287090039932611, ec = 228.171314273305, eW = 272.201363844744, fH = 1.22591971967808, fc = 292.206221787048, fW = 1134.52455512381, gH = 7.41063235334191, gc = 369.867701147817, gW = 3499.15994772992, hH = -0.61489499584011, hc = 271.319024293053, hW = 532.437194483944, iH = -0.35, ic = 183, iW = 100, ax, bx, cx, dx, ex, fx, gx, hx, ix, x, Totx;

if(psi_angle<0)
{
psi_angle=360+psi_angle;
}
x=psi_angle;
ax = aH * exp(-pow((x-(ac)),2.0)/aW);
bx = bH * exp(-pow((x-(bc)),2.0)/bW);
cx = cH * exp(-pow((x-(cc)),2.0)/cW);
dx = dH * exp(-pow((x-(dc)),2.0)/dW);
ex = eH * exp(-pow((x-(ec)),2.0)/eW);
fx = fH * exp(-pow((x-(fc)),2.0)/fW);
gx = gH * exp(-pow((x-(gc)),2.0)/gW);
hx = hH * exp(-pow((x-(hc)),2.0)/hW);
ix = iH * exp(-pow((x-(ic)),2.0)/iW);
Totx = ax + bx + cx + dx + ex + fx + gx + hx + ix;
return Totx;
}

double model::omega_6A_energy(double omega_angle)
{
double x, energy, b, k=0.0025;
        if(omega_angle<0)
        {
        x=360+omega_angle;
        }
        else
        {
        x=omega_angle;
        }
        if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
        {
        b=0.0;
        energy=k*pow((x-60),2)+b;
        }
        else if(x>=120.0 && x<240.0)
        {
        b=0.3;
        energy=k*pow((x-180),2)+b;
        }
        else if(x>=240.0 && x<360.0)
        {
        b=1.0;
        energy=k*pow((x-300),2)+b;
        }
return energy;
}


double model::omega_6E_energy(double omega_angle)
{
double x, energy, b, k=0.0025;
        if(omega_angle<0)
        {
        x=360+omega_angle;
        }
        else
        {
        x=omega_angle;
        }
        if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
        {
        b=0.21;
        energy=k*pow((x-60),2)+b;
        }
        else if(x>=120.0 && x<240.0)
        {
        b=1.39;
        energy=k*pow((x-180),2)+b;
        }
        else if(x>=240.0 && x<360.0)
        {
        b=0.0;
        energy=k*pow((x-300),2)+b;
        }
return energy;
}



fl model::eval         (const precalculate& p, const igrid& ig, const vec& v, const conf& c           ) { // clean up
	set(c);
	fl e = evale(p, ig, v);
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords
	return e;
}

double model::get_torsion_coords_vecs_list(vec A, vec B, vec C, vec D)
{
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

vecv model::get_flexible_coords(){
return coords;
}


fl model::eval_chi(const fl chi_coeff, const fl chi_cutoff)
{
fl e=0.0;
if(chi_coeff!=0){
	double phi_energy=0.0, psi_energy=0.0,omega_energy=0.0, phi=0.0, psi=0.0, omega=0.0, total=0.0,current_energy=0.0;
	std::vector< std::vector<size_t*> > glyco_info=glycan_info_func();
	VINA_FOR(i,glyco_info.size())
	{
	phi=get_torsion_coords_vecs_list(coords[glyco_info[i][0][0]],coords[glyco_info[i][1][0]],coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]]);
	psi=get_torsion_coords_vecs_list(coords[glyco_info[i][1][0]],coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]],coords[glyco_info[i][4][0]]);



	/// GlycoTorch GAG specific scoring functions (START)

	// While there is some cross-over with the two SFs, having the GlycoTorch SFs called first in the IF/ELSE loop
	// means there will be no double-ups. 

	///                       2SO                        4C1	     equitorial = 1
	if ((glyco_info[i][10][0]==3) && (glyco_info[i][11][0]==1) && (glyco_info[i][6][0]==1)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_2SO_1_e4_D_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_2SO_1_e4_D_PSI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}

	///                       4C1                        2SO				alpha = 0              equitorial = 1
	if ((glyco_info[i][10][0]==1) && (glyco_info[i][11][0]==3) && (glyco_info[i][5][0]==0) && (glyco_info[i][6][0]==1)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_D_a1_4_2SO_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_D_a1_4_2SO_PSI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}

	///                       1C4                        4C1				alpha = 0              equitorial = 1
	if ((glyco_info[i][10][0]==2) && (glyco_info[i][11][0]==1) && (glyco_info[i][5][0]==0) && (glyco_info[i][6][0]==1)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_L_a1_e4_D_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_L_a1_e4_D_PHI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}


	///                       4C1                       1C4				alpha = 0              axial = 0
	if ((glyco_info[i][10][0]==1) && (glyco_info[i][11][0]==2) && (glyco_info[i][5][0]==0) && (glyco_info[i][6][0]==0)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_L_a1_e4_D_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_L_a1_e4_D_PHI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}


	///                       4C1                       4C1				alpha = 0              axial = 0
	if ((glyco_info[i][10][0]==1) && (glyco_info[i][11][0]==1) && (glyco_info[i][5][0]==0) && (glyco_info[i][6][0]==0)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_4c1_a1_a4_4c1_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_4c1_a1_a4_4c1_PSI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}

	///                       4C1                       4C1				beta = 1              equitorial = 1
	if ((glyco_info[i][10][0]==1) && (glyco_info[i][11][0]==1) && (glyco_info[i][5][0]==1) && (glyco_info[i][6][0]==1)
	
	//                    linkage
	&& (glyco_info[i][7][0]==3) ) {

				current_energy=energy_4c1_b1_e3_4c1_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_4c1_b1_e3_4c1_PSI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}

	///                       4C1                       4C1				beta = 1              equitorial = 1
	if ((glyco_info[i][10][0]==1) && (glyco_info[i][11][0]==1) && (glyco_info[i][5][0]==1) && (glyco_info[i][6][0]==1)
	
	//                    linkage
	&& (glyco_info[i][7][0]==4) ) {

				current_energy=energy_4c1_b1_4e_4c1_PHI(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}

				current_energy=energy_4c1_b1_4e_4c1_PSI(psi);
				if(current_energy>chi_cutoff){
				psi_energy+=current_energy;
				current_energy=0.0;
				}
	}
	/// GlycoTorch GAG specific scoring functions (END)

	else {

	if(glyco_info[i][10][0]!=0){
		if(glyco_info[i][5][0]==0)
		{//Alpha
			if(glyco_info[i][10][0]==1)
			{//D sugar (Alpha)
			current_energy=phi_alpha_energy(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
			else if(glyco_info[i][10][0]==2)
			{//L-Sugar (Alpha)
			current_energy=phi_alpha_energy(-phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
                current_energy=0.0;
				}
			}
		}
		else
		{//Beta
			if(glyco_info[i][10][0]==1)
			{//D sugar (Beta)
			current_energy=phi_beta_energy(phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
			else if(glyco_info[i][10][0]==2)
			{//L Sugar (Beta)
			current_energy=phi_beta_energy(-phi);
				if(current_energy>chi_cutoff){
				phi_energy+=current_energy;
				current_energy=0.0;
				}
			}
		}
	}
	if(glyco_info[i][11][0]!=0){
	omega=get_torsion_coords_vecs_list(coords[glyco_info[i][2][0]],coords[glyco_info[i][3][0]],coords[glyco_info[i][4][0]],coords[glyco_info[i][9][0]]);
		if(glyco_info[i][7][0]==2 || glyco_info[i][7][0]==4)
		{//2 or 4 linkage
			if(glyco_info[i][6][0]==0)
			{//Axial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2A3E_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
					current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2A3E_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			}
			else
			{//Equatorial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2E3A_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
					current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2E3A_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			
			}
		}
		else if(glyco_info[i][7][0]==3)
		{//3 linkage
			if(glyco_info[i][6][0]==0)
			{//Axial attachment
			
				if(glyco_info[i][11][0]==2)
				{//L sugar
				current_energy=psi_2E3A_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
                                current_energy=psi_2E3A_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			}
			else
			{//Equatorial attachment
				if(glyco_info[i][11][0]==2)
				{//L sugar
                                current_energy=psi_2A3E_energy(-psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
				else if(glyco_info[i][11][0]==1)
				{//D sugar
				current_energy=psi_2A3E_energy(psi);
					if(current_energy>chi_cutoff){
					psi_energy+=current_energy;
                                        current_energy=0.0;
					}
				}
			
			}
		}
	else if(glyco_info[i][7][0]==6)
                {//6 linkage
                        if(glyco_info[i][8][0]==0)
                                {
                                current_energy=psi_6A_energy(psi);
                                        if(current_energy>chi_cutoff){
                                        psi_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                current_energy=omega_6A_energy(omega);
                                        if(current_energy>chi_cutoff){
                                        omega_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                }
                        else if(glyco_info[i][8][0]==1)
                                {
                                current_energy=psi_6E_energy(psi);
                                        if(current_energy>chi_cutoff){
                                        psi_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                current_energy=omega_6E_energy(omega);
                                        if(current_energy>chi_cutoff){
                                        omega_energy+=current_energy;
                                        current_energy=0.0;
                                        }
                                }
                }
	}
}
}
	total=phi_energy+psi_energy+omega_energy; total=total*chi_coeff; double combined_score=total+e; e=combined_score;}
	else
	e=0.0;
// std::cout << e << " ";
return e;
}


fl  model::eval_deriv  (const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g/*, std::vector< std::vector<int> > glyco_info*/,const fl chi_coeff, const fl chi_cutoff ) { // clean up
	set(c);
	fl chi_energy = 0.0;
	chi_energy=eval_chi(chi_coeff,chi_cutoff);
	fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex
        e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
        VINA_FOR_IN(i, ligands){
                e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces
}
        ligands.derivative(coords, minus_forces, g.ligands);
	flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored
	return (e+chi_energy);
}


fl model::eval_intramolecular(const precalculate& p, const vec& v, const conf& c) {
	set(c);
	fl e = 0;

	// internal for each ligand
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

	sz nat = num_atom_types(atom_typing_used());
	const fl cutoff_sqr = p.cutoff_sqr();

	// flex-rigid
	VINA_FOR(i, num_movable_atoms()) {
		if(find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
		const atom_vc& a = atoms[i];
		sz t1 = a.get(atom_typing_used());
		if(t1 >= nat) continue;
		VINA_FOR_IN(j, grid_atoms) {
			const atom_vc& b = grid_atoms[j];
			sz t2 = b.get(atom_typing_used());
			if(t2 >= nat) continue;
			fl r2 = vec_distance_sqr(coords[i], b.coords);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = triangular_matrix_index_permissive(nat, t1, t2);
				fl this_e = p.eval_fast(type_pair_index, r2);
				curl(this_e, v[1]);
				e += this_e;
			}
		}
	}

	// flex-flex
	VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& pair = other_pairs[i];
		if(find_ligand(pair.a) < ligands.size() || find_ligand(pair.b) < ligands.size()) continue; // we only need flex-flex
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		if(r2 < cutoff_sqr) {
			fl this_e = p.eval_fast(pair.type_pair_index, r2);
			curl(this_e, v[2]);
			e += this_e;
		}
	}
	return e;
}


fl model::eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy) {
	fl e = eval(p, ig, v, c); // sets c
	return sf.conf_independent(*this, e - intramolecular_energy);
}

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
	sz n = x.m_num_movable_atoms; 
	VINA_CHECK(n == y.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, n) {
		const atom_vc& a =   x.atoms[i];
		if(a.el != EL_TYPE_H) {
			fl r2 = max_fl;
			VINA_FOR(j, n) {
				const atom_vc& b = y.atoms[j];
				if(a.same_element(b) && !b.is_hydrogen()) {
					fl this_r2 = vec_distance_sqr(x.coords[i], 
					                              y.coords[j]);
					if(this_r2 < r2)
						r2 = this_r2;
				}
			}
			assert(not_max(r2));
			sum += r2;
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_lower_bound(const model& m) const {
	return (std::max)(rmsd_lower_bound_asymmetric(*this,     m),
		            rmsd_lower_bound_asymmetric(    m, *this));
}

fl model::rmsd_upper_bound(const model& m) const {
	VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, m_num_movable_atoms) {
		const atom_vc& a =   atoms[i];
		const atom_vc& b = m.atoms[i];
		assert(a.ad == b.ad);
		assert(a.xs == b.xs);
		if(a.el != EL_TYPE_H) {
			sum += vec_distance_sqr(coords[i], m.coords[i]);
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_ligands_upper_bound(const model& m) const {
	VINA_CHECK(ligands.size() == m.ligands.size());
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR_IN(ligand_i, ligands) {
		const ligand&   lig =   ligands[ligand_i];
		const ligand& m_lig = m.ligands[ligand_i];
		VINA_CHECK(lig.begin == m_lig.begin);
		VINA_CHECK(lig.end   == m_lig.end);
		VINA_RANGE(i, lig.begin, lig.end) {
			const atom_vc& a =   atoms[i];
			const atom_vc& b = m.atoms[i];
			assert(a.ad == b.ad);
			assert(a.xs == b.xs);
			if(a.el != EL_TYPE_H) {
				sum += vec_distance_sqr(coords[i], m.coords[i]);
				++counter;
			}
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}


void model::verify_bond_lengths() const {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		const atom_vc& a = get_atom(ai);
		VINA_FOR_IN(j, a.bonds) {
			const bond_vc& b = a.bonds[j];
			fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
			bool ok = eq(d, b.length);
			if(!ok) {
				VINA_SHOW(d);
				VINA_SHOW(b.length);
			}
			VINA_CHECK(ok);
		}
	}
}

void model::check_internal_pairs() const {
	VINA_FOR_IN(i, ligands) {
		const ligand& lig = ligands[i];
		const interacting_pairs& pairs = lig.pairs;
		VINA_FOR_IN(j, pairs) {
			const interacting_pair& ip = pairs[j];
			VINA_CHECK(ip.a >= lig.begin);
			VINA_CHECK(ip.b  < lig.end);
		}
	}
}

void model::about() const {
	VINA_SHOW(atom_typing_used());
	VINA_SHOW(num_movable_atoms());
	VINA_SHOW(num_internal_pairs());
	VINA_SHOW(num_other_pairs());
	VINA_SHOW(num_ligands());
	VINA_SHOW(num_flex());
}

void model::print_stuff() const {
	std::cout << "coords:\n";
	VINA_FOR_IN(i, coords)
		printnl(coords[i]);

	std::cout << "internal_coords:\n";
	VINA_FOR_IN(i, internal_coords)
		printnl(internal_coords[i]);

	std::cout << "atoms:\n";
	VINA_FOR_IN(i, atoms) {
		const atom_vc& a = atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}

	std::cout << "grid_atoms:\n";
	VINA_FOR_IN(i, grid_atoms) {
		const atom_vc& a = grid_atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}
	about();
}

fl pairwise_clash_penalty(fl r, fl covalent_r) {
	// r = 0          -> max_penalty 
	// r = covalent_r -> 1
	// elsewhere      -> hyperbolic function
	assert(r >= 0);
	assert(covalent_r > epsilon_fl);
	const fl x = r / covalent_r;
	if(x > 2) return 0;
	return 1-x*x/4;
}

fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
		const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
		e += pairwise_clash_penalty(r, covalent_r);
	}
	return e;
}

fl model::clash_penalty() const {
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += clash_penalty_aux(ligands[i].pairs);
	e += clash_penalty_aux(other_pairs);
	return e;}
