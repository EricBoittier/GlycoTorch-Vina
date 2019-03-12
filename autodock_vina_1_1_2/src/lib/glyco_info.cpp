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

struct glyco_info
{
	size_t *S1_C1, *S1_O5, *S1_O1, *S1_C5, *S1_C2, *S1_C3, *S1_C4, *S2_Cxp1, *S2_C2, *S2_O5, *S2_C5, *S2_C1, *S2_O1, *S2_Ox, *S2_Cx, *S2_Cxm1, *S2_O4, *S2_C4, *S2_C3, *sizet_ring1_conf;
	size_t *S1_AB, *S2_AE, *S2_Link, *S2_6AE, *sizet_ring2_conf; 

};


		// glyco_info.push_back(S1_O5); //index 0
		// glyco_info.push_back(S1_C1); //index 1
		// glyco_info.push_back(S2_Ox); //index 2
		// glyco_info.push_back(S2_Cx); //index 3
		// glyco_info.push_back(S2_Cxm1); //index 4
		// //Co-ordinates END
		// glyco_info.push_back(S1_AB); //index 5
		// glyco_info.push_back(S2_AE); //index 6
		// glyco_info.push_back(S2_Link); //index 7 
        // glyco_info.push_back(S2_6AE); //index 8  //0 -> positive; 1-> negative
        // //Co-ordinates START
    	// glyco_info.push_back(S2_O5); //index 9 
		// glyco_info.push_back(sizet_ring1_conf); //index 10
		// glyco_info.push_back(sizet_ring2_conf); //index 11
        // //Co-ordinates ED
		// ligand_glyco_info.push_back(glyco_info);