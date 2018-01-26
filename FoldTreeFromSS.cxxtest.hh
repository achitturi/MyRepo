// (c) Copyright Rosetta Commons Member Institutions.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>


utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
    utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
    int strand_start = -1;
    for ( core::Size j = 0; j < ss_string.size(); ++j ) {
        if ( ss_string[ j ] == 'E' || ss_string[ j ] == 'H'  ) {
            if ( strand_start == -1 ) {
                strand_start = j;
            } else if ( ss_string[j] != ss_string[strand_start] ) {
                ss_boundaries.push_back( std::make_pair( strand_start+1, j ) );
                strand_start = j;
            } 
        }   else {
                if ( strand_start != -1 ) {
                    ss_boundaries.push_back( std::make_pair( strand_start+1, j ) );
                    strand_start = -1;
                }
            }
        }
    if ( strand_start != -1 ) {
        // last residue was part of a ss-eleemnt
        ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
    }    
    return ss_boundaries;
}

core::kinematics::FoldTree fold_tree_from_dssp_string(std::string dssp_string) {
    core::kinematics::FoldTree ft;

    utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
    ss_boundaries = identify_secondary_structure_spans( dssp_string );
    
    core::Size residue_orig;
    core::Size loop_residue_orig;
    core::Size first_orig;
    core::Size jump_counter = 1;
    for ( core::Size j = 1; j <= ss_boundaries.size(); ++j ) {
        residue_orig = static_cast< core::Size >(((ss_boundaries[j].first + ss_boundaries[j].second)/2));
        std::cout << residue_orig << std::endl;
        if ( j != ss_boundaries.size() ) {
            loop_residue_orig = static_cast<core::Size> (((ss_boundaries[j].second + ss_boundaries[j+1].first)/2));
            if ( j == 1 ) {
                first_orig = residue_orig;
                ft.add_edge(first_orig, 1, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(first_orig, ss_boundaries[j].second, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(first_orig, loop_residue_orig, jump_counter);
                ++jump_counter;
                ft.add_edge(loop_residue_orig, ++ss_boundaries[j].second, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(loop_residue_orig, --ss_boundaries[j + 1].first, core::kinematics::Edge::PEPTIDE);
            } else {
                ft.add_edge(first_orig, residue_orig, jump_counter);
                ++jump_counter;
                ft.add_edge(first_orig, loop_residue_orig, jump_counter);
                ++jump_counter;
                ft.add_edge(residue_orig, ss_boundaries[j].first, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(residue_orig, ss_boundaries[j].second, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(loop_residue_orig, ++ss_boundaries[j].second, core::kinematics::Edge::PEPTIDE);
                ft.add_edge(loop_residue_orig, --ss_boundaries[j+1].first, core::kinematics::Edge::PEPTIDE);
            }
        } else {
            ft.add_edge(first_orig, residue_orig, jump_counter);
            ft.add_edge(residue_orig, ss_boundaries[j].first, core::kinematics::Edge::PEPTIDE);
            ft.add_edge(residue_orig, dssp_string.size(), core::kinematics::Edge::PEPTIDE); // to end
        }
    }

    return ft;
}

core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose const & mypose) {
    core::scoring::dssp::Dssp dssp_( mypose );
    return fold_tree_from_dssp_string( dssp_.get_dssp_secstruct() );
}


//  Test Class  //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	//  Fixtures  //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}


	//  Test Cases  //

	void test_string1() {
        utility::vector1< std::pair< core::Size, core::Size > > testing = identify_secondary_structure_spans("   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ");
        TS_ASSERT( testing[1].first == 4 && testing[1].second == 8 );
        TS_ASSERT( testing[2].first == 12 && testing[2].second == 19 );
        TS_ASSERT( testing[3].first == 22 && testing[3].second == 26 );
        TS_ASSERT( testing[4].first == 36 && testing[4].second == 41 );
        TS_ASSERT( testing[5].first == 45 && testing[5].second == 55 );
        TS_ASSERT( testing[6].first == 58 && testing[6].second == 62 );
        TS_ASSERT( testing[7].first == 65 && testing[7].second == 68 );
    }

	void test_string2() {
        utility::vector1< std::pair< core::Size, core::Size > > testing = identify_secondary_structure_spans("HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ");
        TS_ASSERT( testing[1].first == 1 && testing[1].second == 7 );
        TS_ASSERT( testing[2].first == 11 && testing[2].second == 22 );
        TS_ASSERT( testing[3].first == 29 && testing[3].second == 40 );
        TS_ASSERT( testing[4].first == 41 && testing[4].second == 50 );
        TS_ASSERT( testing[5].first == 51 && testing[5].second == 57 );
        TS_ASSERT( testing[6].first == 59 && testing[6].second == 62 );
        TS_ASSERT( testing[7].first == 63 && testing[7].second == 65 );
    }

	void test_string3() {
        utility::vector1< std::pair< core::Size, core::Size > > testing = identify_secondary_structure_spans("EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE");
        TS_ASSERT( testing[1].first == 1 && testing[1].second == 9 );
        TS_ASSERT( testing[2].first == 11 && testing[2].second == 18 );
        TS_ASSERT( testing[3].first == 20 && testing[3].second == 28 );
        TS_ASSERT( testing[4].first == 30 && testing[4].second == 30 );
        TS_ASSERT( testing[5].first == 32 && testing[5].second == 36 );
        TS_ASSERT( testing[6].first == 38 && testing[6].second == 38 );
        TS_ASSERT( testing[7].first == 40 && testing[7].second == 40 );
        TS_ASSERT( testing[8].first == 42 && testing[8].second == 42 );
        TS_ASSERT( testing[9].first == 44 && testing[9].second == 51 );
    }

    void test_output() {
        std::cout << fold_tree_from_dssp_string("   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ");
        TS_ASSERT( true);
    }

};
