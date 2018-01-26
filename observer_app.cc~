// (c) Copyright Rosetta Commons Member Institutions.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <random>

#include <numeric/random/random.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/residue_selector/ClashBasedRepackShellSelector.hh>
#include <core/pack/task/operation/ReplicateTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

using namespace core::scoring;

int main( int argc, char ** argv ) {
    std::cout << "Hello World!" << std::endl;

    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value(); //enter PDB filename here
    if ( filenames.size() > 0 ) {
        std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
    
    ScoreFunctionOP sfxn = get_score_function();
    core::Real score = sfxn->score( *mypose );

    std::cout << "The score is: " << score << std::endl;

    float N = mypose->size(); //calls 'size' method of 'mypose' object; max value for loop
    float it = 0; //counter for loop
    protocols::moves::MonteCarlo::Real temp = 1;

    protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo( *mypose, *sfxn, temp );

    protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );
    the_observer->pymol().apply( *mypose);
    //bool sanity = mypose->residue(randres).is_protein()

    core::kinematics::MoveMap mm; //set backbone/chi-angle values
    mm.set_bb( true );
    mm.set_chi( true );

    core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
    core::optimization::AtomTreeMinimizer atm;

    core::pose::Pose copy_pose = *mypose;

    while (it < N ) {
        double uniform_random_number = ((double) rand() / (RAND_MAX));
        core::Size randres = static_cast< core::Size > ( uniform_random_number * N + 1 );//… code here to pick the index of a random residue in the Pose

        core::Real pert1 = static_cast< core::Size > (numeric::random::gaussian());//… code here to get a random number
        core::Real pert2 = static_cast< core::Size > (numeric::random::gaussian());//… code here to get another random number
        core::Real orig_phi = mypose->phi( randres );
        core::Real orig_psi = mypose->psi( randres );
        mypose->set_phi( randres, orig_phi + pert1 );
        mypose->set_psi( randres, orig_psi + pert2 );

        core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
        repack_task->restrict_to_repacking();
        core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
        atm.run( copy_pose, mm, *sfxn, min_opts );
        *mypose = copy_pose;

        mc.boltzmann( *mypose );
        it++;
    }
    std::cout << "/n Exiting loop now..." << std::endl;
    return 0;
}
