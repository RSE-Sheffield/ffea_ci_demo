#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cctype>
#include <omp.h>

#include <boost/program_options.hpp>
#include "FFEA_return_codes.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
//#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#define MAX_FNAME_LENGTH 255

using namespace std;
namespace b_po = boost::program_options;
namespace b_fs = boost::filesystem; 

int main(int argc, char *argv[])
{
#ifdef USE_MPI
  MPI::Init();
#endif
	cout << "\n\n\n***************************************************\n\tFLUCTUATING FINITE ELEMENT ANALYSIS\n***************************************************\n\n" << endl;
	cout << " Version:\t" << FFEA_VERSION << " [" << FFEA_MASCOT << "]" << endl;
	cout << "Compiled:\t" << __DATE__ " at " << __TIME__ << endl;
	cout << "  Coding:\tAlbert Solernou (a.solernou@leeds.ac.uk), Ben Hanson (py09bh@leeds.ac.uk), Robin Richardson (pyrar@leeds.ac.uk),\n" << endl;
	cout << "  Theory:\tOliver Harlen, Sarah Harris, Robin Oliver, Daniel Read, Robin Richardson, Ben Hanson, Albert Solernou\n" << endl;

	#ifdef FFEA_PARALLEL_WITHIN_BLOB
		cout << "Parallelisation switch: FFEA_PARALLEL_WITHIN_BLOB\n" << endl;
	#endif

	#ifdef FFEA_PARALLEL_PER_BLOB
		 cout << "Parallelisation switch: FFEA_PARALLEL_PER_BLOB\n" << endl;
	#endif
	
	// Get some arguments using boost
	b_po::options_description desc("Allowed options");
	b_po::positional_options_description pdesc;
	b_po::variables_map var_map;

	string script_fname;
	int mode = 0;
	int frames_to_delete = 0;
	int verbose;
#ifdef USE_MPI
  double st,et,st1,et1;

  st1=MPI::Wtime();
#endif
  
	// Options for visible and non-visible cmd line params
	desc.add_options()
		("help,h", "Print usage message")
		("energy,e", "Output energies out at a higher level of detail (to a .ffeaene file)")
		("input-file,i", b_po::value<string>(&script_fname), "Input script filename")
		("delete-frames,l", b_po::value<int>(&frames_to_delete)->default_value(0), "If restarting a simulation, this will delete the final 'arg' frames before restarting")
		("mode,m", b_po::value<int>(&mode)->default_value(0), "Simulation Mode\n\t0 - FFEA\n\t1 - Elastic Network Model\n\t2 - Dynamic Mode Model\n\t3 - Timestep Calculator)\n")
		("verbose,v", b_po::value<int>(&verbose)->default_value(0), "Prints extra details to stdout on what FFEA is doing\n\t0 - Low\n\t1 - Medium\n\t2 - High\n\t3 - Manic")
	;
		
	// 1 input file max! Option invisible (positional arg)
	pdesc.add("input-file", 1);

	// Parse
	try {
		b_po::store(b_po::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), var_map);
		b_po::notify(var_map);

	// User related problems end here
	} catch(exception& e) {
		cerr << e.what() << endl;
	}

	// Boost makes sure options are valid. We check arguments are valid
	
	// Verbosity
	set_verbosity_level(verbose);
	//userInfo::verblevel = verbose;
	cout << userInfo::verblevel << endl;
	// --mode
	if(var_map.count("mode")) {
		if(mode != 0 && mode != 1 && mode != 2 && mode != 3 && mode != 4) {
			FFEA_error_text();		
			cout << "Unrecognised argument to --mode" << endl;
			cout << desc << endl;
			return FFEA_ERROR;
		}
	}

	// --delete_frames
	if(var_map.count("delete_frames")) {
		if(frames_to_delete < 0) {
			FFEA_error_text();		
			cout << "Argument to --delete_frames cannot be less than zero. We can't magic up simulation frames for you! :D" << endl;
			cout << desc << endl;
			return FFEA_ERROR;
		}
	}

	/* Starting actual functionality */

	// Help text is built in to boost	
	if (var_map.count("help")) {  
		cout << desc << endl;
		return 0;
	}

	// Check we have an input script
	if (! var_map.count("input-file")) {  
		cout << "\n\nUsage: ffea [FFEA SCRIPT FILE (.ffea)] [OPTIONS]\n\n\n" << endl;
		cout << desc << endl;
		return FFEA_ERROR;
	}

	// set up a script_fname with the absolute path
	b_fs::path fs_script_fname = script_fname; 
	b_fs::path canonicalPath = b_fs::canonical(fs_script_fname.parent_path());
	fs_script_fname = canonicalPath / fs_script_fname.filename();
	script_fname = fs_script_fname.string(); 
	cout << "Input FFEA script - " << script_fname << "\n";

	// Use it to name the logfile
	b_fs::path fs_log_fname = fs_script_fname;
	fs_log_fname.replace_extension(".log");
	set_log_fname(fs_log_fname.string()); 

	/* Open and begin the logfile
	SimulationParams::checkFileName(userInfo::log_out_fname);
	userInfo::log_out = fopen(userInfo::log_out_fname.c_str(), "w");
	fprintf(userInfo::log_out, "FFEA Log File\n\nScript - %s\n\n", script_fname.c_str()); */

	//The system of all proteins, electrostatics and water
	World *world;

	// Allocate the world
	world = new World();
	//print_normal("lolnormal");
	//print_mid("lolmid");
	//print_high("lolhigh");
	//print_mania("lolmania");
	
	// Initialise the world, loading all blobs, parameters, electrostatics, kinetics etc.
	cout << "Initialising the world:\n" << endl;
	if(world->init(script_fname, frames_to_delete, mode, var_map.count("energy")) == FFEA_ERROR) {
		FFEA_error_text();
		cout << "Errors during initialisation mean World cannot be constructed properly." << endl;

		// Delete the world (oh no!)
		cout << "Deleting world..." << endl;
		delete world;
		cout << "...done. World has been sucessfully destroyed." << endl;

		return FFEA_ERROR;
	}

	/* World is initialised. How shall we run FFEA? */
	int myreturn;
	if(mode == 0) {
		
		// Full FFEA
		cout << "FFEA simulation - Running world...\n" << endl;
		if(world->run() == FFEA_ERROR) {
			FFEA_error_text();
			cout << "Error occurred when running World\n" << endl;
			myreturn = FFEA_ERROR;
		} else {
			cout << "...done\n" << endl;
			myreturn = FFEA_OK;
		}

	} else if(mode == 1 || mode == 2 || mode == 4) {
		
		// Some form of network model
		if(mode == 1) {

			/* Elastic Network Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Elastic Network Model\n***************************************************\n\n";
		} else if (mode == 2) {

			/* Dynamic Mode Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Dynamic Mode Model\n***************************************************\n\n";
		} else if (mode == 4) {

			/* Rotne Prager Dynamic Mode Model */
			cout << "\n\n\n***************************************************\n\tFFEA - Rotne Prager Dynamic Mode Model\n***************************************************\n\n";
		}

		// Read some user input

		// Get a blob list
		cout << "\tFirstly, lets decide which blobs to analyse!" << endl << endl;
		set<int> blobs;
		int ablob;
		int error;
		string buf;
		while(true) {
			error = 0;
			try{
				cout << "\t\tEnter an index for the blob you would like an elastic network model for, or type 'q' to finish?:";

				// Get index as string
				cin >> buf;
				if(buf.compare("q") == 0 or buf.compare("Q") == 0) {
					cout << endl << "\tThat's all the blobs!" << endl;
					break;
				}

				// Convert to int
				ablob = atoi(buf.c_str());

				// Check int
				if(ablob < 0 || ablob >= world->get_num_blobs()) {
					FFEA_error_text();
					cout << "\tPlease enter a valid blob index (0 <= x <" << world->get_num_blobs() << ")" << endl;
					continue;
				} else {
					blobs.insert(ablob);
				}
			} catch (exception &e) {
				FFEA_error_text();
				cout << "Exception occured:" << endl;
				cout << e.what() << endl;
			}
		}

		if(blobs.size() == 0) {
			cout << "No blobs selected. No elastic network models will be found :(" << endl;
			return FFEA_OK;
		}

		// How many modes?
		int num_modes;
		while(true) {
			try{
				cout << "\n\tHow many modes would you like to visualise?:";

				// Get index as string
				cin >> buf;

				// Convert to int
				num_modes = atoi(buf.c_str());

				// Check int
				if(num_modes <= 0) {
					FFEA_error_text();
					cout << "\tYou must choose a whole number i > 0 and i < 3N" << endl;
					continue;
				}
				break;

			} catch(exception &e) {
				FFEA_error_text();
				cout << "Exception occured:" << endl;
				cout << e.what() << endl;
				cout << "\tPlease enter a whole number less than 3 * the total number of nodes in your system and greater than zero." << endl; 
			}
		}
		
		set<int>::iterator it;
		cout << endl << "\tYou've asked for " << num_modes << " modes to be calculated for blob(s):" << endl << "\t";
		for(it = blobs.begin(); it != blobs.end(); ++it) {
			cout << *it << " ";
		}

		if(mode == 1) {

			/* Elastic Network Model */
			cout << endl << endl << "\tBeginning the calculation of the ENMs..." << endl;
			if(world->enm(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating elastic network models." << endl;
				myreturn = FFEA_ERROR;

			} else {
				cout << "...done! " << endl;
				myreturn = FFEA_OK;
			}

		} else if (mode == 2) {

			/* Dynamic Mode Model */
			cout << endl << endl << "\tBeginning the calculation of the DMMs..." << endl;
			if(world->dmm(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating dynamic mode models." << endl;
				myreturn = FFEA_ERROR;
			} else {
				cout << "...done! " << endl;
				myreturn = FFEA_OK;
			}

		} else if (mode == 4) {

			/* Rotne-Prager Dynamic Mode Model */
			cout << endl << endl << "\tBeginning the calculation of the Rotne-Prager DMMs..." << endl;
			if(world->dmm_rp(blobs, num_modes) == FFEA_ERROR) {
				cout << endl;
				FFEA_error_text();
				cout << "Problem when calculating RP dynamic mode models." << endl;
				myreturn = FFEA_ERROR;
			} else {
				cout << "...done! " << endl;
				myreturn = FFEA_OK;
			}
		}

	} else if(mode == 3) {
		
		/* Calculate maximum allowed timestep for system */
		cout << endl << "FFEA - Calculating System Timescales" << endl << endl;
		if(world->get_smallest_time_constants() == FFEA_ERROR) {
			FFEA_error_text();
			cout << "Error occurred in 'get_smallest_time_constants()' in World\n" << endl;
			myreturn = FFEA_ERROR;
		} else {
			cout << "...done\n" << endl;
			myreturn = FFEA_OK;
		}
	} else {
		cout << "Error. Unrecognised argument to --mode \"" << mode << "\"" << endl << endl;
		cout << desc << endl;
	}

	/* Delete the world (oh no!) */
	cout << "Deleting world..." << endl;
#ifdef USE_MPI
  st = MPI::Wtime();
	delete world;
  et = MPI::Wtime()-st;
  et1 = MPI::Wtime()-st1;
  cout<< "benchmarking--------Finalising time of ffea:"<<et<<"seconds"<<endl;
  cout<< "benchmarking--------total executing time:"<<et1<<"seconds"<<endl;
#endif 
	cout << "...done. World has been sucessfully destroyed." << endl;

#ifdef USE_MPI
  MPI::Finalize();
#endif
	return myreturn;
}
