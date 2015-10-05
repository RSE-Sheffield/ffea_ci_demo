#include "World.h"

World::World() {

    // Initialise everything to zero
    blob_array = NULL;
    active_conformation_index = NULL;
    active_state_index = NULL;
    spring_array = NULL;
    kinetic_map_array = NULL;
    kinetic_double_map_array = NULL;
    kinetic_state = NULL;
    kinetic_rate = NULL;
    num_blobs = 0;
    num_conformations = NULL;
    num_springs = 0;
    num_threads = 0;
    rng = NULL;
    phi_Gamma = NULL;
    total_num_surface_faces = 0;
    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;
    trajectory_out = NULL;
    measurement_out = NULL;
}

World::~World() {
    delete[] rng;
    rng = NULL;
    num_threads = 0;

    delete[] blob_array;
    blob_array = NULL;
    num_blobs = 0;
    delete[] num_conformations;
    num_conformations = NULL;

    delete[] active_conformation_index;
    active_conformation_index = NULL;
    delete[] active_state_index;
    active_state_index = NULL;

    delete[] spring_array;
    spring_array = NULL;
    num_springs = 0;
  
    delete[] kinetic_map_array;
    kinetic_map_array = NULL;

    delete[] kinetic_double_map_array;
    kinetic_double_map_array = NULL;

    delete[] kinetic_state;
    kinetic_state = NULL;
    delete[] kinetic_rate;
    kinetic_rate = NULL;

    delete[] phi_Gamma;
    phi_Gamma = NULL;

    total_num_surface_faces = 0;
  
    box_dim.x = 0;
    box_dim.y = 0;
    box_dim.z = 0;
    step_initial = 0;
  
    trajectory_out = NULL;

    delete[] measurement_out;
    measurement_out = NULL;   
}

/**
 * @brief Reads the .ffea file and initialises the World.
 * @param[in] string FFEA_script_filename
 * @details Open and read .ffea file,
 *   parse the <param> block through SimulationParams::extract_params in the 
 *   "SimulationParams params" private attribute,
 *   parse the <blobs> and <springs> blocks through World::read_and_build_system
 * initialise a number of RNG,
 * prepare output files,
 * initialise VdW solver,
 * initialise BEM PBE solver
 * */
int World::init(string FFEA_script_filename, int frames_to_delete) {
	
	// Set some constants and variables
	int i, j;
	const int MAX_BUF_SIZE = 255;
	string buf_string;
	FFEA_input_reader *ffeareader;
	ffeareader = new FFEA_input_reader();

	// Open script
	ifstream fin;
	fin.open(FFEA_script_filename.c_str());

	// Copy entire script into string
	vector<string> script_vector;
	ffeareader->file_to_lines(FFEA_script_filename, &script_vector);

	// Get params section
	cout << "Extracting Parameters..." << endl;
	if(params.extract_params(script_vector) != 0) {
		FFEA_error_text();
		printf("Error parsing parameters in SimulationParams::extract_params()\n");
		return FFEA_ERROR;
	}
	cout << "...done!" << endl;

	// Check for consistency
	cout << "\nVerifying Parameters..." << endl;
	if(params.validate() != 0) {
		FFEA_error_text();
		printf("Parameters found to be inconsistent in SimulationParams::validate()\n");
		return FFEA_ERROR;
	}
	
	// Build kinetic maps if necessary (and rates and binding site matrix)
	if(params.calc_kinetics == 1) {
		kinetic_rate = new scalar**[params.num_blobs];
		kinetic_state = new KineticState*[params.num_blobs];
		kinetic_map_array = new SparseMatrixFixedPattern**[params.num_blobs];
		kinetic_double_map_array = new SparseMatrixFixedPattern***[params.num_blobs];
		for(i = 0; i < params.num_blobs; ++i) {
			kinetic_map_array[i] = new SparseMatrixFixedPattern*[params.num_conformations[i]];
			kinetic_double_map_array[i] = new SparseMatrixFixedPattern**[params.num_conformations[i]];
			for(j = 0; j < params.num_conformations[i]; ++j) {
				kinetic_map_array[i][j] = new SparseMatrixFixedPattern[params.num_conformations[i]];
				kinetic_double_map_array[i][j] = new SparseMatrixFixedPattern*[params.num_conformations[i]];
			}
		}

		kinetic_rng.seed(params.rng_seed);

		if(binding_matrix.init(params.binding_params_fname) == FFEA_ERROR) {
			FFEA_ERROR_MESSG("Error when reading from binding site params file.\n")
		}
		
		
	}

	// Load the vdw forcefield params matrix
    	if (lj_matrix.init(params.vdw_params_fname) == FFEA_ERROR) {
        	FFEA_ERROR_MESSG("Error when reading from vdw forcefeild params file.\n")
    	}

    	// detect how many threads we have for openmp
    	int tid;
	
#ifdef USE_OPENMP
	#pragma omp parallel default(none) private(tid)
	{
        	tid = omp_get_thread_num();
        	if (tid == 0) {
            		num_threads = omp_get_num_threads();
            		printf("\n\tNumber of threads detected: %d\n\n", num_threads);
        	}
    	}
#else
	num_threads = 1;
#endif

    	// We need one rng for each thread, to avoid concurrency problems,
    	// so generate an array of instances of mersenne-twister rngs.
    	rng = new MTRand[num_threads];

    	// Seed each rng differently
    	for (i = 0; i < num_threads; i++) {
        	rng[i].seed(params.rng_seed + i);
	}

	// Build system of blobs, conformations, kinetics etc
	if(read_and_build_system(script_vector) != 0) {
		FFEA_error_text();
		cout << "System could not be built in World::read_and_build_system()" << endl;
		return FFEA_ERROR;
	}

        // If requested, initialise the PreComp_solver. 
        //   Because beads need to be related to elements, it is much easier if 
        //   it is done before moving the blobs to the latest trajectory step in 
        //   case of "restart".
        if (params.calc_preComp ==1)
              pc_solver.init(&pc_params, &params, blob_array);


	// Create measurement files
	measurement_out = new FILE *[params.num_blobs + 1];

	// If not restarting a previous simulation, create new trajectory and measurement files
	if (params.restart == 0) {
		
		// Open the trajectory output file for writing
		if ((trajectory_out = fopen(params.trajectory_out_fname, "w")) == NULL) {
		    FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
		}

		// Open the measurement output file for writing
		for (i = 0; i < params.num_blobs + 1; ++i) {
		    if ((measurement_out[i] = fopen(params.measurement_out_fname[i], "w")) == NULL) {
		        FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
		    }
		}

		// Print initial info stuff
		fprintf(trajectory_out, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\nNumber of Conformations", params.num_blobs);
		for (i = 0; i < params.num_blobs; ++i) {
			fprintf(trajectory_out, " %d", params.num_conformations[i]);
		}
		fprintf(trajectory_out, "\n");

		for (i = 0; i < params.num_blobs; ++i) {
		    fprintf(trajectory_out, "Blob %d:\t", i);
		    for(j = 0; j < params.num_conformations[i]; ++j) {
			    fprintf(trajectory_out, "Conformation %d Nodes %d\t", j, blob_array[i][j].get_num_nodes());
		    }
		    fprintf(trajectory_out, "\n");
		}
		fprintf(trajectory_out, "\n");

		// First line in trajectory data should be an asterisk (used to delimit different steps for easy seek-search in restart code)
		fprintf(trajectory_out, "*\n");

		// First line in measurements file should be a header explaining what quantities are in each column
		for (i = 0; i < params.num_blobs; ++i) {
		    fprintf(measurement_out[i], "# step | KE | PE | CoM x | CoM y | CoM z | L_x | L_y | L_z | rmsd | vdw_area_%d_surface | vdw_force_%d_surface | vdw_energy_%d_surface\n", i, i, i);
		    fflush(measurement_out[i]);
		}
		fprintf(measurement_out[params.num_blobs], "# step ");
		for (i = 0; i < params.num_blobs; ++i) {
		    for (j = i + 1; j < params.num_blobs; ++j) {
		        fprintf(measurement_out[params.num_blobs], "| vdw_area_%d_%d | vdw_force_%d_%d | vdw_energy_%d_%d ", i, j, i, j, i, j);
		    }
		}
		fprintf(measurement_out[params.num_blobs], "\n");
		fflush(measurement_out[params.num_blobs]);

	    } else {
		// Otherwise, seek backwards from the end of the trajectory file looking for '*' character (delimitter for snapshots)
		printf("Restarting from trajectory file %s\n", params.trajectory_out_fname);
		if ((trajectory_out = fopen(params.trajectory_out_fname, "r")) == NULL) {
		    FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
		}

		printf("Reverse searching for 3 asterisks ");
		if(frames_to_delete != 0) {
			printf(", plus an extra %d, ", (frames_to_delete) * 2);
		}
		printf("(denoting %d completely written snapshots)...\n", frames_to_delete + 1);
		if (fseek(trajectory_out, 0, SEEK_END) != 0) {
		    FFEA_ERROR_MESSG("Could not seek to end of file\n")
		}

		// Variable to store position of last asterisk in trajectory file (initialise it at end of file)
		off_t last_asterisk_pos = ftello(trajectory_out);

		int num_asterisks = 0;
		int num_asterisks_to_find = 3 + (frames_to_delete) * 2;
		while (num_asterisks != num_asterisks_to_find) {
		    if (fseek(trajectory_out, -2, SEEK_CUR) != 0) {
		        perror(NULL);
		        FFEA_ERROR_MESSG("It is likely we have reached the begininng of the file whilst trying to delete frames. You can't delete %d frames.\n", frames_to_delete)
		    }
		    char c = fgetc(trajectory_out);
		    if (c == '*') {
		        num_asterisks++;
		        printf("Found %d\n", num_asterisks);

		        // get the position in the file of this last asterisk
		        if (num_asterisks == num_asterisks_to_find - 2) {
		            last_asterisk_pos = ftello(trajectory_out);
		        }
		    }
		}

		char c;
		if ((c = fgetc(trajectory_out)) != '\n') {
		    ungetc(c, trajectory_out);
		}

		printf("Loading Blob position and velocity data from last completely written snapshot \n");
		int blob_id, conformation_id;
		long long rstep;
		for (int b = 0; b < params.num_blobs; b++) {
			if (fscanf(trajectory_out, "Blob %d, Conformation %d, step %lld\n", &blob_id, &conformation_id, &rstep) != 3) {
		        	FFEA_ERROR_MESSG("Error reading header info for Blob %d\n", b)
		    	}
		    if (blob_id != b) {
		        FFEA_ERROR_MESSG("Mismatch in trajectory file - found blob id = %d, was expecting blob id = %d\n", blob_id, b)
		    }
		    printf("Loading node position, velocity and potential from restart trajectory file, for blob %d, step %lld\n", blob_id, rstep);
		    if (active_blob_array[b]->read_nodes_from_file(trajectory_out) == FFEA_ERROR) {
		        FFEA_ERROR_MESSG("Error restarting blob %d\n", blob_id)
		    }
		}
		step_initial = rstep;
		printf("...done. Simulation will commence from step %lld\n", step_initial);
		fclose(trajectory_out);

		// Truncate the trajectory file up to the point of the last asterisk (thereby erasing any half-written time steps that may occur after it)
		printf("Truncating the trajectory file to the last asterisk (to remove any half-written time steps)\n");
		if (truncate(params.trajectory_out_fname, last_asterisk_pos) != 0) {
		    FFEA_ERROR_MESSG("Error when trying to truncate trajectory file %s\n", params.trajectory_out_fname)
		}

		// Open trajectory and measurment files for appending
		printf("Opening trajectory and measurement files for appending.\n");
		if ((trajectory_out = fopen(params.trajectory_out_fname, "a")) == NULL) {
		    FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
		}
		for (i = 0; i < params.num_blobs + 1; ++i) {
		    if ((measurement_out[i] = fopen(params.measurement_out_fname[i], "a")) == NULL) {
		        FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
		    }
		}

		// Append a newline to the end of this truncated trajectory file (to replace the one that may or may not have been there)
		fprintf(trajectory_out, "\n");

		for (i = 0; i < params.num_blobs + 1; ++i) {
		    fprintf(measurement_out[i], "#==RESTART==\n");
		}

	    }


	    // Initialise the Van der Waals solver
	    if(params.calc_vdw == 1 || params.calc_es == 1) {
		vector3 world_centroid, shift;
	        get_system_centroid(&world_centroid);
		if(params.es_N_x < 1 || params.es_N_y < 1 || params.es_N_z < 1) {
			vector3 dimension_vector;
			get_system_dimensions(&dimension_vector);
			
			// Calculate decent box size
			params.es_N_x = 2 * (int)ceil(dimension_vector.x * (params.kappa / params.es_h));
			params.es_N_y = 2 * (int)ceil(dimension_vector.y * (params.kappa / params.es_h));
			params.es_N_z = 2 * (int)ceil(dimension_vector.z * (params.kappa / params.es_h));
			cout << " " << params.es_N_x << " " << params.es_N_y << " " << params.es_N_z << endl;
		}

		// Move to box centre
		box_dim.x = params.es_h * (1.0 / params.kappa) * params.es_N_x;
	        box_dim.y = params.es_h * (1.0 / params.kappa) * params.es_N_y;
	        box_dim.z = params.es_h * (1.0 / params.kappa) * params.es_N_z;
		
		shift.x = box_dim.x / 2.0 - world_centroid.x;
		shift.y = box_dim.y / 2.0 - world_centroid.y;
		shift.z = box_dim.z / 2.0 - world_centroid.z;

		for (i = 0; i < params.num_blobs; i++) {
			active_blob_array[i]->get_centroid(&world_centroid);
			active_blob_array[i]->move(shift.x, shift.y, shift.z);
		}
	    }

	    box_dim.x = params.es_h * (1.0 / params.kappa) * params.es_N_x;
	    box_dim.y = params.es_h * (1.0 / params.kappa) * params.es_N_y;
	    box_dim.z = params.es_h * (1.0 / params.kappa) * params.es_N_z;
	    
	    vdw_solver.init(&lookup, &box_dim, &lj_matrix);

	    // Calculate the total number of vdw interacting faces in the entire system
	    total_num_surface_faces = 0;
	    for (i = 0; i < params.num_blobs; i++) {
		total_num_surface_faces += active_blob_array[i]->get_num_faces();
	    }
	    printf("Total number of surface faces in system: %d\n", total_num_surface_faces);

	    if (params.es_N_x > 0 && params.es_N_y > 0 && params.es_N_z > 0) {

		// Allocate memory for an NxNxN grid, with a 'pool' for the required number of surface faces
		printf("Allocating memory for nearest neighbour lookup grid...\n");
		if (lookup.alloc(params.es_N_x, params.es_N_y, params.es_N_z, total_num_surface_faces) == FFEA_ERROR) {
		    FFEA_error_text();
		    printf("When allocating memory for nearest neighbour lookup grid\n");
		    return FFEA_ERROR;
		}
		printf("...done\n");

		printf("Box has volume %e cubic angstroms\n", (box_dim.x * box_dim.y * box_dim.z) * 1e30);

		// Add all the faces from each Blob to the lookup pool
		printf("Adding all faces to nearest neighbour grid lookup pool\n");
		for (i = 0; i < params.num_blobs; i++) {
		    int num_faces_added = 0;
		    for (j = 0; j < active_blob_array[i]->get_num_faces(); j++) {
		        Face *b_face = active_blob_array[i]->get_face(j);
		        if (b_face != NULL) {
		            if (lookup.add_to_pool(b_face) == FFEA_ERROR) {
		                FFEA_error_text();
		                printf("When attempting to add a face to the lookup pool\n");
		                return FFEA_ERROR;
		            }
		            num_faces_added++;
		        }
		    }
		    printf("%d 'VdW active' faces, from blob %d, added to lookup grid.\n", num_faces_added, i);
		}

		// Initialise the BEM PBE solver
		if (params.calc_es == 1) {
		    printf("Initialising Boundary Element Poisson Boltzmann Solver...\n");
		    PB_solver.init(&lookup);
		    PB_solver.set_kappa(params.kappa);

		    // Initialise the nonsymmetric matrix solver
		    nonsymmetric_solver.init(total_num_surface_faces, params.epsilon2, params.max_iterations_cg);

		    // Allocate memory for surface potential vector
		    phi_Gamma = new scalar[total_num_surface_faces];
		    for (i = 0; i < total_num_surface_faces; i++)
		        phi_Gamma[i] = 0;

		    work_vec = new scalar[total_num_surface_faces];
		    for (i = 0; i < total_num_surface_faces; i++)
		        work_vec[i] = 0;

		    J_Gamma = new scalar[total_num_surface_faces];
		    for (i = 0; i < total_num_surface_faces; i++)
		        J_Gamma[i] = 0;
		}
	    }

	    if (params.restart == 0) {
		// Carry out measurements on the system before carrying out any updates (if this is step 0)
		print_trajectory_and_measurement_files(0, 0);
		print_trajectory_conformation_changes(trajectory_out, 0, NULL, NULL);
	    }

#ifdef FFEA_PARALLEL_WITHIN_BLOB
    printf("Now initialised with 'within-blob parallelisation' (FFEA_PARALLEL_WITHIN_BLOB) on %d threads.\n", num_threads);
#endif

#ifdef FFEA_PARALLEL_PER_BLOB
    printf("Now initialised with 'per-blob parallelisation' (FFEA_PARALLEL_PER_BLOB) on %d threads.\n", num_threads);
#endif

    return FFEA_OK;
}

/* Linearise the system in order to solve exactly for a linear approximation about initial position */
int World::get_smallest_time_constants() {
	
	int blob_index;
	int num_nodes, num_rows;
	double dt_min = INFINITY;
	double dt_max = -1 * INFINITY;
	int dt_max_bin = -1, dt_min_bin = -1;
	cout << "\n\nFFEA mode - Calculate Maximum Allowed Timesteps" << endl << endl;

	// Do active blobs first
	for(blob_index = 0; blob_index < params.num_blobs; ++blob_index) {

		cout << "Blob " << blob_index << ":" << endl << endl;

		// Ignore if we have a static blob
		if(active_blob_array[blob_index]->get_motion_state() == FFEA_BLOB_IS_STATIC) {
			cout << "Blob " << blob_index << " is STATIC. No associated timesteps." << endl;
			continue;
		}

		// Define and reset all required variables for this crazy thing (maybe don't use stack memory??)
		num_nodes = active_blob_array[blob_index]->get_num_linear_nodes();
		num_rows = 3 * num_nodes;

		Eigen::SparseMatrix<double> K(num_rows, num_rows);
		Eigen::SparseMatrix<double> A(num_rows, num_rows);
		Eigen::SparseMatrix<double> K_inv(num_rows, num_rows);
		Eigen::SparseMatrix<double> I(num_rows, num_rows);
		Eigen::SparseMatrix<double> tau_inv(num_rows, num_rows);

		/* Build K */
		cout << "\tCalculating the Global Viscosity Matrix, K...";
		if(active_blob_array[blob_index]->build_linear_node_viscosity_matrix(&K) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_viscosity_matrix' from blob " << blob_index << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		/* Build A */
		cout << "\tCalculating the Global Linearised Elasticity Matrix, A...";
		if(active_blob_array[blob_index]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << blob_index << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		/* Invert K (it's symmetric! Will not work if stokes_visc == 0) */
		cout << "\tAttempting to invert K to form K_inv...";
		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> Cholesky(K); // performs a Cholesky factorization of K
		I.setIdentity();
		K_inv = Cholesky.solve(I);
		if(Cholesky.info() == Eigen::Success) {
			cout << "done! Successful inversion of K!" << endl;
		} else if (Cholesky.info() == Eigen::NumericalIssue) {
			cout << "\nK cannot be inverted via Cholesky factorisation due to numerical issues. You possible don't have an external solvent set." << endl;
			return FFEA_OK;
		} else if (Cholesky.info() == Eigen::NoConvergence) {
			cout << "\nInversion iteration couldn't converge. K must be a crazy matrix. Possibly has zero eigenvalues?" << endl;
			return FFEA_OK;
		}

		/* Apply to A */
		cout << "\tCalculating inverse time constant matrix, tau_inv = K_inv * A...";
		tau_inv = K_inv * A;
		cout << tau_inv << endl;
		tau_inv = A * K_inv;
		cout << endl << tau_inv << endl;
		cout << "done!" << endl;
		exit(0);

		/* Diagonalise */
		cout << "\tDiagonalising tau_inv...";
		Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(tau_inv);
		cout << eigensolver.eigenvalues() << endl;
		exit(0);

		double smallest_val = INFINITY;
		double largest_val = -1 * INFINITY;
		for(int i = 0; i < num_rows; ++i) {

			// Quick fix for weird eigenvalues
			if(eigensolver.eigenvalues()[i].imag() != 0) {
				continue;
			} else {
				if(eigensolver.eigenvalues()[i].real() < smallest_val && eigensolver.eigenvalues()[i].real() > 0) {
					smallest_val = eigensolver.eigenvalues()[i].real();
				} 
				if (eigensolver.eigenvalues()[i].real() > largest_val) {
					largest_val = eigensolver.eigenvalues()[i].real();
				}
			}	
		}
		cout << "done!" << endl;

		// Output
		cout << "\tThe time-constant of the slowest mode in Blob " << blob_index << ", tau_max = " << 1.0 / smallest_val << "s" << endl;
		cout << "\tThe time-constant of the fastest mode in Blob " << blob_index << ", tau_min = " << 1.0 / largest_val << "s" << endl << endl;

		// Global stuff
		if(1.0 / smallest_val > dt_max) {
			dt_max = 1.0 / smallest_val;
			dt_max_bin = blob_index;
		}

		if(1.0 / largest_val < dt_min) {
			dt_min = 1.0 / largest_val;
			dt_min_bin = blob_index;
		}

	}
	cout << "The time-constant of the slowest mode in all blobs, tau_max = " << dt_max << "s, from Blob " << dt_max_bin << endl;
	cout << "The time-constant of the fastest mode in all blobs, tau_min = " << dt_min << "s, from Blob " << dt_min_bin << endl << endl;
	cout << "Remember, the energies in your system will begin to be incorrect long before the dt = tau_min. I'd have dt << tau_min if I were you." << endl;
	return FFEA_OK;
}

/*
 *  Build and output an elastic network model
 */
int World::enm(set<int> blob_indices, int num_modes) {

	set<int> order;
	set<int>::iterator it, it2;
	int i, j, k, num_nodes, num_rows;

	// Calculate the modes for each blob
	for(it = blob_indices.begin(); it != blob_indices.end(); ++it) {

		// Because I'm lazy!
		i = *it;
		
		/* Firstly build the linear elasticity matrix */
		num_nodes = active_blob_array[i]->get_num_linear_nodes();
		num_rows = 3 * num_nodes;
		Eigen::SparseMatrix<double> A(num_rows, num_rows);
		cout << "\t\tCalculating the Global Linearised Elasticity Matrix, A...";
		if(active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		/* Diagonalise it */
		cout << "\t\tDiagonalising A...";
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);
		cout << "done!" << endl;

		/* Find the top 'num_modes' eigenvalues by ordering the list */
		cout << "\t\tOrdering the eigenvalues...";
		int index[num_rows], tempi;
		double evals_ordered[num_rows], tempd;
		for(j = 0; j < num_rows; ++j) {
			index[j] = j;
			evals_ordered[j] = eigensolver.eigenvalues()[j];	
		}

		int fin = 0;
		while(fin == 0) {
			fin = 1;
			for(j = 0; j < num_rows - 1; j++) {
				if(fabs(evals_ordered[j]) < fabs(evals_ordered[j + 1])) {
					fin = 0;
					tempd = evals_ordered[j + 1];
					evals_ordered[j + 1] = evals_ordered[j];
					evals_ordered[j] = tempd;
					tempi = index[j + 1];
					index[j + 1] = index[j];
					index[j] = tempi;
				}	
			}
		}
		cout << "done! " << endl;
		for(j = 0; j < num_rows; ++j) {
			cout << evals_ordered[j] << endl;
		}

		/* Use the eigenvalues to build mode trajectories */

		// Firstly get the overall size of the top mode
		double L = -1 * INFINITY, dx;
		vector3 min, max;
		active_blob_array[i]->get_min_max(&min, &max);
		if(max.x - min.x > L) {
			L = max.x - min.x;
		}
		if(max.y - min.y > L) {
			L = max.y - min.y;
		}
		if(max.z - min.z > L) {
			L = max.z - min.z;
		}

		// L becomes a scaling factor
		L = 0.5 * sqrt(fabs(evals_ordered[0]) / params.kT) * L;
		for(j = 0; j < num_modes; ++j) {

			// See what eigenvector multiplier is
			dx = L * sqrt(params.kT / fabs(evals_ordered[0])) / 25.0;

			// Make 20 frames out of this
			make_trajectory_from_eigenvector(i, j, eigensolver.eigenvectors().col(index[j]), dx);
		}
	}
	return FFEA_OK;
}

/*
 *  Build and output a dynamic mode model
 */
int World::dmm(set<int> blob_indices, int num_modes) {

	set<int> order;
	set<int>::iterator it, it2;
	int i, j, k, num_nodes, num_rows;

	// Calculate the modes for each blob
	for(it = blob_indices.begin(); it != blob_indices.end(); ++it) {

		// Because I'm lazy!
		i = *it;
		
		num_nodes = active_blob_array[i]->get_num_linear_nodes();
		num_rows = 3 * num_nodes;
		Eigen::SparseMatrix<double> K(num_rows, num_rows);
		Eigen::SparseMatrix<double> Q(num_rows, num_rows);
		Eigen::SparseMatrix<double> A(num_rows, num_rows);
		Eigen::MatrixXd Ahat(num_rows, num_rows);
		Eigen::MatrixXd R(num_rows, num_rows);

		/* Build K, the linear viscosity matrix */
		cout << "\tCalculating the Global Viscosity Matrix, K...";
		if(active_blob_array[i]->build_linear_node_viscosity_matrix(&K) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_viscosity_matrix' from blob " << i << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		/* Diagonalise it */
		cout << "\t\tDiagonalising K...";
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esK(K);
		cout << "done!" << endl;
		cout << esK.eigenvalues() << endl;
		//exit(0);

		/* Invert the and sqrt the elements of K evals to form Q */
		std::vector<Eigen::Triplet<double>> vals;
		for(j = 0; j < num_rows; ++j) {

			// This will catch the zero eigenvalues i.e. ridiculously small ones
			if(esK.eigenvalues()[j] < 0) {
				vals.push_back(Eigen::Triplet<double>(j,j, 1.0 / sqrt(-1 * esK.eigenvalues()[j])));
			} else {
				vals.push_back(Eigen::Triplet<double>(j,j, 1.0 / sqrt(esK.eigenvalues()[j])));
			}		
		}
		Q.setFromTriplets(vals.begin(), vals.end());

		/* Build the linear elasticity matrix and get symmetric part (just in case) */
		cout << "\t\tCalculating the Global Linearised Elasticity Matrix, A...";
		if(active_blob_array[i]->build_linear_node_elasticity_matrix(&A) == FFEA_ERROR) {
			cout << endl;
			FFEA_error_text();
			cout << "In function 'Blob::build_linear_node_elasticity_matrix'" << i << endl;
			return FFEA_ERROR;
		}
		cout << "done!" << endl;

		/* Build the transformation matrix Ahat */
		cout << "\t\tBuilding the transformation matrix, Ahat..." << endl;
		Ahat = Q.transpose() * esK.eigenvectors().transpose() * A * esK.eigenvectors() * Q;
		cout << "done!" << endl;

		/* Diagonalise it */
		cout << "\t\tDiagonalising Ahat...";
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esAhat(Ahat);
		cout << "done!" << endl;

		/* Finally, build the dynamical mode matrix R */
		R = esK.eigenvectors() * Q * esAhat.eigenvectors();
 
		/* Find the top 'num_modes' eigenvalues by ordering the list */
		cout << "\t\tOrdering the eigenvalues...";
		int index[num_rows], tempi;
		double evals_ordered[num_rows], tempd;
		for(j = 0; j < num_rows; ++j) {
			index[j] = j;
			evals_ordered[j] = esAhat.eigenvalues()[j];	
		}

		int fin = 0;
		while(fin == 0) {
			fin = 1;
			for(j = 0; j < num_rows - 1; j++) {
				if(evals_ordered[j] < evals_ordered[j + 1]) {
					fin = 0;
					tempd = evals_ordered[j + 1];
					evals_ordered[j + 1] = evals_ordered[j];
					evals_ordered[j] = tempd;
					tempi = index[j + 1];
					index[j + 1] = index[j];
					index[j] = tempi;
				}	
			}
		}
		cout << "done! " << endl;
		for(int j = 0; j < num_rows; ++j) {
			cout << evals_ordered[j] << endl;
		}
		/* Use the eigenvalues to build mode trajectories */

		// Firstly get the overall size of the top mode
		double L = -1 * INFINITY, dx;
		vector3 min, max;
		active_blob_array[i]->get_min_max(&min, &max);
		if(max.x - min.x > L) {
			L = max.x - min.x;
		}
		if(max.y - min.y > L) {
			L = max.y - min.y;
		}
		if(max.z - min.z > L) {
			L = max.z - min.z;
		}

		// L becomes a scaling factor
		L = 0.5 * sqrt(fabs(evals_ordered[0]) / params.kT) * L;
		for(j = 0; j < num_modes; ++j) {

			// See what eigenvector multiplier is
			dx = L * sqrt(params.kT / fabs(evals_ordered[0])) / 25.0;

			// Make 20 frames out of this
			make_trajectory_from_eigenvector(i, j, esAhat.eigenvectors().col(index[j]), dx);
		}
	}
	return FFEA_OK;
}

/*
 * Update entire World for num_steps time steps
 * */
int World::run() {
    int es_count = 1;
    double wtime = omp_get_wtime();
    for (long long step = step_initial; step < params.num_steps; step++) {

        // Zero the force across all blobs
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(runtime)
#endif
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->zero_force();
            if (params.calc_vdw == 1) {
                active_blob_array[i]->zero_vdw_bb_measurement_data();
            }
            if (params.sticky_wall_xz == 1) {
                active_blob_array[i]->zero_vdw_xz_measurement_data();
            }
        }

#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(guided) shared(stderr)
#endif
        for (int i = 0; i < params.num_blobs; i++) {
	
            // If blob centre of mass moves outside simulation box, apply PBC to it
            vector3 com;
            active_blob_array[i]->get_centroid(&com);
            scalar dx = 0, dy = 0, dz = 0;
            int check_move = 0;

            if (com.x < 0) {
                if (params.wall_x_1 == WALL_TYPE_PBC) {
                    dx += box_dim.x;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.x > box_dim.x) {
                if (params.wall_x_2 == WALL_TYPE_PBC) {
                    dx -= box_dim.x;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (com.y < 0) {
                if (params.wall_y_1 == WALL_TYPE_PBC) {
                    dy += box_dim.y;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.y > box_dim.y) {
                if (params.wall_y_2 == WALL_TYPE_PBC) {
                    dy -= box_dim.y;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (com.z < 0) {
                if (params.wall_z_1 == WALL_TYPE_PBC) {
                    dz += box_dim.z;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            } else if (com.z > box_dim.z) {
                if (params.wall_z_2 == WALL_TYPE_PBC) {
                    dz -= box_dim.z;
                    //					printf("fuck\n");
                    check_move = 1;
                }
            }
            if (check_move == 1) {
                active_blob_array[i]->move(dx, dy, dz);
            }

            // If Blob is near a hard wall, prevent it from moving further into it
            active_blob_array[i]->enforce_box_boundaries(&box_dim);
        }



        if (params.calc_es == 1 || params.calc_vdw == 1 || params.sticky_wall_xz == 1) {
            if (es_count == params.es_update) {


#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) schedule(guided)
#endif
                for (int i = 0; i < params.num_blobs; i++) {
                    active_blob_array[i]->calc_centroids_and_normals_of_all_faces();
                    active_blob_array[i]->reset_all_faces();
                }


                // Attempt to place all faces in the nearest neighbour lookup table
                if (lookup.build_nearest_neighbour_lookup(params.es_h * (1.0 / params.kappa)) == FFEA_ERROR) {
                    FFEA_error_text();
                    printf("When trying to place faces in nearest neighbour lookup table.\n");

                    // attempt to print out the final (bad) time step
                    printf("Dumping final step:\n");
                    print_trajectory_and_measurement_files(step, wtime);
		    print_trajectory_conformation_changes(trajectory_out, 0, NULL, NULL);

                    return FFEA_ERROR;
                }

                /*if (params.calc_vdw == 1) {
                    vdw_solver.solve();
                }*/
                if (params.sticky_wall_xz == 1) {
                    vdw_solver.solve_sticky_wall(params.es_h * (1.0 / params.kappa));
                }

                if (params.calc_es == 1) {
                    do_es();
                }

                es_count = 1;
            } else
                es_count++;
        }
        if (params.calc_vdw == 1) vdw_solver.solve();

        // Update all Blobs in the World

        // Set node forces to zero
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->set_forces_to_zero();
        }

        // Apply springs directly to nodes
        apply_springs();

        // if PreComp is required: 
        if (params.calc_preComp == 1) {
          pc_solver.solve();
        }
	
        // Sort internal forces out
        int fatal_errors = 0;
	
#ifdef FFEA_PARALLEL_PER_BLOB
#pragma omp parallel for default(none) shared(step, wtime) reduction(+: fatal_errors) schedule(runtime)
#endif

        for (int i = 0; i < params.num_blobs; i++) {
        	if (active_blob_array[i]->update() == FFEA_ERROR) {
        		FFEA_error_text();
               		printf("A problem occurred when updating Blob %d on step %lld\n", i, step);
                	printf("Simulation ran for %2f seconds (wall clock time) before error ocurred\n", (omp_get_wtime() - wtime));
                	//return FFEA_ERROR;
                	fatal_errors++;
            	}
        }
        if (fatal_errors > 0) {
            FFEA_error_text();
            printf("Detected %d fatal errors in this system update. Exiting now...\n", fatal_errors);

            // attempt to print out the final (bad) time step
            printf("Dumping final step:\n");
            print_trajectory_and_measurement_files(step, wtime);
	    print_trajectory_conformation_changes(trajectory_out, 0, NULL, NULL);

            return FFEA_ERROR;
        }

        if (step % params.check == 0) {
            print_trajectory_and_measurement_files(step + 1, wtime);
        }

	/* Kinetic Part of each step */

	// Get a list of old conformations for writing to traj
	int *from = new int[params.num_blobs];
	for(int i = 0; i < params.num_blobs; ++i) {
		from[i] = active_blob_array[i]->conformation_index;
	}

	if (params.calc_kinetics == 1 && step % params.kinetics_update == 0) {

		// Get base rates
		scalar ***switching_probs;
		switching_probs = new scalar**[params.num_blobs];
		for(int i = 0; i < params.num_blobs; ++i) {
			switching_probs[i] = new scalar*[params.num_states[i]];
			for(int j = 0; j < params.num_states[i]; ++j) {
				switching_probs[i][j] = new scalar[params.num_states[i]];
				for(int k = 0; k < params.num_states[i]; ++k) {
					if(params.num_states[i] > 1) {
						switching_probs[i][j][k] = kinetic_rate[i][j][k];
					} else {
						switching_probs[i][j][k] = 1;
					}	
				}
			}
		}

		// Rescale rates based on energies
		if(rescale_kinetic_rates(switching_probs) != FFEA_OK) {

			delete[] from;
			delete[] switching_probs;
			cout << "Error in 'rescale_kinetic_rates' function." << endl;
			return FFEA_ERROR;
		}
		
		// Change conformation based on these rates
		// Make some bins for the random number generator
		scalar prob_sum;
		scalar **bin_limits = new scalar*[params.num_blobs];
		for(int i = 0; i < params.num_blobs; ++i) {
			bin_limits[i] = new scalar[params.num_states[i] + 1];
			prob_sum = 0.0;
			for(int j = 0; j < params.num_states[i]; ++j) {
				bin_limits[i][j] = prob_sum;
				prob_sum += switching_probs[i][active_state_index[i]][j];
			}
			bin_limits[i][params.num_states[i]] = prob_sum;	
		}

		// Changing states now, if necessary
		scalar switch_check;
		for(int i = 0; i < params.num_blobs; ++i) {

			// Get a random number between 0 and 1
			switch_check = kinetic_rng.rand();
			for(int j = 0; j < params.num_states[i]; ++j) {
				
				// Change state based on this random number!
				if(switch_check >= bin_limits[i][j] && switch_check < bin_limits[i][j + 1]) {
					change_blob_state(i, j);
				}
			}
		}

		delete[] bin_limits;
		delete[] switching_probs;
	}

	// Get a list of new conformations for writing to traj
	int *to = new int[params.num_blobs];
	for(int i = 0; i < params.num_blobs; ++i) {
		to[i] = active_blob_array[i]->conformation_index;
	}

	// Print Trajectory Conf change data and free data
	if(step % params.check == 0) {
		print_trajectory_conformation_changes(trajectory_out, step + 1, from, to);
	}
	
	delete[] from;
	delete[] to;
    }
    printf("Time taken: %2f seconds\n", (omp_get_wtime() - wtime));

    return FFEA_OK;
}

int World::change_blob_state(int blob_index, int new_state_index) {

	// First check if work needs to be done at all
	if(new_state_index == active_state_index[blob_index]) {
		return FFEA_OK;
	}

	// Now check on type of state change
	int switch_type;
	if(kinetic_state[blob_index][active_state_index[blob_index]].conformation_index != kinetic_state[blob_index][new_state_index].conformation_index) {
		
		// Flag for conformation change
		switch_type = FFEA_CONFORMATION_CHANGE;

	} else if (kinetic_state[blob_index][active_state_index[blob_index]].bound == 0 && kinetic_state[blob_index][new_state_index].bound == 1) {
		
		// Flag for a binding event
		switch_type = FFEA_BINDING_EVENT;
	
	} else if (kinetic_state[blob_index][active_state_index[blob_index]].bound == 1 && kinetic_state[blob_index][new_state_index].bound == 0) {
		
		// Flag for an ubinding event
		switch_type = FFEA_UNBINDING_EVENT;
	
	} else {
		
		// Flag for an identity event
		switch_type = FFEA_IDENTITY_EVENT;
	}

	// Let's change states

	// Stuff for all state changes
	active_state_index[blob_index] = new_state_index;

	// Then type specific stuff
	switch(switch_type) {
		
		case FFEA_CONFORMATION_CHANGE:
		{	
			int from_conf = active_conformation_index[blob_index];
			int to_conf = kinetic_state[blob_index][new_state_index].conformation_index;
			vector3 **from_nodes = active_blob_array[blob_index]->get_actual_node_positions();
			vector3 **to_nodes = blob_array[blob_index][to_conf].get_actual_node_positions();

			// Apply map
			kinetic_map_array[blob_index][from_conf][to_conf].block_apply(from_nodes, to_nodes);
			
			// Reassign active_blob and conformation
			active_conformation_index[blob_index] = to_conf;
			active_blob_array[blob_index] = &blob_array[blob_index][to_conf];

			// Activate springs for new conformation
			activate_springs();
			break;
		}
		case FFEA_BINDING_EVENT:
		{
			// Calculate properties for binding sites on the current blob
			int i, j, k;
			scalar average_size;
			vector3 separation;

			// Skip if blob is static
			if(active_blob_array[blob_index]->get_motion_state() == FFEA_BLOB_IS_STATIC) {
				break;
			}


			// Scan over all binding sites
			for(i = 0; i < active_blob_array[blob_index]->num_binding_sites; ++i) {

				// If wrong site type, ignore
				if(active_blob_array[blob_index]->binding_site[i].site_type != kinetic_state[blob_index][new_state_index].binding_site_type_from) {
					continue;
				}

				// Calculate properties for binding sites on the current blob
				active_blob_array[blob_index]->binding_site[i].calc_site_shape();

				// Scan over all binding sites on all other blobs
				for(j = 0; j < params.num_blobs; ++j) {
					
					// Skip binding to self (for now)
					if(j == blob_index) {
						continue;
					}
					
					for(k = 0; k < active_blob_array[j]->num_binding_sites; ++k) {
			
						// If wrong site type, ignore
						if(active_blob_array[j]->binding_site[k].site_type != kinetic_state[blob_index][new_state_index].binding_site_type_to) {
							continue;
						}
		
						// Check if binding sites are compatible
						if(binding_matrix.interaction[active_blob_array[blob_index]->binding_site[i].site_type][active_blob_array[j]->binding_site[k].site_type] != 1) {
							continue;
						}

						// Calculate properties for potential binding site and their coupling
						active_blob_array[j]->binding_site[k].calc_site_shape();
						average_size = (active_blob_array[blob_index]->binding_site[i].length + active_blob_array[j]->binding_site[k].length) / 2.0;
						separation.x = active_blob_array[blob_index]->binding_site[i].centroid.x - active_blob_array[j]->binding_site[k].centroid.x;
						separation.y = active_blob_array[blob_index]->binding_site[i].centroid.y - active_blob_array[j]->binding_site[k].centroid.y;
						separation.z = active_blob_array[blob_index]->binding_site[i].centroid.z - active_blob_array[j]->binding_site[k].centroid.z;

						// Bind if close enough!
						if(mag(&separation) < average_size) {
							active_blob_array[blob_index]->kinetic_bind(i);		
						}
					}
				}
				
				
			}
			
			break;
		}
		case FFEA_UNBINDING_EVENT:
			// Calculate properties for binding sites on the current blob
			int i, j, k;
			scalar average_size;
			vector3 separation;

			// Skip if blob is static
			if(active_blob_array[blob_index]->get_motion_state() == FFEA_BLOB_IS_STATIC) {
				break;
			}


			// Scan over all binding sites
			for(i = 0; i < active_blob_array[blob_index]->num_binding_sites; ++i) {

				// If wrong site type, ignore
				if(active_blob_array[blob_index]->binding_site[i].site_type != kinetic_state[blob_index][new_state_index].binding_site_type_from) {
					continue;
				}

				// Calculate properties for binding sites on the current blob
				active_blob_array[blob_index]->binding_site[i].calc_site_shape();

				// Scan over all binding sites on all other blobs
				for(j = 0; j < params.num_blobs; ++j) {
					
					// Skip binding to self (for now)
					if(j == blob_index) {
						continue;
					}
					
					for(k = 0; k < active_blob_array[j]->num_binding_sites; ++k) {
			
						// If wrong site type, ignore
						if(active_blob_array[j]->binding_site[k].site_type != kinetic_state[blob_index][new_state_index].binding_site_type_to) {
							continue;
						}
		
						// Check if binding sites are compatible
						if(binding_matrix.interaction[active_blob_array[blob_index]->binding_site[i].site_type][active_blob_array[j]->binding_site[k].site_type] != 1) {
							continue;
						}

						// Calculate properties for potential binding site and their coupling
						active_blob_array[j]->binding_site[k].calc_site_shape();
						average_size = (active_blob_array[blob_index]->binding_site[i].length + active_blob_array[j]->binding_site[k].length) / 2.0;
						separation.x = active_blob_array[blob_index]->binding_site[i].centroid.x - active_blob_array[j]->binding_site[k].centroid.x;
						separation.y = active_blob_array[blob_index]->binding_site[i].centroid.y - active_blob_array[j]->binding_site[k].centroid.y;
						separation.z = active_blob_array[blob_index]->binding_site[i].centroid.z - active_blob_array[j]->binding_site[k].centroid.z;

						// Bind if close enough!
						if(mag(&separation) < average_size) {
							active_blob_array[blob_index]->kinetic_unbind(i);		
						}
					}
				}
				
				
			}
			
			break;
	
		case FFEA_IDENTITY_EVENT:

			break;

			
	}

	return FFEA_OK;
}

/** 
 * @brief Parses <blobs>, <springs> and <precomp>. 
 * @param[in] vector<string> script_vector, which is essentially the FFEA input file,
 *            line by line, as it comes out of FFEA_input_reader::file_to_lines
 */
int World::read_and_build_system(vector<string> script_vector) {

	// Create some blobs based on params
	cout << "\tCreating blob array..." << endl;
   	blob_array = new Blob*[params.num_blobs];
	active_blob_array = new Blob*[params.num_blobs];
	active_conformation_index = new int[params.num_blobs];
	active_state_index = new int[params.num_blobs];
	for (int i = 0; i < params.num_blobs; ++i) {
	        blob_array[i] = new Blob[params.num_conformations[i]];
	        active_blob_array[i] = &blob_array[i][0];
		active_conformation_index[i] = 0;
		active_state_index[i] = 0;
	}


	// Reading variables
	FFEA_input_reader *systemreader = new FFEA_input_reader();
	int i, j;
	string tag, lrvalue[2], maplvalue[2];
	vector<string> blob_vector, interactions_vector, conformation_vector, kinetics_vector, map_vector, param_vector, spring_vector, binding_vector;
	vector<string>::iterator it;
	
	vector<string> nodes, topology, surface, material, stokes, vdw, binding, pin, maps, beads;
	string states, rates, map_fname;
	int map_indices[2];
	int set_motion_state = 0, set_nodes = 0, set_top = 0, set_surf = 0, set_mat = 0, set_stokes = 0, set_vdw = 0, set_binding = 0, set_pin = 0, set_solver = 0, set_preComp = 0, set_scale = 0;
	scalar scale = 1;
	int solver = FFEA_NOMASS_CG_SOLVER;
	vector<int> motion_state, maps_conf_index_to, maps_conf_index_from;
	vector<int>::iterator maps_conf_ind_it;

	scalar *centroid = NULL, *velocity = NULL, *rotation = NULL;

	// Get Interactions Lines
	int error_code;
	error_code = systemreader->extract_block("interactions", 0, script_vector, &interactions_vector);
	if(error_code == FFEA_ERROR) {
		return FFEA_ERROR;
	} else if (error_code == FFEA_CAUTION) {

		// Block doesn't exist. Maybe add a protection to ensure binding later
	} else {

		// Get spring data
		systemreader->extract_block("springs", 0, interactions_vector, &spring_vector);

		// Error check
		if (spring_vector.size() > 1) {
			FFEA_error_text();
			cout << "'Spring' block should only have 1 file." << endl;
			return FFEA_ERROR; 
		} else if (spring_vector.size() == 1) {
			systemreader->parse_tag(spring_vector.at(0), lrvalue);

			if(load_springs(lrvalue[1].c_str()) != 0) {
				FFEA_error_text();
				cout << "Problem loading springs from " << lrvalue[1] << "." << endl;
				return FFEA_ERROR; 
			}
		}

	       // Get precomputed data
	       pc_params.dist_to_m = 1;
	       pc_params.E_to_J = 1;
               vector<string> precomp_vector;
               systemreader->extract_block("precomp", 0, script_vector, &precomp_vector);
	
               for (i=0; i<precomp_vector.size(); i++){
                 systemreader->parse_tag(precomp_vector[i], lrvalue);
		 if (lrvalue[0] == "types") {
                   lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
                   boost::trim(lrvalue[1]);
                   systemreader->split_string(lrvalue[1], pc_params.types, ",");
                 } else if (lrvalue[0] == "inputData") {
                   pc_params.inputData = stoi(lrvalue[1]);
                 } else if (lrvalue[0] == "approach") {
                   pc_params.approach = lrvalue[1];
                 } else if (lrvalue[0] == "folder") {
                   pc_params.folder = lrvalue[1];
                 } else if (lrvalue[0] == "dist_to_m") {
                   pc_params.dist_to_m = stod(lrvalue[1]);
                 } else if (lrvalue[0] == "E_to_J") {
                   pc_params.E_to_J = stod(lrvalue[1]);
                 }
               }
        }



	// Read in each blob one at a time
	for(i = 0; i < params.num_blobs; ++i) {

		// Get blob data
		systemreader->extract_block("blob", i, script_vector, &blob_vector);

		// Read all conformations
		for(j = 0; j < params.num_conformations[i]; ++j) {

			// Get conformation data
			systemreader->extract_block("conformation", j, blob_vector, &conformation_vector);

			// Error check
			if(conformation_vector.size() == 0) {
				FFEA_error_text();
				cout << " In 'Blob' block " << i << ", expected at least a single 'conformation' block." << endl;
				return FFEA_ERROR;  
			}

			// Parse conformation data
			for(it = conformation_vector.begin(); it != conformation_vector.end(); ++it) {
				systemreader->parse_tag(*it, lrvalue);

				// Assign if possible
				if(lrvalue[0] == "motion_state") {

					if(lrvalue[1] == "DYNAMIC") {
						motion_state.push_back(FFEA_BLOB_IS_DYNAMIC);
					} else if(lrvalue[1] == "STATIC") {
						motion_state.push_back(FFEA_BLOB_IS_STATIC);
					} else if(lrvalue[1] == "FROZEN") {
						motion_state.push_back(FFEA_BLOB_IS_FROZEN);
					}
					set_motion_state = 1;
				} else if (lrvalue[0] == "nodes") {
					nodes.push_back(lrvalue[1]);
					set_nodes = 1;
				} else if (lrvalue[0] == "topology") {
					topology.push_back(lrvalue[1]);
					set_top = 1;
				} else if (lrvalue[0] == "surface") {
					surface.push_back(lrvalue[1]);
					set_surf = 1;
				} else if (lrvalue[0] == "material") {
					material.push_back(lrvalue[1]);
					set_mat = 1;
				} else if (lrvalue[0] == "stokes") {
					stokes.push_back(lrvalue[1]);
					set_stokes = 1;
				} else if (lrvalue[0] == "vdw") {
					vdw.push_back(lrvalue[1]);
					set_vdw = 1;
				} else if (lrvalue[0] == "binding_sites") {
					binding.push_back(lrvalue[1]);
					set_binding = 1;
				} else if (lrvalue[0] == "pin") {
					pin.push_back(lrvalue[1]);
					set_pin = 1;
                                } else if (lrvalue[0] == "beads") {
					beads.push_back(lrvalue[1]);
					set_preComp = 1;
				} else {
					FFEA_error_text();
					cout << "Unrecognised conformation lvalue" << endl;
					return FFEA_ERROR;
				}
			}
		
			// Error check
			if (set_motion_state == 0) {
				FFEA_error_text();
				cout << "Blob " << i << ", conformation " << j << " must have a motion state set.\nAccepted states: DYNAMIC, STATIC, FROZEN." << endl;
				return FFEA_ERROR;
			} else {
				if(set_nodes == 0 || set_surf == 0 || set_vdw == 0) {
					FFEA_error_text();
					cout << "In blob " << i << ", conformation " << j << ":\nFor any blob conformation, 'nodes', 'surface' and 'vdw' must be set." << endl;
					return FFEA_ERROR; 

				} 
				if(set_preComp == 0) {
                                        beads.push_back("");
                                } 

				if(motion_state.back() == FFEA_BLOB_IS_DYNAMIC) {
					if(set_top == 0 || set_mat == 0 || set_stokes == 0 || set_pin == 0) {
						FFEA_error_text();
						cout << "In blob " << i << ", conformation " << j << ":\nFor a DYNAMIC blob conformation, 'topology', 'material', 'stokes' and 'pin' must be set." << endl;
						return FFEA_ERROR; 
					}
				} else {
					topology.push_back("");
					material.push_back("");
					stokes.push_back("");
					pin.push_back("");
					set_top = 1;
					set_mat = 1;
					set_stokes = 1;
					set_pin = 1;
				}

				if(set_preComp == 0) {
					beads.push_back("");
				}
			}

			if(params.calc_kinetics == 0) {
				if(set_binding == 1) {
					cout << "In blob " << i << ", conformation " << j << ":\nIf 'calc_kinetics' is set to 0, 'binding_sites' will be ignored." << endl;
					binding.pop_back();
				}
				binding.push_back("");
			} else {
				if(set_binding == 0) {
					FFEA_error_text();
					cout << "In blob " << i << ", conformation " << j << ":\nIf 'calc_kinetics' is set to 1, 'binding_sites' must be set (can contain 0 sites if necessary)." << endl;
					return FFEA_ERROR;
				}
			}
			
			// Clear conformation vector and set values for next round
			set_nodes = 0;
			set_top = 0;
			set_surf = 0;
			set_mat = 0;
			set_stokes = 0;
			set_vdw = 0;
			set_binding = 0;
			set_pin = 0;
			set_preComp = 0;
			conformation_vector.clear();
		}

		// Read kinetic info if necessary
		if(params.calc_kinetics == 1 && params.num_states[i] > 1) {

			// Get kinetic data
			systemreader->extract_block("kinetics", 0, blob_vector, &kinetics_vector);
			
			// Get map info if necessary
			if(params.num_conformations[i] > 1) {

				// Get map data
				systemreader->extract_block("maps", 0, kinetics_vector, &map_vector);

				// Parse map data
				for(it = map_vector.begin(); it != map_vector.end(); ++it) {
					systemreader->parse_map_tag(*it, map_indices, &map_fname);
					maps.push_back(map_fname);
					maps_conf_index_from.push_back(map_indices[0]);
					maps_conf_index_to.push_back(map_indices[1]);
				}

				// Error check
				if(maps.size() != params.num_conformations[i] * (params.num_conformations[i] - 1)) {
					FFEA_error_text();
					cout << "In blob " << i << ", expected " << params.num_conformations[i] * (params.num_conformations[i] - 1) << " maps to describe all possible switches.\n Read " << maps.size() << " maps." << endl;
					return FFEA_ERROR;
				}
			}

			// Then, states and rates data
			for(it = kinetics_vector.begin(); it != kinetics_vector.end(); ++it) {
				systemreader->parse_tag(*it, lrvalue);
				if(lrvalue[0] == "maps" || lrvalue[0] == "/maps") {
					continue;
				} else if(lrvalue[0] == "states") {
					states = lrvalue[1];
				} else if (lrvalue[0] == "rates") {
					rates = lrvalue[1];
				}
			}

			// Error check
			if(states.size() == 0 || rates.size() == 0) {
				FFEA_error_text();
				cout << "In blob " << i << ", for kinetics, we require both a 'states' file and a 'rates' file." << endl;
				return FFEA_ERROR;
			}
		}		

		// Finally, get the extra blob data (solver, scale, centroid etc)

		int rotation_type = -1;

		for(it = blob_vector.begin(); it != blob_vector.end(); ++it) {
			systemreader->parse_tag(*it, lrvalue);

			if(lrvalue[0] == "conformation" || lrvalue[0] == "/conformation" || lrvalue[0] == "kinetics" || lrvalue[0] == "/kinetics" || lrvalue[0] == "maps" || lrvalue[0] == "/maps") {
				continue;
			} else if(lrvalue[0] == "solver") {
				if(lrvalue[1] == "CG") {
					solver = FFEA_ITERATIVE_SOLVER;
				} else if (lrvalue[1] == "CG_nomass") {
					solver = FFEA_NOMASS_CG_SOLVER;
				} else if (lrvalue[1] == "direct") {
					solver = FFEA_DIRECT_SOLVER;
				} else if (lrvalue[1] == "masslumped") {
					solver = FFEA_MASSLUMPED_SOLVER;
				} else {
					FFEA_error_text();
					cout << "In blob " << i << ", unrecognised solver type.\nRecognised solvers:CG, CG_nomass, direct, masslumped." << endl;
					return FFEA_ERROR;
				}
				set_solver = 1;

			} else if(lrvalue[0] == "scale") {
				scale = atof(lrvalue[1].c_str());
				set_scale = 1;
			} else if(lrvalue[0] == "centroid" || lrvalue[0] == "centroid_pos") {
				
				centroid = new scalar[3];

				lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
				boost::trim(lrvalue[1]);
				systemreader->split_string(lrvalue[1], centroid, ",");

			} else if(lrvalue[0] == "velocity") {
				
				velocity = new scalar[3];
				
				lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
				boost::trim(lrvalue[1]);
				systemreader->split_string(lrvalue[1], velocity, ",");

			} else if(lrvalue[0] == "rotation") {

				rotation = new scalar[9];
				lrvalue[1] = boost::erase_last_copy(boost::erase_first_copy(lrvalue[1], "("), ")");
				boost::trim(lrvalue[1]);
				if(systemreader->split_string(lrvalue[1], rotation, ",") == 3) {
					rotation_type = 0;
				} else {
					rotation_type = 1;
				}	
			}
		}

		// Error checking
		if(set_solver == 0) {
			cout << "\tBlob " << i << ", solver not set. Defaulting to CG_nomass." << endl;
			solver = FFEA_NOMASS_CG_SOLVER;
			set_solver = 1;
		}

		if(set_scale == 0) {
			FFEA_error_text();
			cout << "Blob " << i << ", scale not set." << endl;
			return FFEA_ERROR;
		}

		// Build blob
		//Build conformations
		for(j = 0; j < params.num_conformations[i]; ++j) {
			cout << "\tInitialising blob " << i << " conformation " << j << "..." << endl;
			if (blob_array[i][j].init(i, j, nodes.at(j).c_str(), topology.at(j).c_str(), surface.at(j).c_str(), material.at(j).c_str(), stokes.at(j).c_str(), vdw.at(j).c_str(), binding.at(j).c_str(), pin.at(j).c_str(), beads.at(j).c_str(), 
                       		scale, solver, motion_state.at(j), &params, &pc_params, &lj_matrix, &binding_matrix, rng, num_threads) == FFEA_ERROR) {
                       		FFEA_error_text();
                        	cout << "\tError when trying to initialise Blob " << i << ", conformation " << j << "." << endl;
                    		return FFEA_ERROR;
               		}

			// if centroid position is set, position the blob's centroid at that position. If vdw is set, move to center of box
                	if (centroid != NULL) {

                    		// Rescale first	
                    		centroid[0] *= scale;
                    		centroid[1] *= scale;
                    		centroid[2] *= scale;
                    		vector3 dv = blob_array[i][j].position(centroid[0], centroid[1], centroid[2]);
                                // if Blob has a number of beads, transform them too:
                                if (blob_array[i][j].get_num_beads() > 0)
                                  blob_array[i][j].position_beads(dv.x, dv.y, dv.z);
                	}
                      
                	if(rotation != NULL) {
				if(rotation_type == 0) {
                                   if (blob_array[i][j].get_num_beads() > 0) {
					blob_array[i][j].rotate(rotation[0], rotation[1], rotation[2]);
                                    } else {
					blob_array[i][j].rotate(rotation[0], rotation[1], rotation[2]);
                                    }
				} else {
                                   if (blob_array[i][j].get_num_beads() > 0) {
	                    		blob_array[i][j].rotate(rotation[0], rotation[1], rotation[2], rotation[3], rotation[4], rotation[5], rotation[6], rotation[7], rotation[8]);
                                   } else { 
                                      // if Blob has a number of beads, transform them too:
                                        blob_array[i][j].rotate(rotation[0], rotation[1], rotation[2], rotation[3], rotation[4], rotation[5], rotation[6], rotation[7], rotation[8], 1);
                                   }
				}        
	        	}

                	if (velocity != NULL)
                    		blob_array[i][j].velocity_all(velocity[0], velocity[1], velocity[2]);

                	// set the current node positions as pos_0 for this blob, so that all rmsd values
                	// are calculated relative to this conformation centred at this point in space.
                	blob_array[i][j].set_rmsd_pos_0();

			cout << "\t...done!" << endl;
		}

		// Build kinetics
		if(params.calc_kinetics == 1 && params.num_states[i] > 1) {
			cout << "\tInitialising kinetics for blob " << i << "..." << endl;
			if(load_kinetic_maps(maps, maps_conf_index_from, maps_conf_index_to, i) == FFEA_ERROR) {
				FFEA_error_text();
				cout << "Problem reading kinetic maps in 'read_kinetic_maps' function" << endl;			
				return FFEA_ERROR;
			}
			
			if(load_kinetic_states(states, i) == FFEA_ERROR) {
				FFEA_error_text();
				cout << "Problem reading kinetic states in 'read_kinetic_states' function" << endl;			
				return FFEA_ERROR;
			}
			
			if(load_kinetic_rates(rates, i) == FFEA_ERROR) {
				FFEA_error_text();
				cout << "Problem reading kinetic rates in 'read_kinetic_rates' function" << endl;			
				return FFEA_ERROR;
			}
			cout << "\t...done!" << endl;
		
			cout << "\tBuilding extra maps for energy comparison..." << endl;
			if(build_kinetic_identity_maps() == FFEA_ERROR) {
				FFEA_error_text();
				cout << "Problem reading kinetic maps in 'build_kinetic_identity_maps' function" << endl;			
				return FFEA_ERROR;
			}

			cout << "\t...done!" << endl;
		}
		
		// Clear blob vector and other vectors for next round
		motion_state.clear();
		nodes.clear();
		topology.clear();
		surface.clear();
		material.clear();
		stokes.clear();
		vdw.clear();
		binding.clear();
		pin.clear();
		maps.clear();
                beads.clear();
		scale = 1;
		solver = FFEA_NOMASS_CG_SOLVER;
		map_vector.clear();
		kinetics_vector.clear();
		blob_vector.clear();
		centroid = NULL;
		velocity = NULL;
		rotation = NULL;
		set_scale = 0;
		set_solver = 0;
	}

	cout << "\t...done!" << endl;


	return FFEA_OK;
}

int World::load_kinetic_maps(vector<string> map_fnames, vector<int> map_from, vector<int> map_to, int blob_index) {
	
	int MAX_BUF_SIZE = 255;
	char buf[MAX_BUF_SIZE];
	unsigned int i, j, num_rows, num_entries;
	string buf_string;
	vector<string> string_vec;
	for(i = 0; i < map_fnames.size(); ++i) {

		cout << "\t\tReading map " << i << ", " << map_fnames.at(i) << ": from conformation " << map_from.at(i) << " to " << map_to.at(i) << endl;
		ifstream fin;
		fin.open(map_fnames.at(i).c_str());
		
		// Check if sparse or dense
		fin.getline(buf, MAX_BUF_SIZE);
		buf_string = string(buf);
		boost::trim(buf_string);
		if(buf_string != "FFEA Kinetic Conformation Mapping File (Sparse)") {
			FFEA_error_text();
			cout << "In " << map_fnames.at(i) << ", expected 'FFEA Kinetic Conformation Mapping File (Sparse)'" << endl;
			cout << "but got " << buf_string << endl;
			return FFEA_ERROR;
		}

		// Get nodes to and from
		fin.getline(buf, MAX_BUF_SIZE);
		fin.getline(buf, MAX_BUF_SIZE);
		buf_string = string(buf);
		boost::split(string_vec, buf_string, boost::is_space());
		num_rows = atoi(string_vec.at(string_vec.size() - 1).c_str());
		
		// Get num_entries
		fin.getline(buf, MAX_BUF_SIZE);
		buf_string = string(buf);
		boost::split(string_vec, buf_string, boost::is_space());
		num_entries = atoi(string_vec.at(string_vec.size() - 1).c_str());

		// Create some memory for the matrix stuff
		scalar *entries = new scalar[num_entries];
		int *key = new int[num_rows + 1];
		int *col_index = new int[num_entries];

		// Get 'map:'
		fin.getline(buf, MAX_BUF_SIZE);
		
		// Read matrix
		// 'entries -'
		fin >> buf_string;
		fin >> buf_string;
		for(j = 0; j < num_entries; j++) {
			fin >> buf_string;
			entries[j] = atof(buf_string.c_str());
		}

		// 'key -'
		fin >> buf_string;
		fin >> buf_string;
		for(j = 0; j < num_rows + 1; j++) {
			fin >> buf_string;
			key[j] = atoi(buf_string.c_str());
		}

		// 'columns -'
		fin >> buf_string;
		fin >> buf_string;
		for(j = 0; j < num_entries; j++) {
			fin >> buf_string;
			col_index[j] = atoi(buf_string.c_str());
		}

		// Close file
		fin.close();

		// Create sparse matrix
		kinetic_map_array[blob_index][map_from[i]][map_to[i]].init(num_rows, num_entries, entries, key, col_index);
	}

	return FFEA_OK;
}

int World::build_kinetic_identity_maps() {
	
	int i, j, k;

	// For each blob, build map_ij*map_ji and map_ji*map_ij so we can compare energies using only the conserved modes. Well clever this, Oliver Harlen's idea.
	// He didn't write this though!
	
	for(i = 0; i < params.num_blobs; ++i) {
		for(j = 0; j < params.num_conformations[i]; ++j) {
			for(k = j + 1; k < params.num_conformations[i]; ++k) {

				kinetic_double_map_array[i][j][k] = kinetic_map_array[i][k][j].apply(&kinetic_map_array[i][j][k]);
				kinetic_double_map_array[i][k][j] = kinetic_map_array[i][j][k].apply(&kinetic_map_array[i][k][j]);
			}
		}
	}
	return FFEA_OK;
}

int World::rescale_kinetic_rates(scalar ***rates) {

	return FFEA_OK;
	int i, j, current_state, case_index;
	scalar total_prob, current_energy, final_energy;
	for(i = 0; i < params.num_blobs; ++i) {
		
		// Get current energy
		current_state = active_state_index[i];
		current_energy = active_blob_array[i]->calculate_strain_energy();
		//cout << i << " " << current_energy / (params.kT * (3 * active_blob_array[i]->get_num_linear_nodes() - 6)) << endl;

		total_prob = 0.0;
		for(j = 0; j < params.num_states[i]; ++j) {

			// Calculate remainder at end
			if(j == current_state) {
				continue;
			}

			// Rates change depending on situation
			if(kinetic_state[i][current_state].bound == kinetic_state[i][j].bound && kinetic_state[i][current_state].conformation_index == kinetic_state[i][j].conformation_index)
			{
				// Null transform
				case_index = 0;

			} else if(kinetic_state[i][current_state].bound == 0 && kinetic_state[i][j].bound == 1) {
				
				// Unbound to bound
				case_index = 1;

			} else if(kinetic_state[i][current_state].bound == 1 && kinetic_state[i][j].bound == 0) {
		
				// Bound to unbound
				case_index = 2;

			} else if (kinetic_state[i][current_state].conformation_index != kinetic_state[i][j].conformation_index) {
			
				// Conformation change
				case_index = 3;
			} else {
				FFEA_ERROR_MESSG("Error when trying to rescale kinetic rate: Blob %d from state %d to state %d\n", i, current_state, j);
				return FFEA_ERROR;
			}

			switch(case_index) {
				
				// Null
				case(0):
					break;

				// Binding
				case(1):

					/*for(k = 0; k < num_binding_sites; ++k) {

						// Binding site on our molecule
						if(kinetic_binding_sites[k].blob_index == i && kinetic_binding_sites[k].conf_index == j) {
							kinetic_binding_sites[k].calc_centroid();
							radius[0] = sqrt(kinetic_binding_sites[k].area / 3.1415926);
							for(l = 0; l < num_binding_sites; ++l) {
								if(l == k) {
									continue;
								}
								kinetic_binding_sites[l].calc_centroid();
								radius[1] = sqrt(kinetic_binding_sites[l].area / 3.1415926);
								distance.x = kinetic_binding_sites[k].centroid.x - kinetic_binding_sites[l].centroid.x;
								distance.y = kinetic_binding_sites[k].centroid.y - kinetic_binding_sites[l].centroid.y;
								distance.z = kinetic_binding_sites[k].centroid.z - kinetic_binding_sites[l].centroid.z;
								
								// Successful binding!
								if(mag(&distance) < radius[0] + radius[1]) {

								}

								
							}
						}	
					}*/				
					break;
		
				// Unbinding
				case(2):

					
					break;

				// Conformation change
				case(3):
				
					
					break;
			}

			total_prob += rates[i][current_state][j];
			
			
		}

		// Prob of staying put
		rates[i][current_state][current_state] = 1 - total_prob;
	}

	return FFEA_OK;
}

int World::load_kinetic_states(string states_fname, int blob_index) {
	
	int i, MAX_BUF_SIZE = 255, num_kinetic_states;
	int conf_index, bound_state, binding_site_type_from, binding_site_type_to;
	char buf[MAX_BUF_SIZE];
	string string_buf;

	cout << "\t\tReading states file '" << states_fname << "'" << endl;
	ifstream fin;
	fin.open(states_fname.c_str());

	// Check filetype
	fin.getline(buf, MAX_BUF_SIZE);
	string_buf = string(buf);
	boost::trim(string_buf);
	if(string_buf != "ffea kinetic states file") {
		FFEA_error_text();
		cout << "In '" << states_fname << "', expected 'ffea kinetic states file' but got " << buf << endl;
		return FFEA_ERROR;
	}
	
	// Get num_states
	fin >> buf >> num_kinetic_states;

	// Check against params
	if(params.num_states[blob_index] != num_kinetic_states) {
		FFEA_error_text();
		cout << "Read 'num_states " << num_kinetic_states << "' from " << states_fname << " but expected " << params.num_states[blob_index] << " from .ffea script." << endl;
		return FFEA_ERROR;
	}

	// Create states objects
	kinetic_state[blob_index] = new KineticState[num_kinetic_states];

	// Get states line
	fin.getline(buf, MAX_BUF_SIZE);
	fin.getline(buf, MAX_BUF_SIZE);

	// Initialise states
	for(i = 0; i < num_kinetic_states; ++i) {

		// Get conformation index
		fin >> conf_index;
		if(conf_index >= params.num_conformations[blob_index]) {
			FFEA_error_text();
			cout << "In kinetic states file '" << states_fname << "', state " << i << " references a conformation index which is out of the range of the number of conformation defined in the script" << endl;
			return FFEA_ERROR;
		}

		// Get bound state
		fin >> string_buf;
		if(string_buf == "u" or string_buf == "unbound") {
			bound_state = 0;

		} else if (string_buf == "b" or string_buf == "bound") {
			bound_state = 1;

		} else {
			FFEA_error_text();
			cout << "Expected bound/unbound, got " << string_buf << endl;
			return FFEA_ERROR;
		}

		// Get site types in binding event
		fin >> binding_site_type_from;
		fin >> binding_site_type_to;
		if(binding_site_type_from >= binding_matrix.num_interaction_types || binding_site_type_to >= binding_matrix.num_interaction_types) {
			FFEA_error_text();
			cout << "In kinetic states file '" << states_fname << "', state " << i << " references a state type which is out of the range of the number of types defined in 'binding_params'" << endl;
			return FFEA_ERROR; 
		}

		// Initialise state
		kinetic_state[blob_index][i].init(conf_index, bound_state, binding_site_type_from, binding_site_type_to);
	}
	return FFEA_OK;
}

int World::load_kinetic_rates(string rates_fname, int blob_index) {
	
	int i, j, MAX_BUF_SIZE = 255, num_kinetic_states;
	char buf[MAX_BUF_SIZE];
	string string_buf;

	cout << "\t\tReading rates file '" << rates_fname << "'" << endl;
	ifstream fin;
	fin.open(rates_fname.c_str());

	// Check filetype
	fin.getline(buf, MAX_BUF_SIZE);
	string_buf = string(buf);
	boost::trim(string_buf);
	if(string_buf != "ffea kinetic rates file") {
		FFEA_error_text();
		cout << "In '" << rates_fname << "', expected 'ffea kinetic rates file' but got " << buf << endl;
		return FFEA_ERROR;
	}
	
	// Get num_states
	fin >> buf >> num_kinetic_states;
	
	// Check against params
	if(params.num_states[blob_index] != num_kinetic_states) {
		FFEA_error_text();
		cout << "Read 'num_states " << num_kinetic_states << "' from " << rates_fname << " but expected " << params.num_states[blob_index] << " from .ffea script." << endl;
		return FFEA_ERROR;
	}

	// Create rates matrix
	kinetic_rate[blob_index] = new scalar*[num_kinetic_states];
	for(i = 0; i < num_kinetic_states; ++i) {
		kinetic_rate[blob_index][i] = new scalar[num_kinetic_states];	
	}

	// Get rates line
	fin.getline(buf, MAX_BUF_SIZE);
	fin.getline(buf, MAX_BUF_SIZE);
	double total_prob;
	for(i = 0; i < num_kinetic_states; ++i) {
		total_prob = 0.0;
		for(j = 0; j < num_kinetic_states; ++j) {
			fin >> kinetic_rate[blob_index][i][j];

			// Change to probabilities
			kinetic_rate[blob_index][i][j] *= params.dt * params.kinetics_update;
			if(i != j) {
				total_prob += kinetic_rate[blob_index][i][j];
			}
		}
		
		// Prob of not switching (for completion)
		if(total_prob > 1) {
			FFEA_error_text();
			cout << "P(switch_state in kinetic update period) = rate(switch_state)(Hz) * dt * kinetics_update" << endl;
			cout << "Due to the size of your rates, your timestep, and your kinetic_update value, the total probability of changing states each kinetic update period is greater than one." << endl;
			cout << "Best solution - Reduce 'kinetics_update' parameter" << endl;
			return FFEA_ERROR;
		}
		kinetic_rate[blob_index][i][i] = 1 - total_prob;
	}

	return FFEA_OK;
}

/* */
void World::get_system_CoM(vector3 *system_CoM) {
    system_CoM->x = 0;
    system_CoM->y = 0;
    system_CoM->z = 0;
    scalar total_mass = 0;
    for (int i = 0; i < params.num_blobs; i++) {
        vector3 com;
        active_blob_array[i]->get_CoM(&com);
        system_CoM->x += com.x * active_blob_array[i]->get_mass();
        system_CoM->y += com.y * active_blob_array[i]->get_mass();
        system_CoM->z += com.z * active_blob_array[i]->get_mass();

        total_mass += active_blob_array[i]->get_mass();
    }
    system_CoM->x /= total_mass;
    system_CoM->y /= total_mass;
    system_CoM->z /= total_mass;
}

/* */
void World::get_system_centroid(vector3 *centroid) {
    centroid->x = 0;
    centroid->y = 0;
    centroid->z = 0;
    scalar total_num_nodes = 0;
    for (int i = 0; i < params.num_blobs; i++) {
        vector3 cen;
        active_blob_array[i]->get_centroid(&cen);
        centroid->x += cen.x * active_blob_array[i]->get_num_nodes();
        centroid->y += cen.y * active_blob_array[i]->get_num_nodes();
        centroid->z += cen.z * active_blob_array[i]->get_num_nodes();

        total_num_nodes += active_blob_array[i]->get_num_nodes();
    }
    centroid->x /= total_num_nodes;
    centroid->y /= total_num_nodes;
    centroid->z /= total_num_nodes;
}

void World::get_system_dimensions(vector3 *dimension) {
	dimension->x = 0;
	dimension->y = 0;
	dimension->z = 0;
	
	vector3 min, max;
	min.x = INFINITY;
	min.y = INFINITY;
	min.z = INFINITY;
	max.x = -1 * INFINITY;
	max.y = -1 * INFINITY;
	max.z = -1 * INFINITY;
	
	vector3 blob_min, blob_max;
	for(int i = 0; i < params.num_blobs; i++) {
		active_blob_array[i]->get_min_max(&blob_min, &blob_max);
		if(blob_min.x < min.x) {
			min.x = blob_min.x;
		}
		if(blob_min.y < min.y) {
			min.y = blob_min.y;
		}
		if(blob_min.z < min.z) {
			min.z = blob_min.z;
		}

		if(blob_max.x > max.x) {
			max.x = blob_max.x;
		}

		if(blob_max.y > max.y) {
			max.y = blob_max.y;
		}

		if(blob_max.z > max.z) {
			max.z = blob_max.z;
		}
	}

	dimension->x = max.x - min.x;
	dimension->y = max.y - min.y;
	dimension->z = max.z - min.z;
}

int World::get_num_blobs() {

	return params.num_blobs;
}

int World::load_springs(const char *fname) {

    int i;
    FILE *in = NULL;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the spring file
    if ((in = fopen(fname, "r")) == NULL) {
        FFEA_FILE_ERROR_MESSG(fname)
    }
    printf("\t\tReading in springs file: %s\n", fname);

    // first line should be the file type "ffea springs file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading first line of spring file\n")
    }
    if (strcmp(line, "ffea spring file\n") != 0) {
        fclose(in);
        FFEA_ERROR_MESSG("This is not a 'ffea spring file' (read '%s') \n", line)
    }

    // read in the number of springs in the file
    if (fscanf(in, "num_springs %d\n", &num_springs) != 1) {
        fclose(in);
        FFEA_ERROR_MESSG("Error reading number of springs\n")
    }
    printf("\t\t\tNumber of springs = %d\n", num_springs);

    // Allocate memory for springs
    spring_array = new Spring[num_springs];

    // Read in next line
    fscanf(in,"springs:\n");
    for (i = 0; i < num_springs; ++i) {
       if (fscanf(in, "%d %d %d %d %d %d %lf %lf\n", &spring_array[i].blob_index[0], &spring_array[i].conformation_index[0], &spring_array[i].node_index[0], &spring_array[i].blob_index[1], &spring_array[i].conformation_index[1], &spring_array[i].node_index[1], &spring_array[i].k, &spring_array[i].l) != 8) {
            FFEA_error_text();
            printf("Problem reading spring data from %s. Format is:\n\n", fname);
            printf("ffea spring file\nnum_springs ?\n");
            printf("blob_index_0 conformation_index_0 node_index_0 blob_index_1 conformation_index_1 node_index_1 k l\n\n");
            return FFEA_ERROR;
        }

	// Error checking
	for(int j = 0; j < 2; ++j) {
		if(spring_array[i].blob_index[j] >= params.num_blobs || spring_array[i].blob_index[j] < 0) {
			FFEA_error_text();
        	    	printf("In spring %d, blob index %d is out of bounds given the number of blobs defined (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, params.num_blobs);
			return FFEA_ERROR;
		}
		if(spring_array[i].conformation_index[j] >= params.num_conformations[spring_array[i].blob_index[j]] || spring_array[i].conformation_index[j] < 0) {
			FFEA_error_text();
        	    	printf("In spring %d, conformation index %d is out of bounds given the number of conformations defined in blob %d (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, spring_array[i].blob_index[j], params.num_conformations[spring_array[i].blob_index[j]]);
			return FFEA_ERROR;
		}
		if(spring_array[i].node_index[j] >= blob_array[spring_array[i].blob_index[j]][spring_array[i].conformation_index[j]].get_num_nodes() || spring_array[i].node_index[j] < 0) {
			FFEA_error_text();
        	    	printf("In spring %d, node index %d is out of bounds given the number of nodes defined in blob %d, conformation %d (%d). Please fix. Remember, indexing starts at ZERO!\n", i, j, spring_array[i].blob_index[j], spring_array[i].conformation_index[j], blob_array[spring_array[i].blob_index[j]][spring_array[i].conformation_index[j]].get_num_nodes());
			return FFEA_ERROR;
		}
		if(spring_array[i].k < 0) {
			FFEA_error_text();
			printf("In spring %d, spring constant, %e, < 0. This is not going to end well for you...\n", i, spring_array[i].k);
			return FFEA_ERROR;
		}
		if(spring_array[i].l < 0) {
			FFEA_error_text();
			printf("In spring %d, spring equilibrium length, %e, < 0. Reverse node definitions for consistency.\n", i, spring_array[i].l);
			return FFEA_ERROR;
		}
	}
	if(spring_array[i].blob_index[0] == spring_array[i].blob_index[1] && spring_array[i].conformation_index[0] == spring_array[i].conformation_index[1] && spring_array[i].node_index[0] == spring_array[i].node_index[1]) {
			FFEA_error_text();
			printf("In spring %d, spring is connected to same node on same blob on same conformation. Will probably cause an force calculation error.\n", i);
			return FFEA_ERROR;
	}
	if(spring_array[i].blob_index[0] == spring_array[i].blob_index[1] && spring_array[i].conformation_index[0] != spring_array[i].conformation_index[1]) {
			FFEA_error_text();
			printf("In spring %d, spring is connected two conformations of the same blob (blob_index = %d). This cannot happen as conformations are mutually exclusive.\n", i, spring_array[i].blob_index[0]);
			return FFEA_ERROR;
	}
    }

    fclose(in);
    printf("\t\t\tRead %d springs from %s\n", i, fname);
    activate_springs();
    return 0;
}

void World::activate_springs() {
    for (int i = 0; i < num_springs; ++i) {
        if (spring_array[i].conformation_index[0] == active_blob_array[spring_array[i].blob_index[0]]->conformation_index && spring_array[i].conformation_index[1] == active_blob_array[spring_array[i].blob_index[1]]->conformation_index) {
            spring_array[i].am_i_active = true;
        } else {
	    spring_array[i].am_i_active = false;
	}
    }
}

void World::apply_springs() {
    scalar force_mag;
    vector3 n1, n0, force0, force1, sep, sep_norm;
    for (int i = 0; i < num_springs; ++i) {
        if (spring_array[i].am_i_active == true) {
            n1 = active_blob_array[spring_array[i].blob_index[1]]->get_node(spring_array[i].node_index[1]);
            n0 = active_blob_array[spring_array[i].blob_index[0]]->get_node(spring_array[i].node_index[0]);
            sep.x = n1.x - n0.x;
            sep.y = n1.y - n0.y;
            sep.z = n1.z - n0.z;
            sep_norm = normalise(&sep);
            force_mag = spring_array[i].k * (mag(&sep) - spring_array[i].l);
            force0.x = force_mag * sep_norm.x;
            force0.y = force_mag * sep_norm.y;
            force0.z = force_mag * sep_norm.z;

            force1.x = -1 * force_mag * sep_norm.x;
            force1.y = -1 * force_mag * sep_norm.y;
            force1.z = -1 * force_mag * sep_norm.z;
            active_blob_array[spring_array[i].blob_index[0]]->add_force_to_node(force0, spring_array[i].node_index[0]);
            active_blob_array[spring_array[i].blob_index[1]]->add_force_to_node(force1, spring_array[i].node_index[1]);
        }
    }
    return;
}

int World::get_next_script_tag(FILE *in, char *buf) {
    if (fscanf(in, "%*[^<]") != 0) {
        printf("White space removal error in get_next_script_tag(). Something odd has happened: this error should never occur...\n");
    }
    if (fscanf(in, "<%255[^>]>", buf) != 1) {
        FFEA_error_text();
        printf("Error reading tag in script file.\n");
        return FFEA_ERROR;
    }

    return FFEA_OK;
}

void World::apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N) {
    int i, j;
    for (i = 0; i < N; i++) {
        y[i] = 0;
        for (j = 0; j < N; j++)
            y[i] += M[i * N + j] * x[j];
    }
}

void World::do_es() {
    printf("Building BEM matrices\n");
    PB_solver.build_BEM_matrices();


    //	PB_solver.print_matrices();

    // Build the poisson matrices for each blob
    printf("Building Poisson matrices\n");
    for (int i = 0; i < params.num_blobs; i++) {
        active_blob_array[i]->build_poisson_matrices();
    }

    printf("Solving\n");
    for (int dual = 0; dual < 30; dual++) {


        // Perform Poisson solve step (K phi + E J_Gamma = rho)
        // Obtain the resulting J_Gamma
        int master_index = 0;
        for (int i = 0; i < params.num_blobs; i++) {
            active_blob_array[i]->solve_poisson(&phi_Gamma[master_index], &J_Gamma[master_index]);
            master_index += active_blob_array[i]->get_num_faces();
        }


        // Apply -D matrix to J_Gamma vector (work_vec = -D J_gamma)
        PB_solver.get_D()->apply(J_Gamma, work_vec);
        for (int i = 0; i < total_num_surface_faces; i++) {
            work_vec[i] *= -1;
        }

        // Solve for C matrix (C phi_Gamma = work_vec)
        nonsymmetric_solver.solve(PB_solver.get_C(), phi_Gamma, work_vec);

        scalar sum = 0.0, sumj = 0.0;
        for (int i = 0; i < total_num_surface_faces; i++) {
            sum += phi_Gamma[i];
            sumj += J_Gamma[i];
        }
        printf("\n");

        printf("<yophi> = %e <J> = %e\n", sum / total_num_surface_faces, sumj / total_num_surface_faces);

        //		for(i = 0; i < total_num_surface_faces; i++) {
        //			printf("<JPBSN> %d %e\n", dual, J_Gamma[i]);
        //		}

    }

    //	scalar sum = 0.0, sumj = 0.0;
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		sum += phi_Gamma[i];
    //		sumj += J_Gamma[i];
    //	}
    //	printf("<phi> = %e <J> = %e\n", sum/total_num_surface_faces, sumj/total_num_surface_faces);
    //
    //	for(i = 0; i < total_num_surface_faces; i++) {
    //		printf("<J> %e\n", J_Gamma[i]);
    //	}

    //	blob_array[0].print_phi();
}

void World::make_trajectory_from_eigenvector(int blob_index, int mode_index, Eigen::VectorXd evec, double step) {

	// Get a filename
	int i, j;
	double dx;
	vector<string> all;
	string traj_out_fname, base, ext;
	ostringstream bi, mi;
	bi << blob_index;
	mi << mode_index;
	boost::split(all, params.trajectory_out_fname, boost::is_any_of("."));
	ext = "." + all.at(all.size() - 1);
	base = boost::erase_last_copy(string(params.trajectory_out_fname), ext);
	traj_out_fname = base + "_blob" + bi.str() + "mode" + mi.str() + ext;
	int from_index = 0;
	int to_index = 0;

	// Convert weird eigen thing into a nice vector3 list
	vector3 node_trans[active_blob_array[blob_index]->get_num_linear_nodes()];

	// Open file
	FILE *fout;
	fout = fopen(traj_out_fname.c_str(), "w");

	// Header Stuff
	fprintf(fout, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs 1\nNumber of Conformations 1\nBlob 0:	Conformation 0 Nodes %d\n\n*\n", active_blob_array[blob_index]->get_num_nodes());

	for(i = 0; i < 21; ++i) {
		
		/* Build a frame */

		// Get eigenvector multiplier
		if(i == 0) {
			dx = 0.0;
		} else if(i > 0 && i < 5) {
			dx = step;
		} else if (i > 5 && i <= 15) {
			dx = -step;
		} else {
			dx = step;
		}

		// Get some node translations
		int num_linear_nodes = active_blob_array[blob_index]->get_num_linear_nodes();
		for(j = 0; j < num_linear_nodes; ++j) {
			node_trans[j].x = evec[0 * num_linear_nodes + j] * dx;
			node_trans[j].y = evec[1 * num_linear_nodes + j] * dx;
			node_trans[j].z = evec[2 * num_linear_nodes + j] * dx;
		}

		// Translate all the nodes
		active_blob_array[blob_index]->translate_linear(node_trans);
		fprintf(fout, "Blob 0, Conformation 0, step %d\n", i);
		active_blob_array[blob_index]->write_nodes_to_file(fout);
		fprintf(fout, "*\n");
		print_trajectory_conformation_changes(fout, i, &from_index, &to_index);
	
	}
	fclose(fout);
}

void World::print_trajectory_and_measurement_files(int step, double wtime) {
    if ((step - 1) % (params.check * 10) != 0) {
        printf("step = %d\n", step);
    } else {
        printf("step = %d (simulation time = %.2es, wall clock time = %.2e hrs)\n", step, (scalar) step * params.dt, (omp_get_wtime() - wtime) / 3600.0);
    }

    vector3 system_CoM;
    get_system_CoM(&system_CoM);

    // Write traj and meas data
    if (measurement_out[params.num_blobs] != NULL) {
        fprintf(measurement_out[params.num_blobs], "%d\t", step);
    }

    for (int i = 0; i < params.num_blobs; i++) {

        // Write the node data for this blob
        fprintf(trajectory_out, "Blob %d, Conformation %d, step %d\n", i, active_conformation_index[i], step);
        active_blob_array[i]->write_nodes_to_file(trajectory_out);

        // Write the measurement data for this blob
        active_blob_array[i]->make_measurements(measurement_out[i], step, &system_CoM);

        // Output interblob_vdw info
	for (int j = i + 1; j < params.num_blobs; ++j) {
            active_blob_array[i]->calculate_vdw_bb_interaction_with_another_blob(measurement_out[params.num_blobs], j);
        }
    }

    fprintf(measurement_out[params.num_blobs], "\n");

    // Mark completed end of step with an asterisk (so that the restart code will know if this is a fully written step or if it was cut off half way through due to interrupt)
    fprintf(trajectory_out, "*\n");

    // Force all output in buffers to be written to the output files now
    fflush(trajectory_out);
    for (int i = 0; i < num_blobs + 1; ++i) {
        fflush(measurement_out[i]);
    }
}

void World::print_trajectory_conformation_changes(FILE *fout, int step, int *from_index, int *to_index) {

	// Check input
	int *to;
	int *from;
	if(to_index == NULL || from_index == NULL) {
		to = new int[params.num_blobs];
		from = new int[params.num_blobs];
		for(int i = 0; i < params.num_blobs; ++i) {
			to[i] = active_blob_array[i]->conformation_index;
			from[i] = active_blob_array[i]->conformation_index;
		}
	} else {
		to = to_index;
		from = from_index;
	}

	// Inform whoever is watching of changes
	if(params.calc_kinetics == 1 && (step - 1) % params.kinetics_update == 0) {
		printf("Conformation Changes:\n");
	}	
	fprintf(fout, "Conformation Changes:\n");
	for(int i = 0; i < params.num_blobs; ++i) {
		if(params.calc_kinetics == 1 && (step - 1) % params.kinetics_update == 0) {
			printf("\tBlob %d - Conformation %d -> Conformation %d\n", i, from[i], to[i]);	
		}

		// Print to file
		fprintf(fout, "Blob %d: Conformation %d -> Conformation %d\n", i, from[i], to[i]);
	}
	fprintf(fout, "*\n");

	// Force print in case of ctrl + c stop
	fflush(fout);

	// Free data
	if(to_index == NULL || from_index == NULL) {
		delete[] to;
		delete[] from;
	}
}
void World::print_static_trajectory(int step, double wtime, int blob_index) {
    printf("Printing single trajectory of Blob %d for viewer\n", blob_index);
    // Write the node data for this blob
    fprintf(trajectory_out, "Blob %d, step %d\n", blob_index, step);
    active_blob_array[blob_index]->write_nodes_to_file(trajectory_out);
}

// Well done for reading this far! Hope this makes you smile. People should smile more

/*
quu..__
 $$$b  `---.__
  "$$b        `--.                          ___.---uuudP
   `$$b           `.__.------.__     __.---'      $$$$"              .
     "$b          -'            `-.-'            $$$"              .'|
       ".                                       d$"             _.'  |
         `.   /                              ..."             .'     |
           `./                           ..::-'            _.'       |
            /                         .:::-'            .-'         .'
           :                          ::''\          _.'            |
          .' .-.             .-.           `.      .'               |
          : /'$$|           .@"$\           `.   .'              _.-'
         .'|$u$$|          |$$,$$|           |  <            _.-'
         | `:$$:'          :$$$$$:           `.  `.       .-'
         :                  `"--'             |    `-.     \
        :##.       ==             .###.       `.      `.    `\
        |##:                      :###:        |        >     >
        |#'     `..'`..'          `###'        x:      /     /
         \                                   xXX|     /    ./
          \                                xXXX'|    /   ./
          /`-.                                  `.  /   /
         :    `-  ...........,                   | /  .'
         |         ``:::::::'       .            |<    `.
         |             ```          |           x| \ `.:``.
         |                         .'    /'   xXX|  `:`M`M':.
         |    |                    ;    /:' xXXX'|  -'MMMMM:'
         `.  .'                   :    /:'       |-'MMMM.-'
          |  |                   .'   /'        .'MMM.-'
          `'`'                   :  ,'          |MMM<
            |                     `'            |tbap\
             \                                  :MM.-'
              \                 |              .''
               \.               `.            /
                /     .:::::::.. :           /
               |     .:::::::::::`.         /
               |   .:::------------\       /
              /   .''               >::'  /
              `',:                 :    .'
*/
