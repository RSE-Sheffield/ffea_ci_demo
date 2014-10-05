#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <omp.h>

#include "MersenneTwister.h"
#include "NearestNeighbourLinkedListCube.h"
#include "BEM_Poisson_Boltzmann.h"
#include "BiCGSTAB_solver.h"
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"
#include "Face.h"
#include "Blob.h"
#include "World.h"
#include "VdW_solver.h"
#include "LJ_matrix.h"

class World
{
	public:
		World()
		{
			// Initialise everything to zero
			blob_array = NULL;
			num_threads = 0;
			rng = NULL;
			phi_Gamma = NULL;
			total_num_surface_faces = 0;
			box_dim.x = 0;
			box_dim.y = 0;
			box_dim.z = 0;
			step_initial = 0;
		}

		~World()
		{
			delete[] rng;
			rng = NULL;
			num_threads = 0;

			delete[] blob_array;
			blob_array = NULL;

			delete[] phi_Gamma;
			phi_Gamma = NULL;

			total_num_surface_faces = 0;

			box_dim.x = 0;
			box_dim.y = 0;
			box_dim.z = 0;
			step_initial = 0;
		}

		/* */
		int init(const char *FFEA_script_filename)
		{
			int i, j;

			// Open the ffea script file
			FILE *script_file;
			if((script_file = fopen(FFEA_script_filename, "r")) == NULL) {
				FFEA_FILE_ERROR_MESSG(FFEA_script_filename);
			}

			const int buffer_size = 101;
			char buf[buffer_size];

			printf("Parsing script file '%s'...\n", FFEA_script_filename);

			// Read first tag
			if(get_next_script_tag(script_file, buf) == FFEA_ERROR) {
				FFEA_error_text();
				printf("Error when reading first tag in %s\n", FFEA_script_filename);
				fclose(script_file);
				return FFEA_ERROR;
			}

			// 1st block should be <param> block
			if(strcmp(buf, "param") != 0) {
				FFEA_error_text();
				printf("In script file, '%s': First block in script file should be <param> block, not '<%s>'\n", FFEA_script_filename, buf);
				fclose(script_file);
				return FFEA_ERROR;
			}

			// Parse each tag in the <param> block
			printf("Entering <param> block...\n");
			for(;;) {
				if(get_next_script_tag(script_file, buf) == FFEA_ERROR) {
					FFEA_error_text();
					printf("\tError when reading tag in <param> block\n");
					fclose(script_file);
					return FFEA_ERROR;
				}
				if(feof(script_file)) {
					FFEA_error_text();
					printf("\tReached end of file before end of <param> block\n");
					fclose(script_file);
					return FFEA_ERROR;
				}

				if(strcmp(buf, "/param") == 0) {
					printf("Exiting <param> block\n");
					break;
				}

				if(params.parse_param_assignment(buf) == FFEA_ERROR) {
					FFEA_error_text();
					printf("Parameter assignment failed.\n");
					fclose(script_file);
					return FFEA_ERROR;
				}
			}

			// Check if the simulation parameters are reasonable
			if(params.validate() == FFEA_ERROR) {
				FFEA_error_text();
				printf("There are some required parameters that are either unset or have bad values.\n");
				fclose(script_file);
				return FFEA_ERROR;
			}

			// Load the vdw forcefield params matrix
			if(lj_matrix.init(params.vdw_params_fname) == FFEA_ERROR) {
				FFEA_ERROR_MESSG("Error when reading from vdw forcefeild params file.\n")
			}

			// detect how many threads we have for openmp
			int tid;
			#pragma omp parallel default(none) private(tid)
        		{
                		tid = omp_get_thread_num();
                		if(tid == 0) {
					num_threads = omp_get_num_threads();
                        		printf("Number of threads detected: %d\n", num_threads);
                		}
        		}

        		// We need one rng for each thread, to avoid concurrency problems,
        		// so generate an array of instances of mersenne-twister rngs.
			rng = new MTRand[num_threads];

			// Seed each rng differently
			for(i = 0; i < num_threads; i++)
				rng[i].seed(params.rng_seed + i);


			// Read next main block tag
			if(get_next_script_tag(script_file, buf) == FFEA_ERROR) {
				FFEA_error_text();
				printf("Error when reading/looking for <system> tag in %s\n", FFEA_script_filename);
				fclose(script_file);
				return FFEA_ERROR;
			}

			// 2nd block should be <system> block
			if(strcmp(buf, "system") != 0) {
				FFEA_error_text();
				printf("In script file, '%s': Second block in script file should be <system> block, not '<%s>'\n", FFEA_script_filename, buf);
				fclose(script_file);
				return FFEA_ERROR;
			}

			// Parse system block
			printf("Entering <system> block...\n");
			if(read_and_build_system(script_file) == FFEA_ERROR) {
				FFEA_error_text();
				printf("\tError when reading <system> block. Aborting.\n");
				fclose(script_file);
				return FFEA_ERROR;
			}

			fclose(script_file);


		/*
			// Calculate the total number of surface elements in the entire system
			int total_num_surface_elements = 0;
			for(i = 0; i < params.num_blobs; i++)
				total_num_surface_elements += blob_array[i].get_num_surface_elements();
			printf("Total number of surface elements in system: %d\n", total_num_surface_elements);

			// Allocate memory for an NxNxN grid, with a 'pool' for the required number of surface elements
			printf("Allocating memory for nearest neighbour lookup grid...\n");
			if(surface_element_lookup.alloc(params.es_N, total_num_surface_elements) == FFEA_ERROR) {
				FFEA_error_text();
				printf("When allocating memory for nearest neighbour surface elements lookup grid\n");
				return FFEA_ERROR;
			}
			printf("...done\n");

			// Add all the surface elements from each Blob to the lookup pool
			printf("Adding all surface elements to nearest neighbour grid lookup pool\n");
			for(i = 0; i < params.num_blobs; i++) {
				for(j = 0; j < blob_array[i].get_num_surface_elements(); j++) {
					if(surface_element_lookup.add_to_pool(blob_array[i].get_element(j)) == FFEA_ERROR) {
						FFEA_error_text();
						printf("When attempting to add a surface element to the surface_element_lookup pool\n");
						return FFEA_ERROR;
					}
				}
			}
		*/

			measurement_out = new FILE *[params.num_blobs + 1];

			// If not restarting a previous simulation, create new trajectory and measurement files
			if(params.restart == 0) {
				// Open the trajectory output file for writing
				if((trajectory_out = fopen(params.trajectory_out_fname, "w")) == NULL) {
					FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
				}

				// Open the measurement output file for writing
				for(i = 0; i < params.num_blobs + 1; ++i) {
					if((measurement_out[i] = fopen(params.measurement_out_fname[i], "w")) == NULL) {
						FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
					}
				}

				// Open stress file for writing, maybe
				if(params.stress_out_fname_set == 1) {
					if(((stress_out = fopen(params.stress_out_fname, "w")) == NULL)) {
						FFEA_FILE_ERROR_MESSG(params.stress_out_fname)
					}					
				}

				// Print initial info stuff
				fprintf(trajectory_out, "FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\n", params.num_blobs);
				for(i = 0; i < params.num_blobs; ++i) {
					fprintf(trajectory_out, "Blob %d Nodes %d\t", i, blob_array[i].get_num_nodes());
				}
				fprintf(trajectory_out, "\n\n");

				// First line in trajectory data should be an asterisk (used to delimit different steps for easy seek-search in restart code)
				fprintf(trajectory_out, "*\n");

				// First line in measurements file should be a header explaining what quantities are in each column
				for(i = 0; i < params.num_blobs; ++i) {
					fprintf(measurement_out[i], "# step | KE | PE | CoM x | CoM y | CoM z | L_x | L_y | L_z | rmsd | vdw_area_%d_surface | vdw_force_%d_surface | vdw_energy_%d_surface\n", i, i, i);
					fflush(measurement_out[i]);
				}
				fprintf(measurement_out[params.num_blobs], "# step ");
				for(i = 0; i < params.num_blobs; ++i) {
					for(j = i + 1; j < params.num_blobs; ++j) {
						fprintf(measurement_out[params.num_blobs], "| vdw_area_%d_%d | vdw_force_%d_%d | vdw_energy_%d_%d ", i, j, i, j, i, j);
					}
				}
				fprintf(measurement_out[params.num_blobs], "\n");
				fflush(measurement_out[params.num_blobs]);

			} else {
				// Otherwise, seek backwards from the end of the trajectory file looking for '*' character (delimitter for snapshots)
				printf("Restarting from trajectory file %s\n", params.trajectory_out_fname);
				if((trajectory_out = fopen(params.trajectory_out_fname, "r")) == NULL) {
					FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
				}

				printf("Reverse searching for 2 asterisks (denoting a completely written snapshot)...\n");
				if(fseek(trajectory_out, 0, SEEK_END) != 0) {
					FFEA_ERROR_MESSG("Could not seek to end of file\n")
				}

				// Variable to store position of last asterisk in trajectory file (initialise it at end of file)
				off_t last_asterisk_pos = ftello(trajectory_out);

				int num_asterisks = 0;
				while(num_asterisks != 2) {
					if(fseek(trajectory_out, -2, SEEK_CUR) != 0) {
						perror(NULL);
						FFEA_ERROR_MESSG("Balls.\n")
					}
					char c = fgetc(trajectory_out);
					if(c == '*') {
						num_asterisks++;
						printf("Found %d\n", num_asterisks);

						// get the position in the file of this last asterisk
						if(num_asterisks == 1) {
							last_asterisk_pos = ftello(trajectory_out);
						}
					}
				}

				char c;
				if((c = fgetc(trajectory_out)) != '\n') {
					ungetc(c, trajectory_out);
				}

				printf("Loading Blob position and velocity data from last completely written snapshot \n");
				int blob_id;
				long long rstep;
				for(int b = 0; b < params.num_blobs; b++) {
					if(fscanf(trajectory_out, "Blob %d, step %lld\n", &blob_id, &rstep) != 2) {
						FFEA_ERROR_MESSG("Error reading header info for Blob %d\n", b)
					}
					if(blob_id != b) {
						FFEA_ERROR_MESSG("Mismatch in trajectory file - found blob id = %d, was expecting blob id = %d\n", blob_id, b)
					}
					printf("Loading node position, velocity and potential from restart trajectory file, for blob %d, step %lld\n", blob_id, rstep);
					if(blob_array[b].read_nodes_from_file(trajectory_out) == FFEA_ERROR) {
						FFEA_ERROR_MESSG("Error restarting blob %d\n", blob_id)
					}
				}
				step_initial = rstep;
				printf("...done. Simulation will commence from step %lld\n", step_initial);
				fclose(trajectory_out);

				// Truncate the trajectory file up to the point of the last asterisk (thereby erasing any half-written time steps that may occur after it)
				printf("Truncating the trajectory file to the last asterisk (to remove any half-written time steps)\n");
				if(truncate(params.trajectory_out_fname, last_asterisk_pos) != 0) {
					FFEA_ERROR_MESSG("Error when trying to truncate trajectory file %s\n", params.trajectory_out_fname)
				}

				// Open trajectory and measurment files for appending
				printf("Opening trajectory and measurement files for appending.\n");
				if((trajectory_out = fopen(params.trajectory_out_fname, "a")) == NULL) {
					FFEA_FILE_ERROR_MESSG(params.trajectory_out_fname)
				}
				for(i = 0; i < params.num_blobs + 1; ++i) {
					if((measurement_out[i] = fopen(params.measurement_out_fname[i], "a")) == NULL) {
						FFEA_FILE_ERROR_MESSG(params.measurement_out_fname[i])
					}
				}

				// Open stress file for appending
				if(params.stress_out_fname_set == 1) {
					if(((stress_out = fopen(params.stress_out_fname, "a")) == NULL)) {
						FFEA_FILE_ERROR_MESSG(params.stress_out_fname)
					}					
				}

				// Append a newline to the end of this truncated trajectory file (to replace the one that may or may not have been there)
				fprintf(trajectory_out, "\n");
			
				for(i = 0; i < params.num_blobs + 1; ++i) {
					fprintf(measurement_out[i], "#==RESTART==\n");
				}

			}



			// Initialise the Van der Waals solver
			box_dim.x = params.es_h * (1.0/params.kappa) * params.es_N_x;
			box_dim.y = params.es_h * (1.0/params.kappa) * params.es_N_y;
			box_dim.z = params.es_h * (1.0/params.kappa) * params.es_N_z;
			vdw_solver.init(&lookup, &box_dim, &lj_matrix);
			
			// Calculate the total number of vdw interacting faces in the entire system
			total_num_surface_faces = 0;
			for(i = 0; i < params.num_blobs; i++) {
				total_num_surface_faces += blob_array[i].get_num_faces();
			}
			printf("Total number of surface faces in system: %d\n", total_num_surface_faces);

			if(params.es_N_x > 0 && params.es_N_y > 0 && params.es_N_z > 0) {

				// Allocate memory for an NxNxN grid, with a 'pool' for the required number of surface faces
				printf("Allocating memory for nearest neighbour lookup grid...\n");
				if(lookup.alloc(params.es_N_x, params.es_N_y, params.es_N_z, total_num_surface_faces) == FFEA_ERROR) {
					FFEA_error_text();
					printf("When allocating memory for nearest neighbour lookup grid\n");
					return FFEA_ERROR;
				}
				printf("...done\n");

				printf("Box has volume %e cubic angstroms\n", (box_dim.x * box_dim.y * box_dim.z) * 1e30);

				// Add all the faces from each Blob to the lookup pool
				printf("Adding all faces to nearest neighbour grid lookup pool\n");
				for(i = 0; i < params.num_blobs; i++) {
					int num_faces_added = 0;
					for(j = 0; j < blob_array[i].get_num_faces(); j++) {
						Face *b_face = blob_array[i].get_face(j);
						if(b_face != NULL) {
							if(lookup.add_to_pool(b_face) == FFEA_ERROR) {
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
				if(params.calc_es == 1) {
					printf("Initialising Boundary Element Poisson Boltzmann Solver...\n");
					PB_solver.init(&lookup);
					PB_solver.set_kappa(params.kappa);

					// Initialise the nonsymmetric matrix solver
					nonsymmetric_solver.init(total_num_surface_faces, params.epsilon2, params.max_iterations_cg);

					// Allocate memory for surface potential vector
					phi_Gamma = new scalar[total_num_surface_faces];
					for(i = 0; i < total_num_surface_faces; i++)
						phi_Gamma[i] = 0;

					work_vec = new scalar[total_num_surface_faces];
					for(i = 0; i < total_num_surface_faces; i++)
						work_vec[i] = 0;

					J_Gamma = new scalar[total_num_surface_faces];
					for(i = 0; i < total_num_surface_faces; i++)
						J_Gamma[i] = 0;
				}
			}

			if(params.restart == 0) {
				// Carry out measurements on the system before carrying out any updates (if this is step 0)
				print_trajectory_and_measurement_files(0, 0);
			}

			#ifdef FFEA_PARALLEL_WITHIN_BLOB
				printf("Now ready to run with 'within-blob parallelisation' (FFEA_PARALLEL_WITHIN_BLOB) on %d threads.\n", num_threads);
			#endif

			#ifdef FFEA_PARALLEL_PER_BLOB
				printf("Now ready to run with 'per-blob parallelisation' (FFEA_PARALLEL_PER_BLOB) on %d threads.\n", num_threads);
			#endif

			return FFEA_OK;
		}

		/* */
		int run()
		{
			// Update entire World for num_steps time steps
			int es_count = 1;
			double wtime = omp_get_wtime();
			for(long long step = step_initial; step < params.num_steps; step++) {

				// Zero the force across all blobs
				#ifdef FFEA_PARALLEL_PER_BLOB
					#pragma omp parallel for default(none) schedule(runtime)
				#endif
				for(int i = 0; i < params.num_blobs; i++) {
					blob_array[i].zero_force();
					if(params.calc_vdw == 1) {
						blob_array[i].zero_vdw_bb_measurement_data();
					}
					if(params.sticky_wall_xz == 1) {
						blob_array[i].zero_vdw_xz_measurement_data();
					}
				}

				#ifdef FFEA_PARALLEL_PER_BLOB
					#pragma omp parallel for default(none) schedule(guided) shared(stderr)
				#endif
				for(int i = 0; i < params.num_blobs; i++) {

				// If blob centre of mass moves outside simulation box, apply PBC to it
					vector3 com;
					blob_array[i].get_centroid(&com);
					scalar dx = 0, dy = 0, dz = 0;
					int check_move = 0;

					if(com.x < 0) {
						if(params.wall_x_1 == WALL_TYPE_PBC) {
							dx += box_dim.x;
		//					printf("fuck\n");
							check_move = 1;
						}
					} else if(com.x > box_dim.x) {
						if(params.wall_x_2 == WALL_TYPE_PBC) {
							dx -= box_dim.x;
		//					printf("fuck\n");
							check_move = 1;
						}
					}
					if(com.y < 0) {
						if(params.wall_y_1 == WALL_TYPE_PBC) {
							dy += box_dim.y;
		//					printf("fuck\n");
							check_move = 1;
						}
					} else if(com.y > box_dim.y) {
						if(params.wall_y_2 == WALL_TYPE_PBC) {
							dy -= box_dim.y;
		//					printf("fuck\n");
							check_move = 1;
						}
					}
					if(com.z < 0) {
						if(params.wall_z_1 == WALL_TYPE_PBC) {
							dz += box_dim.z;
		//					printf("fuck\n");
							check_move = 1;
						}
					} else if(com.z > box_dim.z) {
						if(params.wall_z_2 == WALL_TYPE_PBC) {
							dz -= box_dim.z;
		//					printf("fuck\n");
							check_move = 1;
						}
					}
						if(check_move == 1) {
						blob_array[i].move(dx, dy, dz);
					}

					// If Blob is near a hard wall, prevent it from moving further into it
					blob_array[i].enforce_box_boundaries(&box_dim);
				}



				if(params.calc_es == 1 || params.calc_vdw == 1 || params.sticky_wall_xz == 1) {
					if(es_count == params.es_update) {


						#ifdef FFEA_PARALLEL_PER_BLOB
							#pragma omp parallel for default(none) schedule(guided)
						#endif
						for(int i = 0; i < params.num_blobs; i++) {
							blob_array[i].calc_centroids_and_normals_of_all_faces();
							blob_array[i].reset_all_faces();
						}


						// Attempt to place all faces in the nearest neighbour lookup table
						if(lookup.build_nearest_neighbour_lookup(params.es_h * (1.0/params.kappa)) == FFEA_ERROR) {
							FFEA_error_text();
							printf("When trying to place faces in nearest neighbour lookup table.\n");

							// attempt to print out the final (bad) time step
							printf("Dumping final step:\n");
							print_trajectory_and_measurement_files(step, wtime);

							return FFEA_ERROR;
						}

						if(params.calc_vdw == 1) {
							vdw_solver.solve();
						}
						if(params.sticky_wall_xz == 1) {
							vdw_solver.solve_sticky_wall(params.es_h * (1.0/params.kappa));
						}

						if(params.calc_es == 1) {
							do_es();
						}

						es_count = 1;
					} else
						es_count++;
				}

				// Update all Blobs in the World
				int fatal_errors = 0;
				#ifdef FFEA_PARALLEL_PER_BLOB
					#pragma omp parallel for default(none) shared(step, wtime) reduction(+: fatal_errors) schedule(runtime)
				#endif
				for(int i = 0; i < params.num_blobs; i++) {
					if(blob_array[i].update() == FFEA_ERROR) {
						FFEA_error_text();
						printf("A problem occurred when updating Blob %d on step %lld\n", i, step);
						printf("Simulation ran for %2f seconds (wall clock time) before error ocurred\n", (omp_get_wtime() - wtime));
						//return FFEA_ERROR;
						fatal_errors++;
					}
				}
				if(fatal_errors > 0) {
					FFEA_error_text();
					printf("Detected %d fatal errors in this system update. Exiting now...\n", fatal_errors);

					// attempt to print out the final (bad) time step
					printf("Dumping final step:\n");
					print_trajectory_and_measurement_files(step, wtime);

					return FFEA_ERROR;
				}

				if(step%params.check == 0) {
					print_trajectory_and_measurement_files(step + 1, wtime);
				}
			}
			printf("Time taken: %2f seconds\n", (omp_get_wtime() - wtime));

			return FFEA_OK;
		}

		/* */
		int read_and_build_system(FILE *in)
		{
			const int max_buf_size = 255;
			char buf[max_buf_size];
			char lvalue[max_buf_size];
			char rvalue[max_buf_size];

			char node_filename[max_buf_size];
			char topology_filename[max_buf_size];
			char surface_filename[max_buf_size];
			char material_params_filename[max_buf_size];
			char stokes_filename[max_buf_size];
			char vdw_filename[max_buf_size];
			char pin_filename[max_buf_size];

			int linear_solver = FFEA_ITERATIVE_SOLVER;
			int blob_state = FFEA_BLOB_IS_DYNAMIC;
			float com_x, com_y, com_z;
			float vel_x, vel_y, vel_z;

			int set_node_filename = 0;
			int set_topology_filename = 0;
			int set_surface_filename = 0;
			int set_material_params_filename = 0;
			int set_stokes_filename = 0;
			int set_vdw_filename = 0;
			int set_pin_filename = 0;
			int set_linear_solver = 0;
			int set_blob_state = 0;
			int set_centroid_pos = 0;
			int set_velocity = 0;

			scalar scale = 1;

			// Create the array of Blobs
			printf("\tCreating blob array\n");
			blob_array = new Blob[params.num_blobs];

			// Construct each blob defined in the <system> block
			int i, rv;
			for(i = 0; i < params.num_blobs; i++) {

				if(get_next_script_tag(in, buf) == FFEA_ERROR) {
					FFEA_error_text();
					printf("\tError reading tag in <system> block\n");
					return FFEA_ERROR;
				}

				// Check if we have prematurely reached the end of the <system> block
				if(strcmp(buf, "/system") == 0) {
					FFEA_error_text();
					printf("\tPremature end of <system> block. Expected %d blob declarations; only found %d.\n", params.num_blobs, i);
					return FFEA_ERROR;
				}

				// Only blocks should be <blob> blocks
				if(strcmp(buf, "blob") != 0) {
					FFEA_error_text();
					printf("\tFound unexpected '<%s>' tag. Only allowed tag in <system> block is <blob> tag.\n", buf);
					return FFEA_ERROR;
				}

				// Read all blob initialisation info
				printf("\tEntering <blob> block %d...\n", i);
				for(;;) {
					
					if(get_next_script_tag(in, buf) == FFEA_ERROR) {
						FFEA_error_text();
						printf("\t\tError reading tag in <blob> block\n");
						return FFEA_ERROR;
					}
					
					if(feof(in)) {
						FFEA_error_text();
						printf("\t\tReached end of file before end of <blob> block\n");
						return FFEA_ERROR;
					}
					
					if(strcmp(buf, "/system") == 0) {
						FFEA_error_text();
						printf("\t\t</system> came before </blob>\n");
						break;
					}
					if(strcmp(buf, "/blob") == 0) {
						printf("\tExiting <blob> block %d\n", i);
						break;
					}

					rv = sscanf(buf, "%100[^=]=%s", lvalue, rvalue);
					if(rv != 2) {
						FFEA_error_text();
						printf("\t\tError parsing blob parameter assignment, '%s'\n", buf);
						return FFEA_ERROR;
                        		}

					rv = sscanf(lvalue, "%s", lvalue);
					rv = sscanf(rvalue, "%s", rvalue);

					if(strcmp(lvalue, "nodes") == 0) {
						strcpy(node_filename, rvalue);
						printf("\t\tSetting blob %d, nodes input file = '%s'\n", i, node_filename);
						set_node_filename = 1;
					} else if(strcmp(lvalue, "topology") == 0) {
						strcpy(topology_filename, rvalue);
						printf("\t\tSetting blob %d, topology input file = '%s'\n", i, topology_filename);
						set_topology_filename = 1;
					} else if(strcmp(lvalue, "surface") == 0) {
						strcpy(surface_filename, rvalue);
						printf("\t\tSetting blob %d, surface input file = '%s'\n", i, surface_filename);
						set_surface_filename = 1;
					} else if(strcmp(lvalue, "material") == 0) {
						strcpy(material_params_filename, rvalue);
						printf("\t\tSetting blob %d, material parameters file = '%s'\n", i, material_params_filename);
						set_material_params_filename = 1;
					} else if(strcmp(lvalue, "stokes") == 0) {
						strcpy(stokes_filename, rvalue);
						printf("\t\tSetting blob %d, stokes file = '%s'\n", i, stokes_filename);
						set_stokes_filename = 1;
					} else if(strcmp(lvalue, "vdw") == 0) {
						strcpy(vdw_filename, rvalue);
						printf("\t\tSetting blob %d, vdw file = '%s'\n", i, vdw_filename);
						set_vdw_filename = 1;
					} else if(strcmp(lvalue, "pin") == 0) {
						strcpy(pin_filename, rvalue);
						printf("\t\tSetting blob %d, pin file = '%s'\n", i, pin_filename);
						set_pin_filename = 1;
					} else if(strcmp(lvalue, "solver") == 0) {
						if(strcmp(rvalue, "CG") == 0) {
							printf("\t\tSetting blob %d, linear solver = CG (Preconditioned Jacobi Conjugate Gradient)\n", i);
							linear_solver = FFEA_ITERATIVE_SOLVER;
						} else if(strcmp(rvalue, "direct") == 0) {
							printf("\t\tSetting blob %d, linear solver = direct (Sparse Forward/Backward Substitution)\n", i);
							linear_solver = FFEA_DIRECT_SOLVER;
						} else if(strcmp(rvalue, "masslumped") == 0) {
							printf("\t\tSetting blob %d, linear solver = masslumped (diagonal mass matrix)\n", i);
							linear_solver = FFEA_MASSLUMPED_SOLVER;
						} else if(strcmp(rvalue, "CG_nomass") == 0) {
							printf("\t\tSetting blob %d, linear solver = CG_nomass (Preconditioned Jacobi Conjugate Gradient, no mass within system)\n", i);
							linear_solver = FFEA_NOMASS_CG_SOLVER;
						} else {
							FFEA_error_text();
							printf("\t\tThere is no solver option '%s'. Please use 'CG', 'direct' or 'CG_nomass'.\n", rvalue);
							return FFEA_ERROR;
						}
						set_linear_solver = 1;
					} else if(strcmp(lvalue, "state") == 0) {
						if(strcmp(rvalue, "STATIC") == 0) {
							printf("\t\tSetting blob %d, state = STATIC (Fixed position, no dynamics simulated)\n", i);
							blob_state = FFEA_BLOB_IS_STATIC;
						} else if(strcmp(rvalue, "DYNAMIC") == 0) {
							printf("\t\tSetting blob %d, state = DYNAMIC (Dynamics will be simulated)\n", i);
							blob_state = FFEA_BLOB_IS_DYNAMIC;
						} else if(strcmp(rvalue, "FROZEN") == 0) {
							printf("\t\tSetting blob %d, state = FROZEN (Dynamics will not be simulated, but Blob positions still written to trajectory)\n", i);
							blob_state = FFEA_BLOB_IS_FROZEN;
						} else {
							FFEA_error_text();
							printf("\t\tThere is no state option '%s'. Please use 'STATIC', 'DYNAMIC' or 'FROZEN'.\n", rvalue);
							return FFEA_ERROR;
						}
						set_blob_state = 1;
					} else if(strcmp(lvalue, "centroid_pos") == 0) {
						if(sscanf(rvalue, "(%e,%e,%e)", &com_x, &com_y, &com_z) != 3) {
							FFEA_error_text();
							printf("\t\tCould not carry out centroid position assignment: rvalue '%s' is badly formed. Should have form (x,y,z)\n", rvalue);
							return FFEA_ERROR;
						} else
							printf("\t\tSetting blob %d, initial centroid (%e, %e, %e)\n", i, com_x, com_y, com_z);

						set_centroid_pos = 1;
					} else if(strcmp(lvalue, "velocity") == 0) {
						if(sscanf(rvalue, "(%e,%e,%e)", &vel_x, &vel_y, &vel_z) != 3) {
							FFEA_error_text();
							printf("\t\tCould not carry out velocity assignment: rvalue '%s' is badly formed. Should have form (velx,vely,velz)\n", rvalue);
							return FFEA_ERROR;
						} else
							printf("\t\tSetting blob %d, initial velocity (%e, %e, %e)\n", i, vel_x, vel_y, vel_z);

						set_velocity = 1;
					} else if(strcmp(lvalue, "scale") == 0) {
						scale = atof(rvalue);
						printf("\t\tSetting blob %d, scale factor = %e\n", i, scale);
					} else {
						FFEA_error_text();
						printf("\t\tError: '%s' is not a recognised lvalue\n", lvalue);
						printf("\t\tRecognised lvalues are:\n");
						printf("\t\tnodes\n\ttopology\n\tvdw\n\tmeasurement\n\tsolver\n\tCoM_position\n\tscale\n");
						return FFEA_ERROR;
					}
				}

				// First make sure we have all the information necessary for initialisation
				if(set_blob_state == 0) {
					printf("blob_state option unset. Using 'DYNAMIC' by default.\n");
					blob_state = FFEA_BLOB_IS_DYNAMIC;
				}

				if(set_node_filename == 0) {
					FFEA_ERROR_MESSG("\tIn blob %d: 'nodes' has not been specified\n", i);
				}
				if(set_topology_filename == 0) {
					if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
						FFEA_ERROR_MESSG("\tIn blob %d: 'topology' has not been specified\n", i);
					}
				}
				if(set_surface_filename == 0) {
					FFEA_error_text();
					printf("\tIn blob %d: 'surface' has not been specified\n", i);
					return FFEA_ERROR;
				}

				if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
					if(set_material_params_filename == 0) {
						FFEA_ERROR_MESSG("\tIn blob %d: material params file, 'material', has not been specified\n", i)
					}
					if(set_stokes_filename == 0) {
						FFEA_ERROR_MESSG("\tIn blob %d: stokes radii file, 'stokes', has not been specified\n", i)
					}
					if(set_linear_solver == 0) {
						FFEA_ERROR_MESSG("\tIn blob %d: 'solver' has not been specified\n", i);
					}
				}
				if(set_vdw_filename == 0) {
					FFEA_ERROR_MESSG("\tIn blob %d: vdw parameters file, 'vdw', has not been specified\n", i)
				}
				if(set_pin_filename == 0) {
					if(blob_state == FFEA_BLOB_IS_DYNAMIC) {
						FFEA_ERROR_MESSG("\tIn blob %d: pinned nodes file, 'pin', has not been specified\n", i)
					}
				}

				printf("\tInitialising blob %d...\n", i);
				if(blob_array[i].init(i, node_filename, topology_filename, surface_filename, material_params_filename, stokes_filename, vdw_filename, pin_filename,
							scale, linear_solver, blob_state, &params, &lj_matrix, rng, num_threads) == FFEA_ERROR)
				{
					FFEA_error_text();
					printf("\tError when trying to initialise Blob %d.\n", i);
					return FFEA_ERROR;
				}

				// if centroid position is set, position the blob's centroid at that position. If vdw is set, move to center of box
				if(set_centroid_pos == 1) {
					blob_array[i].position(com_x, com_y, com_z);
					/*if(blob_state == FFEA_BLOB_IS_STATIC) {
						if(blob_array[i].create_viewer_node_file(node_filename, scale) == FFEA_ERROR) {
							FFEA_error_text();
							printf("Viewer node file could not be written. Visualisation will contain errors.\n");
							return FFEA_ERROR;
						}
					}*/
				} /*else if (set_vdw_filename == 1) {
					printf("Translating Blob %d to center of simulation box\n", i);
					box_dim.x = params.es_h * (1.0/params.kappa) * params.es_N_x;
					box_dim.y = params.es_h * (1.0/params.kappa) * params.es_N_y;
					box_dim.z = params.es_h * (1.0/params.kappa) * params.es_N_z;
	
					blob_array[i].move(box_dim.x/2.0, box_dim.y/2.0, box_dim.z/2.0);
					if(blob_state == FFEA_BLOB_IS_STATIC) {
						if(blob_array[i].create_viewer_node_file(node_filename, scale) == FFEA_ERROR) {
							FFEA_error_text();
							printf("Viewer node file could not be written. Visualisation will contain errors.\n");
							return FFEA_ERROR;
						}
					}
				}*/
				if(set_velocity == 1)
					blob_array[i].velocity_all(vel_x, vel_y, vel_z);

				// set the current node positions as pos_0 for this blob, so that all rmsd values
				// are calculated relative to this conformation centred at this point in space.
				blob_array[i].set_rmsd_pos_0();
			}

			// Make sure the next tag is the </system> tag, and not more <blob> tags
			if(get_next_script_tag(in, buf) == FFEA_ERROR) {
				FFEA_error_text();
				printf("\tError reading tag in/at end of <system> block\n");
				return FFEA_ERROR;
			}

			if(strcmp(buf, "blob") == 0) {
				FFEA_error_text();
				printf("\tToo many blob blocks. Expected only %d\n", params.num_blobs);
				return FFEA_ERROR;
			}

			if(strcmp(buf, "/system") == 0) {
				printf("All blobs fully initialised.\n");
				printf("Exiting <system> block.\n");
				return FFEA_OK;
			} else {
				FFEA_error_text();
				printf("Unrecognised tag '<%s>' when expecting '</system>' tag.\n", buf);
				return FFEA_ERROR;
			}
		}

		/* */
		void get_system_CoM(vector3 *system_CoM)
		{
			system_CoM->x = 0;
			system_CoM->y = 0;
			system_CoM->z = 0;
			scalar total_mass = 0;
			for(int i = 0; i < params.num_blobs; i++) {
				vector3 com;
				blob_array[i].get_CoM(&com);
				system_CoM->x += com.x * blob_array[i].get_mass();
				system_CoM->y += com.y * blob_array[i].get_mass();
				system_CoM->z += com.z * blob_array[i].get_mass();

				total_mass += blob_array[i].get_mass();
			}
			system_CoM->x /= total_mass;
			system_CoM->y /= total_mass;
			system_CoM->z /= total_mass;
		}

	private:
		/* How many Blobs populate this world */
		int num_blobs;

		/* Array of Blob objects */
		Blob *blob_array;

		/* How many threads are available for parallelisation */
		int num_threads;

		/* An array of pointers to random number generators (for use in parallel) */
		MTRand *rng;

		/* Parameters being used for this simulation */
		SimulationParams params;

		/*
		 * Data structure keeping track of which `cell' each face lies in (where the world has been discretised into a grid of cells of dimension 1.5 kappa)
		 * so that the BEM matrices may be constructed quickly and sparsely.
		 */
		NearestNeighbourLinkedListCube lookup;

		/*
		 * Output trajectory file
		 */
		FILE *trajectory_out;

		/*
		 * Output measurement file
		 */
		FILE **measurement_out;

		/*
		 * Output stress measurement file
		 */
		FILE *stress_out;
		/*
		 *
		 */
//		SurfaceElementLookup surface_element_lookup;

		/*
		 * BEM solver for the exterior electrostatics
		 */
		BEM_Poisson_Boltzmann PB_solver;

		/*
		 * Number of surface faces in entire system
		 */
		int total_num_surface_faces;

		/*
		 * Vector of the electrostatic potential on each surface in entire system
		 */
		scalar *phi_Gamma;
		scalar *J_Gamma;
		scalar *work_vec;

		/*
		 * Biconjugate gradient stabilised solver for nonsymmetric matrices
		 */
		BiCGSTAB_solver nonsymmetric_solver;

		/*
		 * Van der Waals solver
		 */
		VdW_solver vdw_solver;

		/*
		 * LJ parameters matrix
		 */
		LJ_matrix lj_matrix;

		vector3 box_dim;

		long long step_initial;

		int get_next_script_tag(FILE *in, char *buf)
		{
			if(fscanf(in, "%*[^<]") != 0) {
				printf("White space removal error in get_next_script_tag(). Something odd has happened: this error should never occur...\n");
			}
			if(fscanf(in, "<%255[^>]>", buf) != 1) {
				FFEA_error_text();
				printf("Error reading tag in script file.\n");
				return FFEA_ERROR;
			}

			return FFEA_OK;
		}

		void apply_dense_matrix(scalar *y, scalar *M, scalar *x, int N)
		{
			int i, j;
			for(i = 0; i < N; i++) {
				y[i] = 0;
				for(j = 0; j < N; j++)
					y[i] += M[i * N + j] * x[j];
			}
		}

		void do_es()
		{
			printf("Building BEM matrices\n");
			PB_solver.build_BEM_matrices();


		//	PB_solver.print_matrices();

			// Build the poisson matrices for each blob
			printf("Building Poisson matrices\n");
			for(int i = 0; i < params.num_blobs; i++) {
				blob_array[i].build_poisson_matrices();
			}

			printf("Solving\n");
			for(int dual = 0; dual < 30; dual++) {


				// Perform Poisson solve step (K phi + E J_Gamma = rho)
				// Obtain the resulting J_Gamma
				int master_index = 0;
				for(int i = 0; i < params.num_blobs; i++) {
					blob_array[i].solve_poisson(&phi_Gamma[master_index], &J_Gamma[master_index]);
					master_index += blob_array[i].get_num_faces();
				}


				// Apply -D matrix to J_Gamma vector (work_vec = -D J_gamma)
				PB_solver.get_D()->apply(J_Gamma, work_vec);
				for(int i = 0; i < total_num_surface_faces; i++) {
					work_vec[i] *= -1;
				}

				// Solve for C matrix (C phi_Gamma = work_vec)
				nonsymmetric_solver.solve(PB_solver.get_C(), phi_Gamma, work_vec);

				scalar sum = 0.0, sumj = 0.0;
				for(int i = 0; i < total_num_surface_faces; i++) {
					sum += phi_Gamma[i];
					sumj += J_Gamma[i];
				}
				printf("\n");

				printf("<yophi> = %e <J> = %e\n", sum/total_num_surface_faces, sumj/total_num_surface_faces);

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

		void print_trajectory_and_measurement_files(int step, double wtime)
		{
			if((step - 1)%(params.check * 10) != 0) {
				printf("step = %d\n", step);
			} else {
				printf("step = %d (simulation time = %.2es, wall clock time = %.2e hrs)\n", step, (scalar)step * params.dt, (omp_get_wtime() - wtime)/3600.0);
			}

			vector3 system_CoM;
			get_system_CoM(&system_CoM);

			// Write traj and meas data
			if(measurement_out[params.num_blobs] != NULL) {
				fprintf(measurement_out[params.num_blobs], "%d\t", step);
			}
			if(stress_out != NULL) {
				fprintf(stress_out, "step\t%d\n", step);
			}

			for(int i = 0; i < params.num_blobs; i++) {
				// Write the node data for this blob
				fprintf(trajectory_out, "Blob %d, step %d\n", i, step);
				blob_array[i].write_nodes_to_file(trajectory_out);

				// Write the measurement data for this blob
				blob_array[i].make_measurements(measurement_out[i], step, &system_CoM);

				// Output stress information
				blob_array[i].make_stress_measurements(stress_out, i);

				// Output interblob_vdw info
				for(int j = i + 1; j < params.num_blobs; ++j) { 
					blob_array[i].calculate_vdw_bb_interaction_with_another_blob(measurement_out[params.num_blobs], j);
				}
			}
	
			fprintf(measurement_out[params.num_blobs], "\n");

			// Mark completed end of step with an asterisk (so that the restart code will know if this is a fully written step or if it was cut off half way through due to interrupt)
			fprintf(trajectory_out, "*\n");

			// Force all output in buffers to be written to the output files now
			fflush(trajectory_out);
			for(int i = 0; i < params.num_blobs + 1; ++i) {
				fflush(measurement_out[i]);
			}
		}

		void print_static_trajectory(int step, double wtime, int blob_index)
		{
			printf("Printing single trajectory of Blob %d for viewer\n", blob_index);
			// Write the node data for this blob
			fprintf(trajectory_out, "Blob %d, step %d\n", blob_index, step);
			blob_array[blob_index].write_nodes_to_file(trajectory_out);
		}
};

#endif
