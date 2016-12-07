// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

#include "SimulationParams.h"

SimulationParams::SimulationParams() {

    // Constructor used to give default values where necessary

    // Defaulted (into SI units)
    dt = 1e-14 / mesoDimensions::time;
    num_steps = 1e11;
    check = 10000;
    kT = 4.11e-21 / mesoDimensions::Energy;
    max_iterations_cg = 1000;
    epsilon2 = 0.01;
    rng_seed = time(NULL);
    calc_vdw = 1;
    inc_self_vdw = 1;
    sticky_wall_xz = 0;
    vdw_type = "steric";
    vdw_steric_factor = 1e-2;
    move_into_box = 1;
    es_update = 1;
    kappa = 1e9 * mesoDimensions::length;
    es_h = 3;

    calc_noise = 1;
    calc_es = 0;
    dielec_ext = 1;
    epsilon_0 = 1;
    calc_stokes = 1;
    stokes_visc = 1e-3 / (mesoDimensions::pressure * mesoDimensions::time);
    calc_kinetics = 0;
    kinetics_update = 0;
    calc_preComp = 0;
    force_pbc = 0;
    calc_springs = 0; 
    calc_ctforces = 0;


    wall_x_1 = WALL_TYPE_PBC;
    wall_x_2 = WALL_TYPE_PBC;
    wall_y_1 = WALL_TYPE_PBC;
    wall_y_2 = WALL_TYPE_PBC;
    wall_z_1 = WALL_TYPE_PBC;
    wall_z_2 = WALL_TYPE_PBC;

    // Initialised to zero or equivalent for later initialisation
    restart = -1;
    num_blobs = 0;
    num_conformations = NULL;
    num_states = NULL;
    state_array_size = 0;
    conformation_array_size = 0;
    es_N_x = -1;
    es_N_y = -1;
    es_N_z = -1;

    trajectory_out_fname_set = 0;
    kinetics_out_fname_set = 0;
    measurement_out_fname_set = 0;
    icheckpoint_fname_set = 0;
    ocheckpoint_fname_set = 0;
    bsite_in_fname_set = 0;
    vdw_in_fname_set = 0;

    trajectory_out_fname = "\n";
    kinetics_out_fname = "\n";
    measurement_out_fname = "\n";
    vdw_in_fname = "\n";
    bsite_in_fname = "\n";
    icheckpoint_fname = "\n";
    ocheckpoint_fname = "\n";
    detailed_meas_out_fname = "\n";
    ctforces_fname = "\n";
}

SimulationParams::~SimulationParams() {
    dt = 0;
    num_steps = -1;
    check = 0;
    num_blobs = 0;
    delete[] num_conformations;
    num_conformations = NULL;
    delete[] num_states;
    num_states = NULL;
    state_array_size = 0;
    conformation_array_size = 0;
    rng_seed = 0;
    kT = 0;
    max_iterations_cg = 0;
    epsilon2 = 0;
    es_update = 0;
    kinetics_update = 0;
    es_N_x = -1;
    es_N_y = -1;
    es_N_z = -1;
    move_into_box = 0;
    es_h = 0;
    kappa = 0;
    dielec_ext = 0;
    epsilon_0 = 0;
    restart = 0;
    calc_vdw = -1;
    inc_self_vdw = 0;
    vdw_type = "";
    calc_es = 0;
    calc_noise = 0;
    calc_preComp = 0;
    calc_springs = 0; 
    calc_ctforces = 0;
    calc_kinetics = 0;

    wall_x_1 = -1;
    wall_x_2 = -1;
    wall_y_1 = -1;
    wall_y_2 = -1;
    wall_z_1 = -1;
    wall_z_2 = -1;

    sticky_wall_xz = 0;

    calc_stokes = 0;
    stokes_visc = -1;

    vdw_steric_factor = 0;

    trajectory_out_fname_set = 0;
    kinetics_out_fname_set = 0;
    measurement_out_fname_set = 0;
    icheckpoint_fname_set = 0;
    ocheckpoint_fname_set = 0;
    bsite_in_fname_set = 0;
    vdw_in_fname_set = 0;

    trajectory_out_fname = "\n";
    measurement_out_fname = "\n";
    kinetics_out_fname, "\n";
    icheckpoint_fname = "\n";
    ocheckpoint_fname = "\n";
    bsite_in_fname = "\n";
    vdw_in_fname = "\n";
    detailed_meas_out_fname = "\n";
    ctforces_fname = "\n";
}

int SimulationParams::extract_params(vector<string> script_vector) {

	// Extract param string from script string
	vector<string> param_vector;
	FFEA_input_reader *paramreader = new FFEA_input_reader();
	if ( paramreader->extract_block("param", 0, script_vector, &param_vector) == FFEA_ERROR) {
		return FFEA_ERROR;
	}

	// Parse the section
	vector<string>::iterator it;
	string lrvalue[2];
	for(it = param_vector.begin(); it != param_vector.end(); ++it) {
		if (paramreader->parse_tag(*it, lrvalue) == FFEA_ERROR) return FFEA_ERROR;

		// Assign if possible
		if(assign(lrvalue[0], lrvalue[1]) != 0) {
			FFEA_error_text();
			cout << "Assignment of parameter " << lrvalue[0] << " = " << lrvalue[1] << " failed" << endl;
			return FFEA_ERROR;
		}
	}
	return FFEA_OK;
}

int SimulationParams::assign(string lvalue, string rvalue) {

	b_fs::path ffea_script = FFEA_script_filename;
	FFEA_script_path = ffea_script.parent_path();
	FFEA_script_basename = ffea_script.stem();

	// Carry out parameter assignments
	if (lvalue == "restart") {
		restart = atoi(rvalue.c_str());
		cout << "\tSetting " << lvalue << " = " << restart << endl;

	} else if (lvalue == "dt") {
		dt = atof(rvalue.c_str());
		cout << "\tSetting " << lvalue << " = " << dt << endl;
                dt /= mesoDimensions::time;

	} else if (lvalue == "epsilon") {
		epsilon2 = atof(rvalue.c_str());
		epsilon2 *= epsilon2;
		cout << "\tSetting " << lvalue << " = " << epsilon2 << endl;

	} else if (lvalue == "num_steps") {

		// convert to float first so that user may write numbers of the form 1e4 for 10000 etc
	        num_steps = (long long) (atof(rvalue.c_str()));
        	cout << "\tSetting " << lvalue << " = " << num_steps << endl;

   	} else if (lvalue == "max_iterations_cg") {
		max_iterations_cg = atoi(rvalue.c_str());
		cout << "\tSetting " << lvalue << " = " << max_iterations_cg << endl;

    	} else if (lvalue == "check") {
       		check = (int) atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << check << endl;

	} else if (lvalue == "num_blobs") {
        	num_blobs = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << num_blobs << endl;

    	} else if (lvalue == "num_conformations") {

		// Split rvalue into a vector
		vector<string> conformation_vector;
		rvalue = boost::erase_last_copy(boost::erase_first_copy(rvalue, "("), ")");
		boost::split(conformation_vector, rvalue, boost::is_any_of(","));

		// Now assign them to an array
		vector<string>::iterator it;
		conformation_array_size = conformation_vector.size();
		num_conformations = new int[conformation_array_size];
		int i = 0;
		for(it = conformation_vector.begin(); it != conformation_vector.end(); ++it) {
			num_conformations[i] = atoi((*it).c_str());
			cout << "\tSetting Blob " << i << ", " << lvalue << " = " << num_conformations[i] << endl;
			i++;
		}

        } else if (lvalue == "num_states") {

		// Split rvalue into a vector
		vector<string> state_vector;
		rvalue = boost::erase_last_copy(boost::erase_first_copy(rvalue, "("), ")");
		boost::split(state_vector, rvalue, boost::is_any_of(","));

		// Now assign them to an array
		vector<string>::iterator it;
		state_array_size = state_vector.size();
        	num_states = new int[state_array_size];
		int i = 0;
		for(it = state_vector.begin(); it != state_vector.end(); ++it) {
			num_states[i] = atoi((*it).c_str());
			cout << "\tSetting Blob " << i << ", " << lvalue << " = " << num_states[i] << endl;
			i++;
		}

    	} else if (lvalue == "es_update") {
        	es_update = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << es_update << endl;

    	} else if (lvalue == "kinetics_update") {
        	kinetics_update = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << kinetics_update << endl;

    	} else if (lvalue == "es_N_x") {
        	es_N_x = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << es_N_x << endl;

	} else if (lvalue == "es_N_y") {
        	es_N_y = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << es_N_y << endl;

	} else if (lvalue == "es_N_z") {
        	es_N_z = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << es_N_z << endl;

	} else if (lvalue == "move_into_box") {
        	move_into_box = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << move_into_box << endl;

	} else if (lvalue == "sticky_wall_xz") {
        	sticky_wall_xz = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << sticky_wall_xz << endl;

	} else if (lvalue == "es_h") {
        	es_h = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << es_h << endl;

	} else if (lvalue == "rng_seed") {
		if(rvalue == "time") {
			rng_seed = time(NULL);
		} else {
			rng_seed = atoi(rvalue.c_str());
		}
        	cout << "\tSetting " << lvalue << " = " << rng_seed << endl;

    	} else if (lvalue == "kT") {
        	kT = atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << kT << endl;
                kT /= mesoDimensions::Energy;

	} else if (lvalue == "kappa") {
        	kappa = atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << kappa << endl;
        	kappa *= mesoDimensions::length;

    	} else if (lvalue == "dielec_ext") {
        	dielec_ext = atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << dielec_ext << endl;

  	} else if (lvalue == "epsilon_0") {
        	epsilon_0 = atof(rvalue.c_str()); // relative permittivity
        	cout << "\tSetting " << lvalue << " = " << epsilon_0 << endl;

    	} else if (lvalue == "calc_vdw") {
        	calc_vdw = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_vdw << endl;

    	} else if (lvalue == "inc_self_vdw") {
        	inc_self_vdw = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << inc_self_vdw << endl;

    	} else if (lvalue == "vdw_type") {
        	vdw_type = rvalue;
        	cout << "\tSetting " << lvalue << " = " << vdw_type << endl;

	} else if (lvalue == "calc_es") {
        	calc_es = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_es << endl;

	} else if (lvalue == "calc_noise") {
        	calc_noise = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_noise << endl;

	} else if (lvalue == "calc_preComp") {
        	calc_preComp = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_preComp << endl;

	} else if (lvalue == "calc_springs") {
        	calc_springs = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_springs << endl;

	} else if (lvalue == "force_pbc") {
        	force_pbc = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << force_pbc << endl;

	} else if (lvalue == "calc_kinetics") {
        	calc_kinetics = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_kinetics << endl;

   } else if (lvalue == "calc_ctforces") {
        	calc_ctforces = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_ctforces << endl;

	} else if (lvalue == "calc_stokes" || lvalue == "do_stokes") {
        	calc_stokes = atoi(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << calc_stokes << endl;

	} else if (lvalue == "stokes_visc") {
        	stokes_visc = atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << stokes_visc << endl;
                stokes_visc /= mesoDimensions::pressure * mesoDimensions::time;

	} else if (lvalue == "vdw_steric_factor") {
        	vdw_steric_factor = atof(rvalue.c_str());
        	cout << "\tSetting " << lvalue << " = " << vdw_steric_factor << endl;

	} else if (lvalue == "wall_x_1") {
		if (rvalue == "PBC") {
			wall_x_1 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_x_1 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_x_1 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_x_1 << endl;

	} else if (lvalue == "wall_x_2") {
		if (rvalue == "PBC") {
			wall_x_2 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_x_2 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_x_2 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_x_2 << endl;

	} else if (lvalue == "wall_y_1") {
		if (rvalue == "PBC") {
			wall_y_1 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_y_1 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_y_1 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_y_1 << endl;

	} else if (lvalue == "wall_y_2") {
		if (rvalue == "PBC") {
			wall_y_2 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_y_2 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_y_2 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_y_2 << endl;

	} else if (lvalue == "wall_z_1") {
		if (rvalue == "PBC") {
			wall_z_1 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_z_1 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_z_1 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_z_1 << endl;

	} else if (lvalue == "wall_z_2") {
		if (rvalue == "PBC") {
			wall_z_2 = WALL_TYPE_PBC;
		} else if (rvalue == "HARD") {
			wall_z_2 = WALL_TYPE_HARD;
		} else if (rvalue == "STOP") {
			wall_z_2 = WALL_TYPE_STOP;
		}
		cout << "\tSetting " << lvalue << " = " << wall_z_2 << endl;

    	} else if (lvalue == "trajectory_out_fname") {
			b_fs::path auxpath = FFEA_script_path / rvalue;
			trajectory_out_fname = auxpath.string();
			trajectory_out_fname_set = 1;
			cout << "\tSetting " << lvalue << " = " << trajectory_out_fname << endl;

		} else if (lvalue == "det_measurement_out_fname") {
			b_fs::path auxpath = FFEA_script_path / rvalue;
			detailed_meas_out_fname = auxpath.string();

    	} else if (lvalue == "measurement_out_fname") {
			b_fs::path auxpath = FFEA_script_path / rvalue;
			measurement_out_fname = auxpath.string();
			measurement_out_fname_set = 1;
			cout << "\tSetting " << lvalue << " = " << measurement_out_fname << endl;

		// Break up the meas fname to get default names for optional detailed measurements.
		if (detailed_meas_out_fname == "\n") {
			string meas_basename = measurement_out_fname;
			meas_basename = RemoveFileExtension(meas_basename);
			detailed_meas_out_fname = meas_basename + ".fdm";
		}

    	} else if (lvalue == "kinetics_out_fname") {
			b_fs::path auxpath = FFEA_script_path / rvalue;
			kinetics_out_fname = auxpath.string();
			kinetics_out_fname_set = 1;
			cout << "\tSetting " << lvalue << " = " << kinetics_out_fname << endl;

	} else if (lvalue == "vdw_in_fname" || lvalue == "vdw_forcefield_params") {
		b_fs::path auxpath = FFEA_script_path / rvalue;
		vdw_in_fname = auxpath.string();
		vdw_in_fname_set = 1;
		cout << "\tSetting " << lvalue << " = " << vdw_in_fname << endl;

	} else if (lvalue == "checkpoint_in") {
		b_fs::path auxpath = FFEA_script_path / rvalue;
		icheckpoint_fname = auxpath.string();
		icheckpoint_fname_set = 1;
		cout << "\tSetting " << lvalue << " = " << icheckpoint_fname << endl;

	} else if (lvalue == "checkpoint_out") {
		b_fs::path auxpath = FFEA_script_path / rvalue;
		ocheckpoint_fname = auxpath.string();
		ocheckpoint_fname_set = 1;
		cout << "\tSetting " << lvalue << " = " << ocheckpoint_fname << endl;

	} else if (lvalue == "bsite_in_fname") {
		b_fs::path auxpath = FFEA_script_path / rvalue;
		bsite_in_fname = auxpath.string();
		bsite_in_fname_set = 1;
		cout << "\tSetting " << lvalue << " = " << bsite_in_fname << endl;

	} else if (lvalue == "stress_out_fname") {
		cout << lvalue << " no longer recognised" << endl;
	} else {
        	FFEA_error_text();
        	cout << "Error: '" << lvalue << "' is not a recognised lvalue" << endl;
        	cout << "Recognised lvalues are:" << endl;
        	cout << "\tdt\n\tepsilon\n\tnum_steps\n\tmax_iterations_cg\n\tcheck\n\tes_update\n\ttrajectory_out_fname\n\tmeasurement_out_fname\n\tstress_out_fname\n\tes_N_x\n\tes_N_y\n\tes_N_z\n\tes_h\n\trng_seed\n\tkT\n\tkappa\n\tdielec_ext\n\tepsilon_0\n\trestart\n\tcalc_vdw\n\tcalc_noise\n\tcalc_preComp\n" << endl;
        	return FFEA_ERROR;
   	}

	return FFEA_OK;
}

// rename oFile to something else.
int SimulationParams::checkFileName(string oFile){


    if (b_fs::exists(oFile))    // does oFile actually exist?
    {
        int cnt = 1;
        b_fs::path fs_oFile = oFile;
        string base = "__" + fs_oFile.filename().string() + "__bckp.";
        if (fs_oFile.parent_path().string().size() != 0) {
          b_fs::path fs_base = fs_oFile.parent_path() / base;
          base = fs_base.string();
        }

        string bckp = base + boost::lexical_cast<string>(cnt);
        while (b_fs::exists(bckp)) {
          cnt += 1;
          string s_cnt = boost::lexical_cast<string>(cnt);
          bckp = base + s_cnt;
        }
        FFEA_CAUTION_MESSG("Moving %s to %s\n", oFile.c_str(), bckp.c_str());
        // cout << "FFEA: moving " << oFile << " to " << bckp << "\n";
        b_fs::rename(oFile, bckp.c_str());
    }


    return FFEA_OK;
}

int SimulationParams::validate() {

    if (restart != 0 && restart != 1) {
        FFEA_ERROR_MESSG("Required: Restart flag, 'restart', must be 0 (false) or 1 (true).\n");
    }
    if (restart == 1) {
        if (icheckpoint_fname_set == 0) {
            FFEA_ERROR_MESSG("Checkpoint input file required if restart is required.\n");
        }
    }

    if (num_steps < 0) {
        FFEA_ERROR_MESSG("Required: Number of time steps, 'num_steps', must be greater than or equal to 0.\n");
    }
    if (num_blobs <= 0) {
        FFEA_ERROR_MESSG("\tRequired: Number of Blobs, 'num_blobs', must be greater than 0.\n");
    }

    if (conformation_array_size != num_blobs) {
	FFEA_ERROR_MESSG("\tRequired: Number of Conformations, 'num_conformations', must have 'num_blobs' elements. We read %d elements but only %d blobs\n", conformation_array_size, num_blobs);
    }

    if (kappa < 0) {
        FFEA_ERROR_MESSG("Required: Inverse Debye Screening length, 'kappa', cannot be negative.\n");
    }

    if (dielec_ext <= 0) {
        FFEA_ERROR_MESSG("Required: Exterior dielectric constant, 'dielec_ext', must be greater than 0.\n");
    }

    if (epsilon_0 <= 0) {
        FFEA_ERROR_MESSG("Required: Permittivity of free space, 'epsilon_0', must be greater than 0.\n");
    }

    if (calc_vdw != 0 && calc_vdw != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_vdw', must be 0 (no) or 1 (yes).\n");
    }

    if (inc_self_vdw != 0 && inc_self_vdw != 1) {
        FFEA_ERROR_MESSG("Required: 'inc_self_vdw', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_vdw == 1) {

      if (vdw_type != "lennard-jones" && vdw_type != "steric" && vdw_type != "ljsteric") {
          FFEA_ERROR_MESSG("Optional: 'vdw_type', must be either 'steric' (default), 'lennard-jones' or 'ljsteric' (both methods combined).\n");
      }

      if (vdw_type != "steric") {
        if (vdw_in_fname_set == 0) {
            FFEA_ERROR_MESSG("VdW forcefield params file name required (vdw_forcefield_params).\n");
        }
      }
    } else {
      if (inc_self_vdw == 1) {
        printf("\tFRIENDLY WARNING: No face-face interactions will be computed if calc_vdw = 0.\n");
      }
    }

    if (calc_preComp != 0 && calc_preComp != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_preComp', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_springs != 0 && calc_springs != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_springs', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_ctforces != 0 && calc_ctforces != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_ctforces', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_es != 0 && calc_es != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_es', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_kinetics != 0 && calc_kinetics != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_kinetics', must be 0 (no) or 1 (yes).\n");
    }

    if (calc_vdw == 1 or calc_es == 1) {
    	if (es_N_x < 1) {
        	printf("\tFRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_x', is less than 1. Will assign default value to encompass whole system.\n");
    	} else if (es_N_y < 1) {
        	printf("\tFRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_y', is less than 1. Will assign default value to encompass whole system.\n");
    	} else if (es_N_z < 1) {
        	printf("\tFRIENDLY WARNING: Length of the nearest neighbour lookup grid, 'es_N_z', is less than 1. Will assign default value to encompass whole system.\n");
    	} else {
        	if (es_h <= 0) {
            		FFEA_ERROR_MESSG("Required: Nearest neighbour lookup grid cell dimension, 'es_h', must be greater than 0.\n");
        	}
    	}

	if (move_into_box != 0 && move_into_box != 1) {
		FFEA_ERROR_MESSG("'move_into_box' must all be either 0 (system centroid preserved) or 1 (system centroid will move to centroid of simulation box).\n")
	}

    } else {
		printf("\tFRIENDLY WARNING: No electrostatic or vdw interactions will be simulated\n");
		es_N_x = 0;
		es_N_y = 0;
		es_N_z = 0;
    }

    if (calc_noise != 0 && calc_noise != 1) {
        FFEA_ERROR_MESSG("Required: 'calc_noise', must be 0 (no) or 1 (yes).\n");
    }

    if (trajectory_out_fname_set == 0) {
        b_fs::path auxpath = FFEA_script_path / FFEA_script_basename / ".ftj";
        trajectory_out_fname = auxpath.string();
    }

    if (measurement_out_fname_set == 0) {
		b_fs::path auxpath = FFEA_script_path / FFEA_script_basename / ".fm";
		measurement_out_fname = auxpath.string();
    }

    // Three checkings for checkpoint files:
    // CPT.1 - If we don't have a name for checkpoint_out we're assigning one.
    if (ocheckpoint_fname_set == 0) {
      b_fs::path fs_ocpt_fname = FFEA_script_filename;
      fs_ocpt_fname.replace_extension(".fcp");
      ocheckpoint_fname_set = 1;
      ocheckpoint_fname = fs_ocpt_fname.string();
      printf("\tFRIENDLY WARNING: Checkpoint output file name was not specified, so it will be set to %s\n", ocheckpoint_fname.c_str() );
    }
    // CPT.2 - checkpoint_out must differ from checkpoint_in
    if (ocheckpoint_fname.compare(icheckpoint_fname) == 0) {
        FFEA_ERROR_MESSG("it is not allowed to set up checkpoint_in and checkpoint_out with the same file names\n");
    }
    // CPT.3 - checkpoint_out will be backed up if it exists
    checkFileName(ocheckpoint_fname);


    // check if the output files exists, and if so, rename it.
    if (restart == 0) {
      checkFileName(measurement_out_fname);
      checkFileName(detailed_meas_out_fname);
      checkFileName(trajectory_out_fname);
      checkFileName(kinetics_out_fname);
    }

    if (calc_stokes == 1 && stokes_visc <= 0) {
        FFEA_ERROR_MESSG("calc_stokes flag is set, so stokes_visc must be set to a value greater than 0.\n");
    }

    if (vdw_steric_factor < 0) {
       printf("\tFRIENDLY WARNING: Beware, vdw_steric_factor is negative.\n");
    }

    if (calc_kinetics == 1) {
	if(kinetics_update <= 0) {
		//FFEA_ERROR_MESSG("\tRequired: If 'calc_kinetics' = 1, then 'kinetics_update' must be greater than 0.\n");
		cout << "\tDefaulting 'kinetics_update' to " << check << endl;
		kinetics_update = check;
	} //else if (kinetics_update <= check) {

		// This could be fixed by changing how the output files are printed. Fix later, busy now!
		//FFEA_CAUTION_MESSG("\t'kinetics_update' < 'check'. A kinetic switch therefore maybe missed i.e. not printed to the output files.\n")
	//}

	// num_conformations[i] can be > num_states[i], so long as none of the states reference an out of bounds conformation

    } else {
	for(int i = 0; i < num_blobs; ++i) {
		if(num_conformations[i] != 1) {
			FFEA_CAUTION_MESSG("\tNumber of Conformations, 'num_conformations[%d]', not equal to 1. Only first conformation will be loaded.\n", i);
			num_conformations[i] = 1;
		}
	}
    }
    printf("...done\n");

    printf("Parameters:\n");
    printf("\trestart = %d\n", restart);
    printf("\tdt = %e\n", dt*mesoDimensions::time);
    printf("\tnum_steps = %lld\n", num_steps);
    printf("\tcheck = %d\n", check);
    printf("\tnum_blobs = %d\n", num_blobs);
    for (int i = 0; i < num_blobs; ++i) {
        printf("\tBlob %d:\n", i);
        printf("\t\tnum_conformations = %d\n", num_conformations[i]);
	if (calc_kinetics == 1) {
	        printf("\t\tnum_states = %d\n", num_states[i]);
    	}
    }
    printf("\trng_seed = %d\n", rng_seed);
    printf("\tkT = %e\n", kT*mesoDimensions::Energy);
    printf("\ttrajectory_out_fname = %s\n", trajectory_out_fname.c_str());
    printf("\tmeasurement_out_fname = %s\n", measurement_out_fname.c_str());
    printf("\tkinetics_out_fname = %s\n", kinetics_out_fname.c_str());
    if (icheckpoint_fname_set == 1) {
      printf("\tcheckpoint_in = %s\n", icheckpoint_fname.c_str());
    }
    printf("\tcheckpoint_out = %s\n", ocheckpoint_fname.c_str());
    printf("\tvdw_forcefield_params = %s\n", vdw_in_fname.c_str());
    printf("\tmax_iterations_cg = %d\n", max_iterations_cg);
    printf("\tepsilon2 = %e\n", epsilon2);
    printf("\tes_update = %d\n", es_update);
    printf("\tes_N_x = %d\n", es_N_x);
    printf("\tes_N_y = %d\n", es_N_y);
    printf("\tes_N_z = %d\n", es_N_z);
    printf("\tes_h = %e x inverse kappa\n", es_h);
    printf("\tkappa = %e\n", kappa/mesoDimensions::length);
    printf("\tepsilon_0 = %e\n", epsilon_0);
    printf("\tdielec_ext = %e\n", dielec_ext);
    printf("\tcalc_vdw = %d\n", calc_vdw);
    printf("\tvdw_type = %s\n", vdw_type.c_str());
    printf("\tinc_self_vdw = %d\n", inc_self_vdw);
    printf("\tcalc_es = %d\n", calc_es);
    printf("\tcalc_noise = %d\n", calc_noise);
    printf("\tcalc_kinetics = %d\n", calc_kinetics);
    printf("\tcalc_preComp = %d\n", calc_preComp);
    printf("\tcalc_springs = %d\n", calc_springs);
    printf("\tforce_pbc = %d\n", force_pbc);
    printf("\tcalc_stokes = %d\n", calc_stokes);
    printf("\tstokes_visc = %f\n", stokes_visc*mesoDimensions::pressure*mesoDimensions::time);
    printf("\tcalc_kinetics = %d\n", calc_kinetics);

    if(calc_kinetics == 1 && bsite_in_fname_set == 1) {
	printf("\tbsite_in_fname = %s\n", bsite_in_fname.c_str());
    }
    if(calc_vdw == 1 && (vdw_type == "steric" || vdw_type == "ljsteric")) {
        printf("\tvdw_steric_factor = %e\n", vdw_steric_factor);
    }
    return FFEA_OK;
}

int SimulationParams::get_max_num_states() {

    int i, max_num_states = 0;
    for(i = 0; i < num_blobs; ++i) {
	if(num_states[i] > max_num_states) {
	    max_num_states = num_states[i];
	}
    }
    return max_num_states;
}

void SimulationParams::write_to_file(FILE *fout) {

	// This should be getting added to the top of the measurement file!!

	// Then, every parameter (if anyone hates this goddamn class in the future, please make a Param struct / dictionary / vector / list thing so we can just loop over all parameters)
	fprintf(fout, "Parameters:\n");
    	fprintf(fout, "\trestart = %d\n", restart);
    	fprintf(fout, "\tdt = %e\n", dt*mesoDimensions::time);
    	fprintf(fout, "\tnum_steps = %lld\n", num_steps);
    	fprintf(fout, "\tcheck = %d\n", check);
	fprintf(fout, "\trng_seed = %d\n", rng_seed);
	fprintf(fout, "\tkT = %e\n", kT*mesoDimensions::Energy);

    	fprintf(fout, "\tcalc_vdw = %d\n", calc_vdw);
    	fprintf(fout, "\tvdw_type = %s\n", vdw_type.c_str());
      fprintf(fout, "\tinc_self_vdw = %d\n", inc_self_vdw);
    	fprintf(fout, "\tcalc_es = %d\n", calc_es);
    	fprintf(fout, "\tcalc_noise = %d\n", calc_noise);
    	fprintf(fout, "\tcalc_kinetics = %d\n", calc_kinetics);
    	fprintf(fout, "\tcalc_preComp = %d\n", calc_preComp);
      fprintf(fout, "\tcalc_springs = %d\n", calc_springs);
      fprintf(fout, "\tforce_pbc = %d\n", force_pbc);
    	fprintf(fout, "\tcalc_stokes = %d\n", calc_stokes);
    	fprintf(fout, "\tstokes_visc = %e\n", stokes_visc*mesoDimensions::pressure*mesoDimensions::time);

    	fprintf(fout, "\tes_update = %d\n", es_update);
    	fprintf(fout, "\tes_N_x = %d\n", es_N_x);
   	fprintf(fout, "\tes_N_y = %d\n", es_N_y);
    	fprintf(fout, "\tes_N_z = %d\n", es_N_z);
    	fprintf(fout, "\tes_h = %f\n", es_h);
	fprintf(fout, "\tmove_into_box = %d\n", move_into_box);

    	fprintf(fout, "\tkappa = %e\n", kappa/mesoDimensions::length);
    	fprintf(fout, "\tepsilon_0 = %e\n", epsilon_0);
    	fprintf(fout, "\tdielec_ext = %e\n", dielec_ext);

	fprintf(fout, "\tmax_iterations_cg = %d\n", max_iterations_cg);
    	fprintf(fout, "\tepsilon2 = %e\n", epsilon2);
	fprintf(fout, "\tnum_blobs = %d\n", num_blobs);
	fprintf(fout, "\tnum_conformations = (");
    	for (int i = 0; i < num_blobs; ++i) {
        	fprintf(fout, "%d", num_conformations[i]);
		if(i == num_blobs - 1) {
			fprintf(fout, ")");
		} else {
			fprintf(fout, ",");
		}
    	}
	fprintf(fout, "\n");

	if (calc_kinetics == 1) {
		fprintf(fout, "num_states = (");
		for (int i = 0; i < num_blobs; ++i) {
			fprintf(fout, "%d", num_states[i]);
			if(i == num_blobs - 1) {
				fprintf(fout, ")");
			} else {
				fprintf(fout, ",");
			}
		}
		fprintf(fout, "\n");
	}

    	fprintf(fout, "\n\n");
}

string SimulationParams::RemoveFileExtension(const string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(0, FileName.find_last_of("."));
    return "";
}
