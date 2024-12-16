import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np


def final_frame_index(rodtraj_path: str):
    """Note that rod trajectories, as of 16/03/23, seem to have two extra frames
    compared to 'Nframes / check' computed from the .ffea file."""
    n = -1
    with open(rodtraj_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "FRAME" in line:
                n = int(line.split(" ")[1])
    return n


def write_failed_info(fail_info_path: str, config_names, return_codes, output_dir):
    with open(f"{fail_info_path:s}", "w") as f:
        f.write("config,last_frame,return_code\n")
        for name, code in zip(config_names, return_codes):
            frame = final_frame_index(f"{output_dir:s}/{name:s}_1.rodtraj")
            f.write(f"{name:s},{frame:d},{code:d}\n")


def copy_failed(fail_info_path: str, input_dir: str, output_dir : str):

    if Path("input_failed").exists():
        shutil.rmtree("input_failed")
    if Path("output_failed").exists():
        shutil.rmtree("output_failed")
    Path("input_failed").mkdir()
    Path("output_failed").mkdir()

    with open(fail_info_path, "r") as f:
        lines = f.readlines()[1:]
        for i, line in enumerate(lines):
            lines[i] = line.split(",")[0] + ".ffea"
            shutil.copyfile(f"{input_dir:s}/{lines[i]:s}", f"input_failed/{lines[i]:s}")

    check_copied = len(glob.glob(f"input_failed/*.ffea"))
    check_failed = len(lines)

    if check_copied != check_failed:
        raise Exception(
            f"Inconsistent copying of failed FFEA files to 'input_failed'"
        )

    ffea_files = glob.glob("input_failed/*.ffea")
    for path_in in ffea_files:
        prefix = path_in.split("/")[-1].split(".")[0]
        out_files = glob.glob(f"{output_dir}/{prefix}*")
        for path_out in out_files:
            f_out = path_out.split("/")[-1]
            shutil.copyfile(f"{path_out}", f"output_failed/{f_out}")



def main():

    params = OmegaConf.load("params.yml")

    if Path(params["in_dir"]).exists():
        shutil.rmtree(params["in_dir"])
    Path(params["in_dir"]).mkdir()

    if Path(params["out_dir"]).exists():
        shutil.rmtree(params["out_dir"])
    Path(params["out_dir"]).mkdir()
    Path(f"{params['out_dir']:s}/delete").mkdir()

    # Set up simulations
    subprocess.run([f"python", "create_straight_rod.py"])
    subprocess.run(["python", "ffea_from_template.py"])

    if params["run_failed_only"]:
        ffea_files = []
        with open(params["fail_info_path"], "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                name = line.split(",")[0]
                ffea_files.append(f"{params['in_dir']:s}/{name:s}.ffea")

        print("Running failed configurations only\n")
    else:
        ffea_files = glob.glob(f"{params['in_dir']:s}/*.ffea")

    count = 0
    num_configs = len(ffea_files)
    failed_configs = []
    failed_return_codes = []

    # ! - Tests assume steric_force_factor = 20 (rod_structure.h)

    for i, path in enumerate(sorted(ffea_files), 1):

        name = path.split("/")[-1].split(".")[0]

        print(f"{name:s}\t({i:d}/{num_configs:d})")
        print("FFEA simulation...")
        ffea_result = subprocess.run(
            ["ffea", f"{name:s}.ffea"],
            capture_output=True,
            text=True,
            cwd=params["in_dir"],
        )

        print("Writing stdout and stderr...")
        with open(f"{params['out_dir']:s}/{name:s}.stdout", "w") as f:
            f.write(ffea_result.stdout)
        with open(f"{params['out_dir']:s}/{name:s}.stderr", "w") as f:
            f.write(ffea_result.stderr)

        if "thruFail" in name and ffea_result.returncode != 0:
            count += 1
            print("Passed\n")
            continue

        print("Results analysis...")
        subprocess.run(
            [
                f"python",
                "node_node_distance.py",
                "--rod_1_traj",
                f"{params['out_dir']:s}/{name:s}_1.rodtraj",
                "--rod_2_traj",
                f"{params['out_dir']:s}/{name:s}_2.rodtraj",
                "--radius",
                str(params["radius"]),
                "--length",
                str(params["length"]),
                "--out_dir",
                params["out_dir"],
            ]
        )

        dist = np.loadtxt(f"{params['out_dir']:s}/{name:s}_distanceHeatmap.txt")
        # print("Node-node distances (nm):")
        # with np.printoptions(precision=2, suppress=False):
        #     print(dist * 1e9)
        if dist[dist <= 0.99 * (2 * params["radius"])].size == 0 and ffea_result.returncode == 0:
            count += 1
            print("Passed\n")
        else:
            print("Failed\n")
            failed_configs.append(name)
            failed_return_codes.append(ffea_result.returncode)

    print(f"Total passed: {count:d}/{num_configs:d}")

    write_failed_info(
        params["fail_info_path"], failed_configs, failed_return_codes, params["out_dir"]
    )

    if Path("FFEA_meta.json").exists():
        Path("FFEA_meta.json").unlink()

    if len(failed_configs) > 0:
        fail_dir = f"{params['in_dir']:s}_failed"
        copy_failed(
            fail_info_path=params["fail_info_path"],
            input_dir=params["in_dir"],
            output_dir=params["out_dir"]
        )

    if count == num_configs:
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
