# Generate rod configurations and create FFEA input files from them
from pathlib import Path
import shutil
from omegaconf import OmegaConf
import sys
import math


def not_found(file_path: str):
    if not Path(file_path).exists():
        raise FileNotFoundError(f"'{file_path:s}' not found.")


def copy_template(filename: str):
    template = Path("template.ffea")
    if not template.exists():
        raise Exception(f"{template} not found")
    shutil.copy(template.name, filename)
    if not Path(filename).exists():
        raise Exception(f"Could not copy {template.name} to {filename}")
    return Path(filename)


def write_test_config_yaml(R: str, L: str, config_filename: str):
    """
    Write translations and rotations of a rod to a .yml file and return as a dict.

    x: left (-ve) / right (+ve)
    y: up/down (relative to rod axis)
    z: rod axis

    rotation: anticlockwise (-ve) / clockwise (+ve) about an axis, with the origin
    at the rod midpoint
    """

    a_small = math.pi * 1 / 36
    a_big = math.pi * 1 / 3
    a_90 = math.pi * 1 / 2
    half_L_R = 0.5 * (L + R)
    half_L = 0.5 * L

    # avoid using underscores or periods in dict keys
    config = {
        "parallel": {
            "+x": [R, 0, 0, 0, 0, 0],
            "-x": [-R, 0, 0, 0, 0, 0],
            "+y": [0, R, 0, 0, 0, 0],
            "-y": [0, -R, 0, 0, 0, 0],
            "+z": [0, 0, L + R, 0, 0, 0],
            "-z": [0, 0, -L - R, 0, 0, 0],
            "+xzHalf": [R, 0, half_L_R, 0, 0, 0],
            "-xzHalf": [-R, 0, -half_L_R, 0, 0, 0],
            "+xyz": [R, R, L + R, 0, 0, 0],
            "-xyz": [-R, -R, -L - R, 0, 0, 0],
        },
        "perp": {
            # crossed
            "Cross+x+Rx": [R, 0, 0, a_90, 0, 0],
            "Cross+x-Rx": [R, 0, 0, -a_90, 0, 0],
            "Cross-x+Rx": [-R, 0, 0, a_90, 0, 0],
            "Cross-x-Rx": [-R, 0, 0, -a_90, 0, 0],
            "Cross+y+Ry": [0, R, 0, 0, a_90, 0],
            "Cross+y-Ry": [0, R, 0, 0, -a_90, 0],
            "Cross-y+Ry": [0, -R, 0, 0, a_90, 0],
            "Cross-y-Ry": [0, -R, 0, 0, -a_90, 0],
            # T-shape
            "T+z+Rx": [0, 0, half_L_R, a_90, 0, 0],
            "T-z+Rx": [0, 0, -half_L_R, a_90, 0, 0],
            "T+z+Ry": [0, 0, half_L_R, 0, a_90, 0],
            "T-z+Ry": [0, 0, -half_L_R, 0, a_90, 0],
            # L-shape
            "L+x+z+Ry": [half_L, 0, half_L_R, 0, a_90, 0],
            "L-x+z+Ry": [-half_L, 0, half_L_R, 0, a_90, 0],
            "L+x-z+Ry": [half_L, 0, -half_L_R, 0, a_90, 0],
            "L-x-z+Ry": [-half_L, 0, -half_L_R, 0, a_90, 0],
            "L+y+z+Rx": [0, half_L, half_L_R, a_90, 0, 0],
            "L-y+z+Rx": [0, -half_L, half_L_R, a_90, 0, 0],
            "L+y-z+Rx": [0, half_L, -half_L_R, a_90, 0, 0],
            "L-y-z+Rx": [0, -half_L, -half_L_R, a_90, 0, 0],
        },
        "oblique": {
            "Cross+x+RxSmall": [R, 0, 0, a_small, 0, 0],
            "Cross-x+RxSmall": [-R, 0, 0, a_small, 0, 0],
            "Cross+x-RxSmall": [R, 0, 0, -a_small, 0, 0],
            "Cross-x-RxSmall": [-R, 0, 0, -a_small, 0, 0],
            "Cross+x+RxBig": [R, 0, 0, a_big, 0, 0],
            "Cross-x+RxBig": [-R, 0, 0, a_big, 0, 0],
            "Cross+x-RxBig": [R, 0, 0, -a_big, 0, 0],
            "Cross-x-RxBig": [-R, 0, 0, -a_big, 0, 0],
            "Cross+y+RySmall": [0, R, 0, 0, a_small, 0],
            "Cross-y+RySmall": [0, -R, 0, 0, a_small, 0],
            "Cross+y-RySmall": [0, R, 0, 0, -a_small, 0],
            "Cross-y-RySmall": [0, -R, 0, 0, -a_small, 0],
            "Cross+y+RyBig": [0, R, 0, 0, a_big, 0],
            "Cross-y+RyBig": [0, -R, 0, 0, a_big, 0],
            "Cross+y-RyBig": [0, R, 0, 0, -a_big, 0],
            "Cross-y-RyBig": [0, -R, 0, 0, -a_big, 0],
            "T+z+RxSmall": [
                0,
                0,
                half_L_R,
                a_90 - a_small,
                0,
                0,
            ],
            "T+z+RxBig": [
                0,
                0,
                half_L_R,
                a_90 - a_big,
                0,
                0,
            ],
        },
        "thruFail": {
            "parallel": [0, 0, 0, 0, 0, 0],
            "parallelHalf": [0, 0, half_L, 0, 0, 0],
            "perpC": [0, 0, 0, a_90, 0, 0],
            "perpT": [0, 0, half_L, a_90, 0, 0],
            "perpL": [half_L, 0, half_L, a_90, 0, 0],
            "obliqueSmall": [0, 0, 0, a_small, 0, 0],
            "obliqueBig": [0, 0, 0, a_big, 0, 0],
        },
    }

    conf = OmegaConf.create(config)
    yaml_str = OmegaConf.to_yaml(conf)
    with open(config_filename, "w") as f:
        f.write(yaml_str)

    return config


def edit_ffea_field(file_line: str, val_old: str, val_new: str) -> str:
    """Edits a string of format "<foo = bar>" read from a .ffea file."""

    whitespace = file_line.split("<")[0]
    field = file_line.split("<")[-1].split(">")[0]  # assumes no comments in file
    field = field.replace(val_old, val_new)
    return f"{whitespace:s}<{field:s}>\n"


def edit_ffea_file(
    ffea_file_path: str,
    trans_old: str,
    trans_new: "tuple[float]",
    rot_old: str,
    rot_new: "tuple[float]",
):
    """
    Edits in-place the rod traj, centroid, and rotation fields of a .ffea file.

    The first rod occurring in the file is transformed.
    """

    if len(trans_new) and len(rot_new) != 3:
        raise Exception("Translation and rotation vectors must be length 3")

    count = 1
    name = ffea_file_path.split("/")[-1].split(".")[0]

    with open(ffea_file_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "output" and "rodtraj" in line:
            lines[i] = edit_ffea_field(
                line, f"{count:d}.rodtraj", f"{name:s}_{count:d}.rodtraj"
            )
            count += 1

        # Only apply transformation to a single rod of the pair
        if "centroid_pos" in line and count == 2:
            lines[i] = edit_ffea_field(
                line,
                trans_old,
                f"({trans_new[0]:.2f}, {trans_new[1]:.2f}, {trans_new[2]:.2f})",
            )

        if "rotation" in line and count == 2:
            lines[i] = edit_ffea_field(
                line, rot_old, f"({rot_new[0]:.2f}, {rot_new[1]:.2f}, {rot_new[2]:.2f})"
            )

    with open(ffea_file_path, "w") as f:
        f.writelines(lines)

    return lines


def main():

    not_found("params.yml")
    not_found("template.ffea")

    # Define rod configurations
    params = OmegaConf.load("params.yml")
    radius = params["radius"]
    length = params["length"]
    test_config = write_test_config_yaml(
        radius / params["ffea_length"],
        length / params["ffea_length"],
        f"{params['in_dir']:s}/test_config.yml",
    )

    not_found(params["in_dir"])
    not_found(f"{params['in_dir']:s}/test_config.yml")

    # Generate .ffea files for colliding rod pairs
    for key1, val1 in test_config.items():
        for key2, val2 in val1.items():

            path = params["in_dir"] + f"/{key1}_{key2}.ffea"
            copy_template(path)
            edit_ffea_file(
                ffea_file_path=path,
                trans_old="(0.0, 0.0, 0.0)",
                trans_new=val2[:3],
                rot_old="(0.0, 0.0, 0.0)",
                rot_new=val2[3:],
            )


if __name__ == "__main__":
    main()
    sys.exit()
