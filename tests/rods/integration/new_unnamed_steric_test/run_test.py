import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
from ffea_from_template import not_found


def main():

    params = OmegaConf.load("params.yml")

    not_found(params["in_dir"])
    # not_found(params["out_dir"])
    if Path(params["out_dir"]).exists():
        shutil.rmtree(params["out_dir"])
    Path(params["out_dir"]).mkdir()
    Path(f"{params['out_dir']:s}/delete").mkdir()

    ffea_files = sorted(glob.glob(f"{params['in_dir']:s}/*.ffea"))
    count = 0

    # run create_straight_rod.py, unless there's a ModuleNotFoundError, in which
    # case just use the file provided

    for i, path in enumerate(ffea_files, 1):

        name = path.split("/")[-1].split(".")[0]

        print(f"{path:s}\t({i:d}/{len(ffea_files):d})")

        result = subprocess.run(
            ["ffea", f"{name:s}.ffea"],
            capture_output=True,
            text=True,
            cwd=params["in_dir"],
        )
        with open(f"{params['out_dir']:s}/{name:s}.stdout", "w") as f:
            f.write(result.stdout)
        with open(f"{params['out_dir']:s}/{name:s}.stderr", "w") as f:
            f.write(result.stderr)

        # analyse results and determine if configuration has passed here
        passed = False

        if passed:
            count += 1

    print(f"Configurations passed: {count:d}")

    if count == len(ffea_files):
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
