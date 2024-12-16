import sys
import os
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np
import pprint as pp
import matplotlib.pyplot as plt

import ffeatools.ffea_rod as ffea_rod


def plot_displacement_vs_time_all_nodes(rod, rod_no: int, show=False):

    frames = np.arange(rod.num_frames, dtype='int') + 1
    dp = np.linalg.norm(rod.current_r[:, :] - rod.current_r[0, :], axis=2) * 1e9

    plt.title("Rod {0:d}".format(rod_no))
    plt.xlabel("Frame index")
    plt.ylabel("|Node displacement| (nm)")
    plt.xlim(frames[0], frames[-1])
    for node_ind in range(dp.shape[1]):
        plt.plot(frames, dp[:, node_ind], "-", label="Node {0:d}".format(node_ind))
    plt.legend()

    if show:
        plt.show()

    plt.savefig("displacement_vs_time_nodes_rod{0:d}.png".format(rod_no), dpi=300)
    plt.close()

def plot_displacement_vs_time_all_vdw_sites(rod, rod_no: int, show=False):

    frames = np.arange(rod.num_frames, dtype='int') + 1
    dp = np.linalg.norm(rod.vdw_site_pos[:, :] - rod.vdw_site_pos[0, :], axis=2) * 1e9

    plt.title("Rod {0:d}".format(rod_no))
    plt.xlabel("Frame index")
    plt.ylabel("|Node displacement| (nm)")
    plt.xlim(frames[0], frames[-1])
    for node_ind in range(dp.shape[1]):
        plt.plot(frames, dp[:, node_ind], "-", label="Node {0:d}".format(node_ind))
    plt.legend()

    if show:
        plt.show()

    plt.savefig("displacement_vs_time_vdw_sites_rod{0:d}.png".format(rod_no), dpi=300)
    plt.close()

def plot_rod_with_vdw_sites(rod, rod_no: int, frame=0, show=False):

    fig2 = plt.figure(2)
    ax3d = fig2.add_subplot(111, projection='3d')


    ax3d.plot(rod.current_r[frame, 0, 0].flatten(), rod.current_r[frame, 0, 1].flatten(), rod.current_r[frame, 0, 2].flatten(), 'bo')
    ax3d.plot(rod.current_r[frame, 1:, 0].flatten(), rod.current_r[frame, 1:, 1].flatten(), rod.current_r[frame, 1:, 2].flatten(), 'go', label=f"Node")
    ax3d.plot(rod.vdw_site_pos[frame, :, 0].flatten(), rod.vdw_site_pos[frame, :, 1].flatten(), rod.vdw_site_pos[frame, :, 2].flatten(), 'r*', label=f"VDW Site")

    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")

    if show:
        plt.show()
    else:
        plt.savefig(f"rod_{rod_no:d}_preview_with_sites.png", dpi=300)
    plt.close()


def analysis(rod_fnames):

    for i, fname in enumerate(rod_fnames):
        print(f"\nROD {i:d}")
        print(fname)
        rod = ffea_rod.ffea_rod(filename=fname, vdw_filename=fname.split(".rod")[0] + ".rodvdw")
        plot_displacement_vs_time_all_nodes(rod, i)
        plot_displacement_vs_time_all_vdw_sites(rod, i)
        plot_rod_with_vdw_sites(rod, i, 0, True)

        # print("Site positions: ")
        # pp.pprint(np.array(site_pos))
        # print("Node positions: ")
        # pp.pprint(self.current_r[frame_index, :, :])
        # print("Site relative positions: ")
        # pp.pprint((site_pos[:, :] - self.current_r[frame_index, 0, :]) / self.get_contour_length(frame_index))
        # print("Node relative positions: ")
        # pp.pprint((self.current_r[frame_index, :, :] - self.current_r[frame_index, 0, :]) / self.get_contour_length(frame_index))


def main():

    params = OmegaConf.load("params.yml")

    os.system("./cleanup.sh")
    subprocess.run([f"python", "create_rod_with_vdw_sites.py"])

    ffea_files = glob.glob(f"*.ffea")
    count = 0
    num_configs = len(ffea_files)

    for i, path in enumerate(sorted(ffea_files), 1):

        name = path.split("/")[-1].split(".")[0]

        # print(f"{name:s}\t({i:d}/{num_configs:d})")
        # name = "ljtest"
        # print("FFEA simulation...")
        # ffea_result = subprocess.run(
        #     ["ffea", "{}}.ffea"],
        #     capture_output=True,
        #     text=True,
        # )

        os.system("ffea ljtest.ffea")

        rod = ffea_rod.ffea_rod(filename="rod1.rodtraj", vdw_filename="rod1.rodvdw")
        rod.load_trajectory()
        pp.pprint(rod.vdw_site_pos_from_rodtraj[1])
        pp.pprint(rod.current_r[1])

        # print("Writing stdout and stderr...")
        # with open(f"{name:s}.stdout", "w") as f:
        #     f.write(ffea_result.stdout)
        # with open(f"{name:s}.stderr", "w") as f:
        #     f.write(ffea_result.stderr)

        rod_fnames = glob.glob(f"{params['out_dir']}*.rodtraj")

        # analysis(rod_fnames)

        if False:
            count += 1
            print("Passed\n")
        else:
            print("Failed\n")

    print(f"Total passed: {count:d}/{num_configs:d}")

    if Path("FFEA_meta.json").exists():
        Path("FFEA_meta.json").unlink()

    if count == num_configs:
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
