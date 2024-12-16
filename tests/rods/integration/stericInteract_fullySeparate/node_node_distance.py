# Compute and write node to node distances between two rods in a N^2 loop
import sys
import argparse
import numpy as np
import matplotlib as mpl

mpl.use("Agg")  # Allows plotting without an X-server
import matplotlib.pyplot as plt

import ffeatools.ffea_rod as ffea_rod

def node_node_distance(node_pos_a, node_pos_b):
    """Arguments are [N, 3] numpy arrays. Return N^2 node-node distance array."""

    if node_pos_a.shape[1] != 3:
        raise Exception("Unexpected node position shape:" + str(node_pos_a.shape))

    result = np.zeros((node_pos_a.shape[0], node_pos_b.shape[0]))
    for i, ra in enumerate(node_pos_a):
        for j, rb in enumerate(node_pos_b):
            result[i, j] = np.linalg.norm(rb - ra)

    return result


def plot(args, dist, frame, name, rod1, rod2):
    # Discrete colorbar
    cmap = plt.cm.RdBu
    bounds = np.array([0, args.radius * 1e9, 0.99 * args.radius * 2 * 1e9, dist.max() * 1e9])
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    plt.imshow(dist * 1e9, cmap=cmap, norm=norm, interpolation="nearest")
    plt.xlabel("Node index (rod 1)")
    plt.ylabel("Node index (rod 2)")
    plt.title("Node-node distance (nm), frame " + str(frame))
    plt.colorbar()

    plt.savefig(
        "{0:s}/{1:s}_distanceHeatmap.png".format(args.out_dir, name),
        dpi=300,
    )
    plt.close()

    frames = np.arange(rod1.num_frames, dtype='int') + 1
    disp = [
        np.linalg.norm(rod1.current_r[:, :] - rod1.current_r[0, :], axis=2) * 1e9,
        np.linalg.norm(rod2.current_r[:, :] - rod2.current_r[0, :], axis=2) * 1e9,
    ]

    for rod_no, dp in enumerate(disp, 1):
        plt.title("Rod {0:d}".format(rod_no))
        plt.xlabel("Frame index")
        plt.ylabel("|Node displacement| (nm)")
        plt.xlim(frames[0], frames[-1])
        for node_ind in range(dp.shape[1]):
            plt.plot(frames, dp[:, node_ind], "-", label="Node {0:d}".format(node_ind))
        plt.legend()
        plt.savefig(
            "{0:s}/{1:s}_absDisplacementNodes_rod{2:d}.png".format(
                args.out_dir, name, rod_no
            ),
            dpi=300,
        )
        plt.close()


def main():

    parser = argparse.ArgumentParser(description="Paths to two .rodtraj files.")
    parser.add_argument("-t1", "--rod_1_traj", type=str, required=True)
    parser.add_argument("-t2", "--rod_2_traj", type=str, required=True)
    parser.add_argument("-r", "--radius", type=float, required=True)
    parser.add_argument("-l", "--length", type=float, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    args = parser.parse_args()

    rod1 = ffea_rod.ffea_rod(filename=args.rod_1_traj)
    rod2 = ffea_rod.ffea_rod(filename=args.rod_2_traj)

    frame = rod1.num_frames
    dist = node_node_distance(rod1.current_r[frame - 1], rod2.current_r[frame - 1])
    name = args.rod_1_traj.split("/")[-1].split(".")[0][:-2]

    np.savetxt(
        "{0:s}/{1:s}_distanceHeatmap.txt".format(args.out_dir, name),
        dist,
    )

    plot(args, dist, frame, name, rod1, rod2)


if __name__ == "__main__":
    main()
    sys.exit()
