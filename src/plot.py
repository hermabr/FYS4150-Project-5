import argparse
import numpy as np
import pyarma as pa

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Plotter:
    # config
    fontsize = 12
    t_min = 0
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1

    def __init__(self, filename):
        self.filename = filename
        self.dt = float(filename[filename.index("_dt") + 4 : filename.index(".bin")])

        U = pa.cx_mat()
        U.load(self.filename, pa.arma_binary)
        U = np.array(U)
        U_w_h = int(np.sqrt(U.shape[1]))
        U.resize(U.shape[0], U_w_h, U_w_h)
        U = U.transpose((0, 2, 1))

        self.probabilities = np.real(U) ** 2 + np.imag(U) ** 2

        self.fig = plt.figure()
        self.ax = plt.gca()

    def animation(self, i):
        norm = matplotlib.cm.colors.Normalize(
            vmin=0.0, vmax=np.max(self.probabilities[i])
        )
        self.img.set_norm(norm)
        self.img.set_data(self.probabilities[i])
        current_time = self.t_min + i * self.dt
        self.time_txt.set_text(f"t = {current_time:.3e}")
        return self.img

    def create_animation(self, show=False):
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(self.probabilities))

        # Plot the first frame
        self.img = self.ax.imshow(
            self.probabilities[0],
            extent=[self.x_min, self.x_max, self.y_min, self.y_max],
            cmap=plt.get_cmap("viridis"),
            norm=norm,
        )

        # Axis labels
        plt.xlabel("x", fontsize=self.fontsize)
        plt.ylabel("y", fontsize=self.fontsize)
        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)

        # Add a colourbar
        cbar = self.fig.colorbar(self.img, ax=self.ax)
        cbar.set_label("z(x,y,t)", fontsize=self.fontsize)
        cbar.ax.tick_params(labelsize=self.fontsize)

        # Add a text element showing the time
        self.time_txt = plt.text(
            0.95 * self.x_max,
            0.95 * self.y_max,
            f"t = {self.t_min:.3e}",
            color="white",
            horizontalalignment="right",
            verticalalignment="top",
            fontsize=self.fontsize,
        )

        # Use matplotlib.animation.FuncAnimation to put it all together
        anim = FuncAnimation(
            self.fig,
            self.animation,
            interval=1,
            frames=np.arange(0, len(self.probabilities), 2),
            repeat=False,
            blit=0,
        )

        # Run the animation!
        if show:
            plt.show()

        #  anim.save(self.filename + ".mp4", writer="ffmpeg")
        anim.save(
            self.filename.replace("/data/", "/animations/").replace(".bin", ".mp4"),
            writer="ffmpeg",
            bitrate=-1,
            fps=30,
            dpi=300,
        )


if __name__ == "__main__":
    # TODO: Do we want to run the c++ code from python?
    parser = argparse.ArgumentParser(
        description="To run the python plotting and the c++ simulation"
    )
    parser.add_argument(
        "-s",
        "--states",
        help="To get information about different states",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="To generate plots",
        action="store_true",
    )
    #  parser.add_argument(
    #      "-an",
    #      "--analytical",
    #      help="To generate analytical values",
    #      action="store_true",
    #  )
    #  parser.add_argument(
    #      "-z",
    #      "--zoom",
    #      help="Find maximum values and zoom",
    #      action="store_true",
    #  )
    #  parser.add_argument(
    #      "-r",
    #      "--reproduce",
    #      help="Reproduce the experiment as done in the report, using the same seed",
    #      action="store_true",
    #  )
    parser.add_argument(
        "-a",
        "--all",
        help="To run everything",
        action="store_true",
    )
    args = parser.parse_args()
    if args.plot or args.all or True:
        plot = Plotter("output/data/double_slit_dt_0.000025.bin")
        plot.create_animation()
