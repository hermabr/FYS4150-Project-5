import argparse
import numpy as np
import pyarma as pa

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Plotter:
    # config variables
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

    def get_norm(self, frame_idx):
        return matplotlib.cm.colors.Normalize(
            vmin=0.0, vmax=np.max(self.probabilities[frame_idx])
        )

    def current_time(self, frame_idx):
        return self.t_min + frame_idx * self.dt

    def make_frame_plot(self, frame_idx=0):
        if isinstance(frame_idx, float):
            if frame_idx != int(frame_idx):
                raise ValueError("frame_idx must be an integer")
            frame_idx = int(frame_idx)

        plt.close()
        self.fig = plt.figure()
        self.ax = plt.gca()

        # Plot the first frame
        self.img = self.ax.imshow(
            self.probabilities[frame_idx],
            extent=[self.x_min, self.x_max, self.y_min, self.y_max],
            cmap=plt.get_cmap("viridis"),
            norm=self.get_norm(frame_idx),
        )

        # Axis labels
        plt.xlabel("x", fontsize=self.fontsize)
        plt.ylabel("y", fontsize=self.fontsize)
        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)

        # Add a colourbar
        cbar = self.fig.colorbar(self.img, ax=self.ax)
        cbar.set_label("p(x,y,t)", fontsize=self.fontsize)
        cbar.ax.tick_params(labelsize=self.fontsize)

        # Add a text element showing the time
        self.time_txt = plt.text(
            0.95 * self.x_max,
            0.95 * self.y_max,
            #  f"t = {self.t_min:.3e}",
            f"t = {self.current_time(frame_idx):.3e}",
            color="white",
            horizontalalignment="right",
            verticalalignment="top",
            fontsize=self.fontsize,
        )

    def animation(self, frame_idx):
        self.img.set_norm(self.get_norm(frame_idx))
        self.img.set_data(self.probabilities[frame_idx])
        self.time_txt.set_text(f"t = {self.current_time(frame_idx):.3e}")
        return self.img

    def make_time_plots(self):
        for t in [0, 0.001, 0.002]:
            self.make_frame_plot(t / self.dt)
            plt.title("AAA")
            #  plt.savefig(
            #      self.filename.replace("/data/", "/plots/").replace(
            #          ".bin", f"_{t:.3e}.png"
            #      )
            #  )
            plt.show()

    def create_animation(self, show=False, save=True):
        self.make_frame_plot()

        # Use matplotlib.animation.FuncAnimation to put it all together
        anim = FuncAnimation(
            self.fig,
            self.animation,
            interval=1,
            frames=np.arange(0, len(self.probabilities), 2),
            repeat=False,
            blit=0,
        )

        if show:
            # Run the animation!
            plt.show()

        if save:
            # Save the animation as an mp4. This requires ffmpeg or mencoder to be
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
        plot.make_time_plots()
        #  plot.create_animation(show=True, save=False)
