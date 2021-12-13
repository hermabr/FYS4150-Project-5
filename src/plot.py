import re
import argparse
import matplotlib
import tikzplotlib
import numpy as np
import pyarma as pa
<<<<<<< HEAD
import os

import matplotlib
=======
import seaborn as sns
from enum import Enum
>>>>>>> 55543811878e0ea8111923b9bfe27a9e9895db82
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

sns.set_theme()


class PlotType(Enum):
    REAL = 1
    IMAGINARY = 2
    PROBABILITY = 3


class Plotter:
    # config variables
    fontsize = 12
    t_min = 0
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1

    def __init__(self, filename):
        self.filename = filename
        self.dt = float(filename[filename.index("_dt") + 4 : filename.index(".bin")])

<<<<<<< HEAD
        U = pa.cx_mat()
        U.load(self.filename, pa.arma_binary)
        U = np.array(U)
        U_w_h = int(np.sqrt(U.shape[1]))
        U.resize(U.shape[0], U_w_h, U_w_h)
        U = U.transpose((0, 2, 1))
        self.U = U
=======
        self.U = pa.cx_mat()
        self.U.load(self.filename, pa.arma_binary)
        self.U = np.array(self.U)
        U_w_h = int(np.sqrt(self.U.shape[1]))
        self.U.resize(self.U.shape[0], U_w_h, U_w_h)
        self.U = self.U.transpose((0, 2, 1))
>>>>>>> 55543811878e0ea8111923b9bfe27a9e9895db82

        self.probabilities = np.real(self.U) ** 2 + np.imag(self.U) ** 2

    def get_norm(self, frame_idx):
        return matplotlib.cm.colors.Normalize(
            vmin=0.0, vmax=np.max(self.probabilities[frame_idx])
        )

    def current_time(self, frame_idx):
        return self.t_min + frame_idx * self.dt

    def make_frame_plot(
        self,
        frame_idx=0,
        plot_type: PlotType = PlotType.PROBABILITY,
        time_label=True,
    ):
        if isinstance(frame_idx, float):
            if frame_idx != int(frame_idx):
                raise ValueError("frame_idx must be an integer")
            frame_idx = int(frame_idx)

        plt.close()
        self.fig = plt.figure()
        self.ax = plt.gca()

        if plot_type == PlotType.PROBABILITY:
            self.img = self.ax.imshow(
                self.probabilities[frame_idx],
                cmap=plt.get_cmap("viridis"),
                extent=[self.x_min, self.x_max, self.y_min, self.y_max],
                norm=self.get_norm(frame_idx),
            )
            cbar_label = "p(x,y,z)"
        elif plot_type == PlotType.REAL:
            self.img = self.ax.imshow(
                np.real(self.U[frame_idx]),
                cmap=plt.get_cmap("viridis"),
                extent=[self.x_min, self.x_max, self.y_min, self.y_max],
                norm=self.get_norm(frame_idx),
            )
            cbar_label = "Re(U(x,y,z))"
        elif plot_type == PlotType.IMAGINARY:
            self.img = self.ax.imshow(
                np.imag(self.U[frame_idx]),
                cmap=plt.get_cmap("viridis"),
                extent=[self.x_min, self.x_max, self.y_min, self.y_max],
                norm=self.get_norm(frame_idx),
            )
            cbar_label = "Im(U(x,y,z))"

        if time_label:
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

        #  self.time_txt.set_text("")

        # Axis labels
        plt.xlabel("x", fontsize=self.fontsize)
        plt.ylabel("y", fontsize=self.fontsize)
        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)

        # Add a colourbar
        cbar = self.fig.colorbar(self.img, ax=self.ax)
        cbar.set_label(cbar_label, fontsize=self.fontsize)
        cbar.ax.tick_params(labelsize=self.fontsize)

        # Add a text element showing the time

    def animation(self, frame_idx):
        self.img.set_norm(self.get_norm(frame_idx))
        self.img.set_data(self.probabilities[frame_idx])
        self.time_txt.set_text(f"t = {self.current_time(frame_idx):.3e}")
        return self.img

    def make_time_plots(self):
        for t in [0, 0.001, 0.002]:
            self.make_frame_plot(
                t / self.dt, plot_type=PlotType.PROBABILITY, time_label=False
            )
            plt.title("AAA")  # TODO: Title
            self.save_tikz(
                f"{self.filename[:-4].replace('/data/', '/plots/')}_t{t:.3e}.tex"
            )
            exit()
            #  plt.show()

        self.make_frame_plot(0.002 / self.dt, PlotType.REAL)
        plt.title("AAA")
        plt.show()

        self.make_frame_plot(0.002 / self.dt, PlotType.IMAGINARY)
        plt.title("AAA")
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

<<<<<<< HEAD
    def detect(self):
        U_w_h = self.U.shape[1]
        x = int(U_w_h * 0.8)
        detection = self.probabilities[-1, :, x].copy()
        detection /= np.sum(detection)
        y = np.linspace(0, 1, len(detection))
        plt.ylabel(f"$p(y | x=0.8; t={(len(self.U) - 1) * self.dt})$")
        plt.xlabel("x")
        plt.plot(y, detection)
        plt.show()
=======
    def tweak_tikz_plots(self, filename):
        """Tweaks the tikz plots to make them look better

        Parameters
        ----------
            filename : str
                The filename of the tikz plot to be tweaked
        """
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            for line in lines:
                if ".png" in line:
                    line = re.sub(
                        r"{(.*\.png)}",
                        r"{Plots/\1}",
                        line,
                    )
                    f.write(line)
                elif "majorticks" in line:
                    f.write(line.replace("false", "true"))
                else:
                    f.write(line)

    def save_tikz(self, filename):
        """Saves the plot as a tikz-tex file

        Parameters
        ----------
            filename : str
                The filename of the tikz plot to be saved
        """
        #  plt.grid(True)
        #  tikzplotlib.clean_figure()
        tikzplotlib.save(filename)
        self.tweak_tikz_plots(filename)
        plt.close()


def detect(filename):
    U = pa.cx_mat()
    U.load(filename, pa.arma_binary)
    U = np.array(U)
    U_w_h = int(np.sqrt(U.shape[1]))
    U.resize(U.shape[0], U_w_h, U_w_h)
    U = U.transpose((0, 2, 1))
    x = int(U_w_h * 0.8)
    probabilities = np.real(U[-1, :, x]) ** 2 + np.imag(U[-1, :, x]) ** 2
    probabilities /= np.sum(probabilities)
    y = np.linspace(0, 1, len(probabilities))
    plt.plot(y, probabilities)
    plt.show()
>>>>>>> 55543811878e0ea8111923b9bfe27a9e9895db82


if __name__ == "__main__":
    # TODO: Do we want to run the c++ code from python?
    parser = argparse.ArgumentParser(description="To run the python plotting")

<<<<<<< HEAD

=======
    parser.add_argument(
        "-f",
        "--filename",
        type=str,
        help="The filename of the binary file to plot",
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
>>>>>>> 55543811878e0ea8111923b9bfe27a9e9895db82
    parser.add_argument(
        "-d",
        "--detect",
        help="Detect and plot at x=0.8 at the end of the simulation. Include filename",
        action="store_true",
    )

    parser.add_argument(
        "-f",
        "--filename",
        help="filename of input",
        type=str,
        default="output/data/*"
    )

    parser.add_argument(
        "-a",
        "--all",
        help="To run everything",
        action="store_true"
    )

    args = parser.parse_args()

<<<<<<< HEAD
    filenames = [args.filename[:-1] + filename for filename in os.listdir(args.filename[:-1])] if args.filename[-1] == "*" else [args.filename]

    if args.all:
        for filename in filenames:
            print(filename)
            plot = Plotter(filename)
            plot.detect()
            plot.make_time_plots()
            plot.create_animation(show=True, save=False)

=======
    #  if not any(vars(args).values()):
    #      parser.print_help()
    if args.all or True:
        plot = Plotter("output/data/double_slit_dt_0.000025.bin")
        plot.make_time_plots()
        plot.create_animation(show=True, save=False)
    if args.detect or args.all:
        detect("output/data/UBER.bin")
>>>>>>> 55543811878e0ea8111923b9bfe27a9e9895db82
