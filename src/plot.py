import re
import argparse
import matplotlib
from matplotlib import animation
import tikzplotlib
import numpy as np
import pyarma as pa
import os

import matplotlib
import seaborn as sns
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

sns.set_theme()


class PlotType(Enum):
    """Plot type enum."""

    REAL = 1
    IMAGINARY = 2
    PROBABILITY = 3


class Plotter:
    """
    Class for handling plotting of the results from a system
    """

    # Config variables
    fontsize = 12
    t_min = 0
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1

    def __init__(self, filename):
        """Parses the armadillo matrix and make it ready for plotting

        Parameters
        ----------
            filename : str
                The filename of the binary representation of the matrix to be plotted. The file must contain the substring "_dt={dt as used in the system}"
        """
        self.filename = "output/data/" + filename
        self.dt = float(filename[filename.index("_dt") + 4 : filename.index(".bin")])

        self.U = pa.cx_mat()
        self.U.load(self.filename, pa.arma_binary)
        self.U = np.array(self.U)
        U_w_h = int(np.sqrt(self.U.shape[1]))
        self.U.resize(self.U.shape[0], U_w_h, U_w_h)
        self.U = self.U.transpose((0, 2, 1))
        # Add in boundary points
        U_boundary = np.zeros((self.U.shape[0], U_w_h + 2, U_w_h + 2), dtype=complex)
        U_boundary[:,1:-1,1:-1] = self.U
        self.U = U_boundary

        self.probabilities = np.real(self.U) ** 2 + np.imag(self.U) ** 2

    def get_norm(self, frame_idx, data=None):
        """Get a norm for the colormap of the probability distribution

        Parameters
        ----------
            frame_idx : int
                The index of the frame to be plotted
            data : ndarray
                The data to be plotted. If not specified, the probabilities are used

        Returns
        -------
            norm : matplotlib.colors.Normalize
                The norm for the colormap
        """
        if data is None:
            data = self.probabilities[frame_idx]
        return matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data))

    def current_time(self, frame_idx):
        """Return the current time of the system

        Parameters
        ----------
            frame_idx : int
                The index of the frame to be plotted

        Returns
        -------
            current_time : float
                The current time of the system
        """
        return self.t_min + frame_idx * self.dt

    def make_frame_plot(
        self,
        frame_idx=0,
        plot_type: PlotType = PlotType.PROBABILITY,
        is_animation=True,
    ):
        """Make a plot of the system at a given frame

        Parameters
        ----------
            frame_idx : int
                The index of the frame to be plotted
            plot_type : PlotType
                The type of plot to be made
            is_animation : bool
                Whether to make an animation or not

        Raises
        ------
            ValueError :
                If the frame is not a valid index
        """
        if isinstance(frame_idx, float):
            if frame_idx != int(frame_idx):
                raise ValueError("frame_idx must be an integer")
            frame_idx = int(frame_idx)

        plt.close()
        self.fig, self.ax = plt.subplots()

        if plot_type == PlotType.PROBABILITY:
            data = self.probabilities[frame_idx]
            cbar_label = r"$p(x, y; t)$"
        elif plot_type == PlotType.REAL:
            data = np.real(self.U[frame_idx])
            cbar_label = r"Re$(u(x, y; t))$"
        elif plot_type == PlotType.IMAGINARY:
            data = np.imag(self.U[frame_idx])
            cbar_label = r"Im$(u(x, y; t))$"

        self.img = plt.imshow(
            data,
            cmap="viridis",
            norm=self.get_norm(frame_idx, data),
            extent=[self.x_min, self.x_max, self.y_min, self.y_max],
        )

        if is_animation:
            self.time_txt = plt.text(
                0.95 * self.x_max,
                0.95 * self.y_max,
                f"t = {self.current_time(frame_idx):.3e}",
                color="white",
                horizontalalignment="right",
                verticalalignment="top",
                fontsize=self.fontsize,
            )

        # Axis labels
        plt.xlabel("x", fontsize=self.fontsize)
        plt.ylabel("y", fontsize=self.fontsize)
        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)

        # Add a colourbar
        cbar = self.fig.colorbar(self.img, ax=self.ax)
        cbar.set_label(cbar_label, fontsize=self.fontsize)
        cbar.ax.tick_params(labelsize=self.fontsize)

    def animation(self, frame_idx):
        """Make one step of the animation for the system

        Parameters
        ----------
            frame_idx : int
                The index of the frame to be plotted

        Returns
        -------
            self.img : matplotlib.image.AxesImage
                The image of the plot
        """
        self.img.set_norm(self.get_norm(frame_idx))
        self.img.set_data(self.probabilities[frame_idx])
        self.time_txt.set_text(f"t = {self.current_time(frame_idx):.3e}")
        return self.img

    def latex_plot_name(self, filename, plot_type, t):
        """Create a latex object with the code for plotting the image

        Parameters
        ----------
            filename : str
                The filename of the binary representation of the matrix to be plotted. The file must contain the substring "_dt={dt as used in the system}"
            plot_type : PlotType
                The type of plot to be made
            t : float
                The time of the plot

        Returns
        -------
            title : str
                The title of the plot
        """
        text = """\\begin{figure}\n    \\centering\n    """
        caption = f"A plot of the "
        title = ""
        if plot_type == PlotType.PROBABILITY:
            title = f"Probability distribution for t={t}"
            caption += f"probability distribution of p(x,y;t)"
        elif plot_type == PlotType.REAL:
            title = f"Real part of the wavefunction for t={t}"
            caption += f"real part of the wavefunction u(x,y;t)"
        elif plot_type == PlotType.IMAGINARY:
            title = f"Imaginary part of the wavefunction for t={t}"
            caption += f"imaginary part of the wavefunction u(x,y;t)"
        caption += rf" for t={t} using setup 1"
        text += f"\\input{{{'Plots/' + filename.split('/')[-1]}}}\n    "
        text += f"\\caption{{{caption}}} "
        text += f"\\label{{fig:{str(plot_type)}_{t}}}\n\\end{{figure}}\n"
        print(text)

        return title

    def make_time_plots(self):
        """Produce time plots for the times t=0, t=0.01 and t=0.002"""
        for t in [0, 0.001, 0.002]:
            for plot_type in [PlotType.PROBABILITY, PlotType.REAL, PlotType.IMAGINARY]:
                self.make_frame_plot(
                    t / self.dt, plot_type=plot_type, is_animation=False
                )

                filename = (
                    "output/plots/"
                    + self.filename.split("/")[-1][:-4]
                    + f"_t{t:.3e}_{plot_type}.tex"
                )

                title = self.latex_plot_name(filename, plot_type, t)
                plt.title(title)
                self.save_tikz(filename, heat_plot=True)

    def create_animation(self, show=False, save=True):
        """Create an animation of the system

        Parameters
        ----------
            show : bool
                Whether to show the animation or not. Default is False
            save : bool
                Whether to save the animation or not. Default is True
        """
        plt.grid(False)

        self.make_frame_plot()

        anim = FuncAnimation(
            self.fig,
            self.animation,
            interval=1,
            frames=np.arange(0, len(self.probabilities), 2),
            repeat=False,
            blit=0,
        )

        if show:
            plt.show()

        if save:
            anim.save(
                "output/animations/" + self.filename.split("/")[-1][:-4] + ".mp4",
                writer="ffmpeg",
                bitrate=-1,
                fps=30,
                dpi=300,
            )

    def detect(self, t=0.002):
        """Find the 1d probavility distributions given x=0.8

        Parameters
        ----------
            t : float
                The time at which to find the distribution. Default is 0.002
        """
        U_w_h = self.U.shape[1]
        x = int(U_w_h * 0.8)
        detection = self.probabilities[int(t / self.dt), :, x].copy()
        detection /= np.sum(detection)
        y = np.linspace(0, 1, len(detection))
        plt.ylabel(f"$p(y | x=0.8; t={t})$")
        plt.xlabel("x")
        plt.plot(y, detection)
        plt.title("1-d probability distribution double slit")
        self.save_tikz("output/plots/heyutherocksteadycrew.tex", line_plot=True)
        plt.show()

    def deviation_plot(self):
        """Plot the deviation of the wavefunction from the ground state"""
        total_probabilities = np.sum(self.probabilities, axis=(1, 2))
        deviations = np.abs(1 - total_probabilities)
        t = np.linspace(self.t_min, self.dt * (len(deviations) - 1), len(deviations))

        #  t = np.mean(t.reshape(-1, 5), axis=1)
        #  deviations = np.mean(deviations.reshape(-1, 5), axis=1)
        #  plt.scatter(t, deviations)

        #  plt.plot(t, deviations)
        plt.scatter(t, deviations)
        #  sns.scatterplot(t, deviations)
        self.save_tikz(
            "output/plots/" + self.filename.split("/")[-1][:-4] + "_deviations.tex",
            scatter_plot=True,
        )
        #  plt.show()

    def tweak_tikz_plots(self, filename, heat_plot, scatter_plot, line_plot):
        """Tweaks the tikz plots to make them look better by modifying some of the parameters of the tikz plots

        Parameters
        ----------
            filename : str
                The filename of the tikz plot to be tweaked
        """
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            should_write = True
            for line in lines:
                if not should_write and "end{tikzpicture}" not in line:
                    continue
                if ".png" in line:
                    line = re.sub(
                        r"{(.*\.png)}",
                        r"{Plots/\1}",
                        line,
                    )
                    f.write(line)
                elif "majorticks" in line:
                    f.write(line.replace("false", "true"))
                elif "begin{axis}[" in line:
                    f.write(line)
                    if heat_plot:
                        f.write("width=7cm,\n")
                    elif line_plot:
                        f.write("scaled y ticks=false,\n")
                    elif scatter_plot:
                        f.write("mark size=1.25,")
                elif "end{axis}" in line:
                    f.write(line)
                    should_write = False
                elif "title=" in line:
                    f.write(line.replace("{", "\\textbf{"))
                else:
                    f.write(line)

    def save_tikz(self, filename, heat_plot=False, scatter_plot=False, line_plot=False):
        """Saves the plot as a tikz-tex file

        Parameters
        ----------
            filename : str
                The filename of the tikz plot to be saved
            heat_plot : bool
                Whether to save the heat plot or not. Default is True
            scatter_plot : bool
                Whether to save the scatter plot or not. Default is False
        """
        if not sum([heat_plot, scatter_plot, line_plot]) == 1:
            raise ValueError(
                "The plot must be of one of the types heat, scatter or line"
            )
        tikzplotlib.save(filename)
        self.tweak_tikz_plots(filename, heat_plot, scatter_plot, line_plot)
        plt.close()


if __name__ == "__main__":
    # TODO: Do we want to run the c++ code from python?
    parser = argparse.ArgumentParser(description="For running the plotting")

    # TODO: RUNNING ARGUMENTS DO BE LOOKIN KINDA SUS
    parser.add_argument(
        "-f",
        "--filename",
        type=str,
        help="The filename of the binary file to plot",
        default="output/data/*",
    )

    parser.add_argument(
        "-t",
        "--time_plots",
        help="To plot the probablity and the real and imaginary parts of the wave function for different time steps",
        action="store_true",
    )
    parser.add_argument(
        "-an",
        "--animation",
        help="To create an animation of probablity",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--detect",
        help="Detect and plot at x=0.8 at the end of the simulation. Include filename",
        action="store_true",
    )
    parser.add_argument(
        "-de",
        "--deviation",
        help="Plot deviation from 1 of the total probability",
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--all",
        help="Produce all plots and animation",
        action="store_true",
    )
    args = parser.parse_args()

    filenames = (
        [args.filename[:-1] + filename for filename in os.listdir(args.filename[:-1])]
        if args.filename[-1] == "*"
        else [args.filename]
    )

    if (
        not args.time_plots
        and not args.animation
        and not args.detect
        and not args.all
        and not args.deviation
    ):
        parser.print_help()
        exit()

    for filename in filenames:
        plot = Plotter(filename)
        if args.time_plots or args.all:
            print(f"Plotting time plots for {filename}")
            plot.make_time_plots()
        if args.animation or args.all:
            print(f"Creating animation for {filename}")
            plot.create_animation()
        if args.detect or args.all:
            print(f"Plotting detect for {filename}")
            plot.detect()
        if args.deviation or args.all:
            print(f"Plotting error for {filename}")
            plot.deviation_plot()
