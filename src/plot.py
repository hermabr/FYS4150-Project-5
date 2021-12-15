import re
import argparse
import matplotlib
import tikzplotlib
import numpy as np
import pyarma as pa
import os

import matplotlib
import seaborn as sns
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
                The filename of the binary representation of the matrix to be plotted. The file must contain the filename (excluding suffix) of the config used in the simulation"
        """
        self.filename = filename
        self.dt = self.get_dt_from_config()

        self.U = pa.cx_mat()
        # parse the binary file using pyarma
        self.U.load(self.filename, pa.arma_binary)
        self.U = np.array(self.U, dtype=np.clongdouble)
        U_w_h = int(np.sqrt(self.U.shape[1]))
        # reshape from 2d to 3d
        self.U.resize(self.U.shape[0], U_w_h, U_w_h)
        # swap the x and y axis, since they are indexed differently in armadillo and numpy
        self.U = self.U.transpose((0, 2, 1))

        # pad the matrix with zeros on all sides (boundary points)
        self.U = np.pad(
            self.U,
            pad_width=((0, 0), (1, 1), (1, 1)),
            mode="constant",
            constant_values=0,
        )

        self.probabilities = np.real(self.U) ** 2 + np.imag(self.U) ** 2

    def get_dt_from_config(self):
        """Get the dt from the config file

        Returns
        -------
            dt : float
                The dt from the config file

        Raises
        ------
            ValueError
                If the dt is not found in the config file
        """
        # search the filename for either simple, double or triple, and then take the next part (which is the name of the config file)
        config = (
            re.search(r"(?:simple|double|triple)_(.*)\..*", filename).group(1) + ".in"
        )
        with open(config) as f:
            for line in f:
                if line.startswith("dt"):
                    return float(line.split(" ")[1])
        raise ValueError("Could not find dt in config file")

    def get_slit_type(self):
        """Get the slit type from the config file

        Returns
        -------
            slit_type : str
                The slit type from the config file
        """
        # extract either simple, double or triple from the filename
        return re.search(r"(simple|double|triple)", filename).group(1)

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
        # use probability data if no data is specified
        if data is None:
            data = self.probabilities[frame_idx]
        # Returns the normalizer. The normalizer will have have vmin is either 0
        # or min(data). This is because we want negative values be shown for
        # the real and imaginary parts, but we never want the heat plots to
        # start at a higher value than 0
        return matplotlib.cm.colors.Normalize(
            vmin=min(0, np.min(data)), vmax=np.max(data)
        )

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
        # check that the frame is a whole number if it is a float
        if isinstance(frame_idx, float):
            if frame_idx != int(frame_idx):
                raise ValueError("frame_idx must be an integer")
            frame_idx = int(frame_idx)

        # clear old plots before making the new one
        plt.close()
        self.fig, self.ax = plt.subplots()

        # get the data and label according to the plot type
        if plot_type == PlotType.PROBABILITY:
            data = self.probabilities[frame_idx]
            cbar_label = r"$p(x, y; t)$"
        elif plot_type == PlotType.REAL:
            data = np.real(self.U[frame_idx])
            cbar_label = r"Re$(u(x, y; t))$"
        elif plot_type == PlotType.IMAGINARY:
            data = np.imag(self.U[frame_idx])
            cbar_label = r"Im$(u(x, y; t))$"

        # make a plot of the actual data to be plotted
        self.img = plt.imshow(
            data,
            cmap="viridis",
            norm=self.get_norm(frame_idx, data),
            extent=[self.x_min, self.x_max, self.y_min, self.y_max],
        )

        # add the time label if it is an animation
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
        # update the norm to fit the data of the current frame
        self.img.set_norm(self.get_norm(frame_idx))
        # update the actual data to be the currect frame
        self.img.set_data(self.probabilities[frame_idx])
        # update the time label
        self.time_txt.set_text(f"t = {self.current_time(frame_idx):.3e}")
        return self.img

    def latex_plot_name(self, filename, plot_type, t):
        """Create a latex object with the code for plotting the image

        Parameters
        ----------
            filename : str
                The filename of the binary representation of the matrix to be plotted. The file must contain the filename (excluding suffix) of the config used in the simulation"
            plot_type : PlotType
                The type of plot to be made
            t : float
                The time of the plot

        Returns
        -------
            title : str
                The title of the plot
        """
        # initialize the text, caption and title
        text = """\\begin{figure}\n    \\centering\n    """
        caption = f"A plot of the "
        title = ""
        # add a title and caption for the plot depending on the plot type
        if plot_type == PlotType.PROBABILITY:
            title = f"Probability distribution for t={t}"
            caption += f"probability distribution of p(x,y;t)"
        elif plot_type == PlotType.REAL:
            title = f"Real part of u for t={t}"
            caption += f"real part of the wavefunction u(x,y;t)"
        elif plot_type == PlotType.IMAGINARY:
            title = f"Imaginary part of u for t={t}"
            caption += f"imaginary part of the wavefunction u(x,y;t)"

        # add the lower part of the latex plots and print
        caption += rf" for t={t} using setup 3"
        text += f"\\input{{{'Plots/' + filename.split('/')[-1][:-4]}}}\n    "
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

                # the file should be saved to the folder 'output/plots', and have information for both the time and plot type
                filename = (
                    "output/plots/"
                    + self.filename.split("/")[-1][:-4]
                    + f"_t{t:.3e}_{plot_type}.tex"
                )

                # add the correct title and save the plot to file
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
        # initialize the heat plot and hide the grid
        self.make_frame_plot()
        plt.grid(False)

        # create the actual animation part
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
            # save the animation to file, using a high dpi and 30 fps
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
        # get the width of the matrix
        U_w_h = self.U.shape[1]
        # set x to be 0.8 * the width of the matrix (x=0.8)
        x = int(U_w_h * 0.8)
        # cut out the interesting part of the matrix (where x=0.8)
        detection = self.probabilities[int(t / self.dt), :, x].copy()
        # normalize the results
        detection /= np.sum(detection)
        y = np.linspace(0, 1, len(detection))
        # add labels to the plot
        plt.ylabel(f"$p(y | x=0.8; t={t})$")
        plt.xlabel("y")
        plt.plot(y, detection)
        plt.title(f"1-d probability distribution {self.get_slit_type()} slit")
        # save the plot to file
        self.save_tikz(
            "output/plots/" + self.filename.split("/")[-1][:-4] + "_detect.tex",
            line_plot=True,
        )
        plt.show()

    def deviation_plot(self):
        """Plot the deviation of the wavefunction from the ground state"""
        # get the total probability per time step by summing over all the x and y values, per t
        total_probabilities = np.sum(self.probabilities, axis=(1, 2))
        # get the deviation by subtracting the ground state probability from the total probability
        deviations = np.abs(1 - total_probabilities)
        t = np.linspace(self.t_min, self.dt * (len(deviations) - 1), len(deviations))
        # add labels to the plot
        plt.plot(t, deviations)
        plt.title("Deviation from 1 of the total probability")
        plt.xlabel("t")
        plt.ylabel(r"$|1 - \sum_{i,j} p(x_i,y_j;t)|$")
        # save the plot to file
        self.save_tikz(
            "output/plots/" + self.filename.split("/")[-1][:-4] + "_deviations.tex",
            line_plot=True,
        )

    def tweak_tikz_plots(self, filename, heat_plot, scatter_line_plot, line_plot):
        """Tweaks the tikz plots to make them look better by modifying some of the parameters of the tikz plots

        Parameters
        ----------
            filename : str
                The filename of the tikz plot to be tweaked
            heat_plot: bool
                True if heat plot
            scatter_line_plot: bool
                True if scatter and line plot
            line_plot: bool
                True if line plot
        """
        # read the entire plot file, which will be used when tweaking the plot
        with open(filename, "r") as f:
            lines = f.readlines()

        with open(filename, "w") as f:
            # is true only while writing the first axis (to avoid having duplicate plots)
            should_write = True
            for line in lines:
                if not should_write and "end{tikzpicture}" not in line:
                    continue
                if ".png" in line:
                    # add the folder of the plots if the plot contains a png (used for the heat plots)
                    line = re.sub(
                        r"{(.*\.png)}",
                        r"{Plots/\1}",
                        line,
                    )
                    f.write(line)
                elif "majorticks" in line:
                    # show the major ticks
                    f.write(line.replace("false", "true"))
                elif "begin{axis}[" in line:
                    f.write(line)
                    if heat_plot:
                        # set the width of the plot, to avoid having the heat plot be too wide
                        f.write("width=7cm,\n")
                    elif line_plot or scatter_line_plot:
                        # show the scaled y ticks on the side, not on top to avoid collisions with the title
                        f.write("scaled y ticks=false,\n")
                elif "end{axis}" in line:
                    f.write(line)
                    # stop writing of the axis was finished
                    should_write = False
                elif "title=" in line:
                    # make the title bold
                    f.write(line.replace("{", "\\textbf{"))
                else:
                    # the default case, just write the line
                    f.write(line)

    def save_tikz(
        self, filename, heat_plot=False, scatter_line_plot=False, line_plot=False
    ):
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
        # check that exactly on of the parameters are 1
        if not sum([heat_plot, scatter_line_plot, line_plot]) == 1:
            raise ValueError(
                "The plot must be of one of the types heat, scatter or line"
            )
        # actually save the plot as a tikz-tex file
        tikzplotlib.save(filename)
        # tweak the plot to make it look better
        self.tweak_tikz_plots(filename, heat_plot, scatter_line_plot, line_plot)
        plt.close()


if __name__ == "__main__":
    # define the arg parser, used for taking command line arguments
    parser = argparse.ArgumentParser(description="For running the plotting")

    # add the arguments
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

    # get all the filenames specified. if no filename is specified, use all the files in the folder
    filenames = (
        [args.filename[:-1] + filename for filename in os.listdir(args.filename[:-1])]
        if args.filename[-1] == "*"
        else [args.filename]
    )

    # check that at least one of the plots is specified. if not: write out a help message
    if (
        not args.time_plots
        and not args.animation
        and not args.detect
        and not args.all
        and not args.deviation
    ):
        parser.print_help()
        exit()

    # loop through the files and make the plots specified by the arguments
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
