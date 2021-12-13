import re
import argparse
import matplotlib
import tikzplotlib
import numpy as np
import pyarma as pa
import seaborn as sns
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

        self.U = pa.cx_mat()
        self.U.load(self.filename, pa.arma_binary)
        self.U = np.array(self.U)
        U_w_h = int(np.sqrt(self.U.shape[1]))
        self.U.resize(self.U.shape[0], U_w_h, U_w_h)
        self.U = self.U.transpose((0, 2, 1))

        self.probabilities = np.real(self.U) ** 2 + np.imag(self.U) ** 2

    def get_norm(self, frame_idx, data=None):
        if data is None:
            data = self.probabilities[frame_idx]
        return matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(data))

    def current_time(self, frame_idx):
        return self.t_min + frame_idx * self.dt

    def make_frame_plot(
        self,
        frame_idx=0,
        plot_type: PlotType = PlotType.PROBABILITY,
        is_animation=True,
    ):
        if isinstance(frame_idx, float):
            if frame_idx != int(frame_idx):
                raise ValueError("frame_idx must be an integer")
            frame_idx = int(frame_idx)

        plt.close()
        #  self.fig = plt.figure()
        #  self.ax = plt.gca()
        self.fig, self.ax = plt.subplots()
        #  self.ax.grid(False)
        #  self.ax.grid(False)

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
        self.img.set_norm(self.get_norm(frame_idx))
        self.img.set_data(self.probabilities[frame_idx])
        self.time_txt.set_text(f"t = {self.current_time(frame_idx):.3e}")
        return self.img

    def make_time_plots(self):
        for t in [0, 0.001, 0.002]:
            for plot_type in [PlotType.REAL, PlotType.IMAGINARY, PlotType.PROBABILITY]:
                self.make_frame_plot(
                    t / self.dt, plot_type=plot_type, is_animation=False
                )

                filename = f"{self.filename[:-4].replace('/data/', '/plots/')}_t{t:.3e}_{plot_type}.tex"

                text = """\\begin{figure}\n    \\centering\n    """
                title = ""
                if plot_type == PlotType.PROBABILITY:
                    title = f"Probability distribution for p(x,y,t={t})"
                elif plot_type == PlotType.REAL:
                    title = f"Real part of the wavefunction u(x,y,t={t})"
                elif plot_type == PlotType.IMAGINARY:
                    title = f"Imaginary part of the wavefunction u(x,y,t={t})"
                text += f"\\text{{{title}}}\\par\\medskip\n    "
                text += (
                    f"\\input{{{filename.replace('output/plots/', 'Plots/')}}}\n    "
                )
                text += "\\caption{An example plot} \\label{fig:example_plot}\n\\end{figure}\n"
                print(text)

                self.save_tikz(filename)

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
                #  elif "begin{tikzpicture}" in line:
                #      f.write(line[:-1] + "[scale=0.80]" + "\n")
                elif "begin{axis}[" in line:
                    f.write(line)
                    f.write("width=7cm,\n")
                elif "end{axis}" in line:
                    f.write(line)
                    should_write = False
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


if __name__ == "__main__":
    # TODO: Do we want to run the c++ code from python?
    parser = argparse.ArgumentParser(description="To run the python plotting")

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
    parser.add_argument(
        "-d",
        "--detect",
        help="Detect and plot at x=0.8 at the end of the simulation",
        action="store_true",
    )

    parser.add_argument(
        "-a",
        "--all",
        help="To run everything",
        action="store_true",
    )
    args = parser.parse_args()

    #  if not any(vars(args).values()):
    #      parser.print_help()
    if args.all or True:
        plot = Plotter("output/data/double_slit_dt_0.000025.bin")
        plot.make_time_plots()
        #  plot.create_animation(show=True, save=False)
    if args.detect or args.all:
        detect("output/data/UBER.bin")
