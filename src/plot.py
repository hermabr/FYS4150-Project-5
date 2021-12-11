import numpy as np
import pyarma as pa

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#  U = pa.cx_cube()
U = pa.cx_mat()  # format (t, x, y)

U.load("UBER.bin", pa.arma_binary)

#TODO: add boundary points (all with value 0)
U = np.array(U)
U = U.reshape(320, 199, 199)  # format (t, y, x)

probabilities = np.real(U) ** 2 + np.imag(U) ** 2

#
# Now the list z_data_list contains a series of "frames" of z(x,y,t),
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = 0
x_min, x_max = 0, 5
y_min, y_max = 0, 5

# Create figure
fig = plt.figure()
ax = plt.gca()

# Create a colour scale normalization according to the max z value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(probabilities))

# Plot the first frame
img = ax.imshow(
    probabilities[0],
    extent=[x_min, x_max, y_min, y_max],
    cmap=plt.get_cmap("viridis"),
    norm=norm,
)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label("z(x,y,t)", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(
    0.95,
    0.95,
    f"t = {t_min:.3e}",
    color="white",
    horizontalalignment="right",
    verticalalignment="top",
    fontsize=fontsize,
)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(probabilities[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(probabilities[i])
    dt = 1
    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img


# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(
    fig,
    animation,
    interval=1,
    frames=np.arange(0, len(probabilities), 2),
    repeat=False,
    blit=0,
)

# Run the animation!
plt.show()

# # Save the animation
anim.save("./animation.mp4", writer="ffmpeg", bitrate=-1, fps=30)
