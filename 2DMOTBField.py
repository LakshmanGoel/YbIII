import numpy as np
import magpylib as mag
import matplotlib.pyplot as plt

# Z-axis: Axis of the push beam
# The magnets exist in the x-z plane
# Magnet 1: (-5.3 cm, 0 cm, -4 cm)
# Magnet 2: (-5.3 cm, 0 cm, 4 cm)
# Magnet 3: (5.3 cm, 0 cm, -4 cm)
# Magnet 4: (5.3 cm, 0 cm, 4 cm)

# Strength of the dipole = Magnetization of the material (B_r) * Volume of the magnet (V)
# Strength of the dipole = 3.2e5 G*cm^3
# Desired B-field gradient = 65 G/cm

# V = 2cm * 2cm * 5.5cm = 22 cm^3
# B_r ~ 14.2 - 14.8 kG (Using a kind of neodymium magnet)
# B_r * V ~ 3.124e5 - 3.256e5 G*cm^3

# At the lower limit of B_r = 14.2 kG: dBx/dx = 65.31 G/cm, dBy/dy = -65.31 G/cm, dBz/dz = 0 G/cm
# At the upper limit of B_r = 14.8 kG: dBx/dx = 68.07 G/cm, dBy/dy = -68.07 G/cm, dBz/dz = 0 G/cm

Br = 1.42  # Tesla
mu0 = 4*np.pi*1e-7
M = Br / mu0

# Create a Matplotlib figure
fig, ax = plt.subplots()

# Create an observer grid in the xy-symmetry plane
ts = np.linspace(-0.01, 0.01, 100)
grid = np.array([[(x, y, 0) for x in ts] for y in ts])
X, Y, Z = np.moveaxis(grid, -1, 0)
ax.set_xlim(-0.01, 0.01)
ax.set_ylim(-0.01, 0.01)
ax.set_aspect('equal')

ticks_cm = np.arange(-1, 1.1, 0.5)   # -1, -0.5, 0, 0.5, 1 cm

ax.set_xticks(ticks_cm * 1e-2)
ax.set_yticks(ticks_cm * 1e-2)

ax.set_xticklabels(ticks_cm)
ax.set_yticklabels(ticks_cm)

magnet_1 = mag.magnet.Cuboid(
    position = (-5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, -M, 0))
magnet_2 = mag.magnet.Cuboid(
    position = (-5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, -M, 0))
magnet_3 = mag.magnet.Cuboid(
    position = (5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, M, 0))
magnet_4 = mag.magnet.Cuboid(
    position = (5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 2*0.01, 5.5*0.01), magnetization = (0, M, 0))

magnets = mag.Collection(magnet_1, magnet_2, magnet_3, magnet_4)

# Compute the B-field of the magnets on the grid
B = magnets.getB(grid)
print(B.shape)
Bx, By, Bz = np.moveaxis(B, -1, 0)
# Field magnitude in G
# Field magnitude in Gauss
Bnorm_G = np.linalg.norm(B, axis=2) * 1e4  # Tesla → Gauss

# --- Colormap (magnitude) ---
pcm = ax.pcolormesh(
    X, Y, Bnorm_G,
    shading="auto",
    cmap="inferno",
)

cb = fig.colorbar(pcm, ax=ax)
cb.set_label("|B| (G)")

# --- Streamlines (direction) ---
ax.streamplot(
    X, Y, Bx, By,
    color="white",
    density=1.2,
    linewidth=0.6,
)

# --- Axes styling ---
ax.set(
    xlabel="x (cm)",
    ylabel="y (cm)",
)
ax.set_aspect("equal")

plt.tight_layout()
plt.show()

#mag.show(magnets)