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

# V = 2cm * 5.5cm * 2cm = 22 cm^3
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
ts = np.linspace(-0.01*np.sqrt(2), 0.01*np.sqrt(2), 100) # Use (-0.01*np.sqrt(2), 0.01*np.sqrt(2), 100) when rotating the axes, otherwise use (-0.01, 0.01, 100)
grid = np.array([[(x, y, 0) for x in ts] for y in ts]) # Put z = 0.01 m to plot B_z 1 cm above the z=0 plane
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
    position = (-5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
magnet_2 = mag.magnet.Cuboid(
    position = (-5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, -M, 0))
magnet_3 = mag.magnet.Cuboid(
    position = (5.3*0.01, 0, -4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))
magnet_4 = mag.magnet.Cuboid(
    position = (5.3*0.01, 0, 4*0.01), dimension = (2*0.01, 5.5*0.01, 2*0.01), magnetization = (0, M, 0))

magnets = mag.Collection(magnet_1, magnet_2, magnet_3, magnet_4)

# Compute the B-field of the magnets on the grid
B = magnets.getB(grid)
Bx, By, Bz = np.moveaxis(B, -1, 0)
# Field magnitude in G
# Field magnitude in Gauss
Bnorm_G = np.linalg.norm(B, axis=2) * 1e4  # Tesla → Gauss

#-------------------------------------For the colormap plots, uncomment from here-------------------------------------#
# Bx_p = (Bx - By) / np.sqrt(2)
# By_p = (Bx + By) / np.sqrt(2)

# dx = ts[1] - ts[0] # Since rotation preserves distance, the new dx is also the same as the old one

# dBx_p_dx_p, dBx_p_dy_p = np.gradient(Bx_p, dx, dx) # MOT gradient along beam 1
# dBy_p_dx_p, dBy_p_dy_p = np.gradient(By_p, dx, dx) # MOT gradient along beam 2

# grad_G_per_m = dBy_p_dy_p * 1e4  # Tesla/m → G/m
# grad_G_per_cm = grad_G_per_m * 1e-2

# # --- Colormap (magnitude) ---
# pcm = ax.pcolormesh(
#     (X-Y)/np.sqrt(2), (X+Y)/np.sqrt(2), grad_G_per_cm, # Use {(X-Y)/np.sqrt(2), (X+Y)/np.sqrt(2), ((Bx-By)/np.sqrt(2))*1e4} when rotating the axes, otherwise use {X, Y, (Bx)*1e4}
#     shading="auto",
#     cmap="inferno",
# )

# cb = fig.colorbar(pcm, ax=ax)

# # # --- Streamlines (direction) ---
# # ax.streamplot(
# #     X, Y, Bx, By,
# #     color="white",
# #     density=1.2,
# #     linewidth=0.6,
# # )

# # --- Axes styling ---
# ax.set(
#     xlabel="MOT Beam 1 (y=x) (cm)", # Change this to "MOT Beam 1 (y=x) (cm)" in case of rotated axes, otherwise use "x (cm)"
#     ylabel="MOT Beam 2 (y=-x) (cm)", # Change this to "MOT Beam 2 (y=-x) (cm)" in case of rotated axes, otherwise use "y (cm)"
# )
# ax.set_aspect("equal")

# plt.title("dBy/dy [G/cm]") # Change this for different B-field component
# plt.tight_layout()
# plt.show()

# # i0 = len(ts) // 2

# # grad_xp_G_per_m = dBx_p_dx_p[i0, i0] * 1e4   # T/m → G/m
# # grad_xp_G_per_cm = grad_xp_G_per_m * 1e-2   # G/m → G/cm

# # grad_yp_G_per_m = dBy_p_dy_p[i0, i0] * 1e4
# # grad_yp_G_per_cm = grad_yp_G_per_m * 1e-2

# # print("dBx'/dx' at center (G/cm):", grad_xp_G_per_cm)
# # print("dBy'/dy' at center (G/cm):", grad_yp_G_per_cm)

# #mag.show(magnets)
#-------------------------------------To here-------------------------------------#

#------------------------------------- Code for the B-field profile along the x-axis from here -------------------------------------#
xline = np.linspace(-0.2, 0.2, 1000)  # ±10 cm
line = np.array([(x, 0, 0) for x in xline])

B_line = magnets.getB(line)   # shape (N, 3)
Bx_line, By_line, Bz_line = np.moveaxis(B_line, -1, 0)
Bmag = np.sqrt(Bx_line**2 + By_line**2 + Bz_line**2) * 1e4  # T→G

plt.figure()
# plt.plot(xline * 100, Bmag)  # m→cm
# plt.plot(xline * 100, Bx_line*1e4)
plt.plot(xline * 100, By_line*1e4)
# plt.plot(xline * 100, Bz_line*1e4)
plt.xlabel("x (cm)")
plt.ylabel("B (G)")
plt.title("B along physical x-axis (y = 0, z = 0)")
plt.grid()
plt.show()

#------------------------------------- To here -------------------------------------#
