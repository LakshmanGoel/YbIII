import magpylib as magpy
import numpy as np
import matplotlib as plt

# Goal: Simulate B fields from Zeeman Slower + 2D MOT + 3D MOT
# Zeeman Slower parameters: Get from Prof. Bishof's Mathematica notebooks
# 2D and 3D MOT paramters: Get from Yb II's magpy simulations

# Code for Zeeman Slower B-field:

def zeeman_slower(I=340, config="H", plot_obj=False):

    if config=="H":
        curr_up = -I
        curr_down = -I
    else:
        if config=="AH":
            curr_up = -I
            curr_down = I
        else:
            raise ValueError("Invalid 3D MOT Coil Configuration")

    # Turns (Number of wires along the z axis)
    T = 7
    # Windings (Number of wires on the x-y plane)
    W = 5

    # Spacing between windings
    # s = 5.588 # in mm
    s = 0.00508
    
    # Width of coil
    # L = T * (6) # in mm
    L = T * (s)

    # Diameter
    # d - Inner diameter
    d = 0.10 # 0.27305 # in mm

    # D - Outer diameter
    D = d + 2 * W * s # in mm
    # D = 370.84 # in mm


    # Error
    e = 0 # in mm

    # Position of coils
    z_1 = 0
    z_2 = -z_1

    coil = magpy.Collection()

    for i in range(0, 2*T, 2):
        for n in range(W):

            # Upper Coils
            winding1 = magpy.current.Circle(
                current = curr_up,
                diameter = d + (2*n + 1) * s,
                position = (0, 0, z_1 + (s)*((i - (T-1))/2)),
            )

            coil.add(winding1)

            # Lower Coils
            winding2 = magpy.current.Circle(
                current = curr_down,
                diameter = d + (2*n + 1) * s,
                position = (0, 0, z_2 + (s)*((i - (T-1))/2)),
            )

            coil.add(winding2)

    if plot_obj:
        coil.show(backend='plotly')

    return coil

coil = zeeman_slower(config="H", plot_obj=True)

def plot_plane(plane="xy", min=[-0.050, -0.050], max=[0.050, 0.050], res=100, title="Magnetic Field Streamplot (XZ plane)", save_plot=False, fname="mot3d_ah_i_200a", dir="results\\"):

    fig, axs = plt.subplots(1, 1, figsize=(26,10))

    # create grid
    ax1_ts = np.linspace(min[0], max[0], res)
    ax2_ts = np.linspace(min[0], max[1], res)
    comp = [0, 1]

    if plane=="xy":
        grid = np.array([[(x, y, 0) for x in ax1_ts] for y in ax2_ts])
        comp = [0, 1]
    else:
        if plane=="yz":
            grid = np.array([[(0, y, z) for y in ax1_ts] for z in ax2_ts])
            comp = [1, 2]
        else:
            if plane=="xz":
                grid = np.array([[(x, 0, z) for x in ax1_ts] for z in ax2_ts])
                comp = [0, 2]
            else:
                raise ValueError("Invalid Plane")

    # compute and plot field of coil
    B = magpy.getB(coil, grid)
    Bamp = np.linalg.norm(B, axis=2)
    Bamp /= np.amax(Bamp)

    sp = axs.streamplot(
        grid[:,:,comp[0]], grid[:,:,comp[1]], B[:,:,comp[0]], B[:,:,comp[1]],
        density=2,
        color=np.where(np.linalg.norm(B, axis=2) > 100, 1000, 10*np.linalg.norm(B, axis=2)),
        linewidth=np.sqrt(Bamp)*3,
        cmap='jet',
    )
    # add horizontal center line of 100 mm
    center_y = 0
        # Parameters
    line_length = 0.10          # length along x-axis in meters
    z_positions = [-0.017, 0, 0.017]  # z values for bottom, center, top lines

    x_start = -line_length/2
    x_end = line_length/2

    # Draw lines parallel to x-axis at different z
    for z in z_positions:
        axs.plot([x_start, x_end], [z, z], color='black', linewidth=2, linestyle='--')

    # figure styling
    axs.set(
        title='Magnetic field of coils [T]',
        xlabel=plane[0] + '-position [m]',
        ylabel=plane[1] + '-position [m]',
        aspect=1,
    )

    plt.colorbar(sp.lines, ax=axs, label='[T]')

    if save_plot:
        plt.savefig(dir + fname + '.png', bbox_inches='tight')

    plt.tight_layout()
    plt.show()

plot_plane(plane="xz", fname="mot3d_ah_i_200a_xz", save_plot=False)

# TODO: Write code for 2D MOT
# def mot2d_coil(I=200, config="H", plot_obj="False"):
#     coil_up = magpy.Collection()
#     coil_down = magpy.Collection()
#     return coil_up, coil_down

# TODO: Write code for 3D MOT
# def mot3d_coil(I=200, config="H", plot_obj="False"):
#     coil_up = magpy.Collection()
#     coil_down = magpy.Collection()
#     return coil_up, coil_down


# Archive
# def zeeman_slower():
#     Dwire = 0.00137 # 16 gauge magnet with "heavy" or "dual" insulation
#     Rwire = Dwire/2
#     Layer = 2*Dwire # the height of one "layer" (one winding there and one winding back)
#     R0 = 0.0254 # The outer radius of the slower support. (The smallest it can be is 1.33"/2 plus some thickness. At 1/4" this is like 0.9" or something so lets just make it 1")
#     pushdist = 0.1016 + 0.0127 # Distance from the end of the slower coil to 2D MOT center
#     Rpush = 0.0667 + 0.010 # Where to set the slower field end, add 12.7 mm to center of MOT coil for better fit

#     n = 1/Dwire # Current Density

#     SlowerLength = 0.3794

#     coil = magpy.Collection()

#     for i in range(0, 2*T, 2):
#         for n in range(W):

#             # Upper Coils
#             winding1 = magpy.current.Circle(
#                 current = curr_up,
#                 diameter = d + (2*n + 1) * s,
#                 position = (0, 0, z_1 + (s)*((i - (T-1))/2)),
#             )

#             coil.add(winding1)

#             # Lower Coils
#             winding2 = magpy.current.Circle(
#                 current = curr_down,
#                 diameter = d + (2*n + 1) * s,
#                 position = (0, 0, z_2 + (s)*((i - (T-1))/2)),
#             )

#             coil.add(winding2)

#     if plot_obj:
#         coil.show(backend='plotly')
#     return coil