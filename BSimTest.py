import magpylib as mag
import numpy as np
import matplotlib.pyplot as plt

# ZEEMAN SLOWER COIL PARAMETERS
# All dimensions in mm, current in Amps

coil_params = [
    # z_center, length, inner_radius, outer_radius, current
    (19, 38, 84/2, 111/2, 4.0),
    (56, 112, 78/2, 84/2, 4.0),
    (80, 160, 73/2, 78/2, 4.0),
    (101, 202, 67.2/2, 73/2, 4.0),
    (118.5, 237, 62/2, 67.2/2, 4.0),
    (139.5, 279, 45.009/2, 62/2, 2.5),
]

# FUNCTION TO BUILD A SOLENOID FROM LOOPS
def make_solenoid(z0, length, r_inner, r_outer, current):
    loops = []
    d_wire = 1.37 # Diameter of wires (in mm)
    n_radial = np.floor((r_outer - r_inner)/d_wire).astype(int)  # radial layers
    n_axial = np.floor((length)/d_wire).astype(int) # axial layers

    # radii = np.linspace(r_inner, r_outer, n_radial)
    # z_positions = np.linspace(
    #     z0 - length/2,
    #     z0 + length/2,
    #     n_axial
    # )

    radii = r_inner + (np.arange(n_radial) + 0.5)*d_wire
    z_positions = z0 - length/2 + (np.arange(n_axial) + 0.5)*d_wire

    for r in radii:
        for z in z_positions:
            loop = mag.current.Circle(
                current=current,
                diameter=2*r,
                position=(0, 0, z)
            )
            loops.append(loop)

    return mag.Collection(loops)

# BUILD ZEEMAN SLOWER

import matplotlib.cm as cm

colors = cm.viridis(np.linspace(0, 1, len(coil_params)))

slower = mag.Collection()

for params, color in zip(coil_params, colors):
    solenoid = make_solenoid(*params)
    solenoid.style.color = tuple(color[:3])
    slower.add(solenoid)

# FIELD CALCULATION
z = np.linspace(-300, 600, 10000)
positions = np.column_stack([np.zeros_like(z),
                              np.zeros_like(z),
                              z])

B = slower.getB(positions)
Bz = B[:, 2]  # axial field

# PLOT
plt.figure(figsize=(9,5))
plt.plot(z, Bz*1e4)  # Tesla → Gauss
plt.xlabel("z (mm)")
plt.ylabel("Bz (Gauss)")
plt.title("Zeeman Slower Magnetic Field")
plt.grid()
plt.tight_layout()
plt.show()

mag.show(
    slower,
    backend="plotly",
    show_axes=True,
)

# # === Create the individual coils separately ===
# def mot2d_individual(I=200, config="H"):
#     if config == "H":
#         curr_up = I
#         curr_down = I
#     elif config == "AH":
#         curr_up = I
#         curr_down = -I
#     else:
#         raise ValueError("Invalid config")

#     T = 1
#     W = 1
#     s = 0.00
#     L = T * s
#     d = 0.130
#     D = d + 2 * W * s
#     z_1 = 279 + 250
#     z_2 = -z_1

#     coil_up = mag.Collection()
#     coil_down = mag.Collection()

#     for i in range(0, 2*T, 2):
#         for n in range(W):
#             winding_up = mag.current.Circle(
#                 current=curr_up,
#                 diameter=d + (2*n + 1)*s,
#                 position=(0, 0, z_1 + (s)*((i - (T-1))/2)),
#             )
#             winding_down = mag.current.Circle(
#                 current=curr_down,
#                 diameter=d + (2*n + 1)*s,
#                 position=(0, 0, z_2 + (s)*((i - (T-1))/2)),
#             )
#             coil_up.add(winding_up)
#             coil_down.add(winding_down)

#     return coil_up, coil_down


# # === Compute and plot the field ===
# I = 210
# config = "H"  # or "H" for Helmholtz-like
# coil_up, coil_down = mot3d_individual(I=I, config=config)

# # Total field (combined system)
# total_coil = mag.Collection(coil_up, coil_down)

# # Sampling points along z-axis
# z = np.linspace(-0.1, 0.1, 300)
# B_up = np.array([coil_up.getB((0, 0, zi)) for zi in z])
# B_down = np.array([coil_down.getB((0, 0, zi)) for zi in z])
# B_total = np.array([total_coil.getB((0, 0, zi)) for zi in z])

# # Extract z-components
# Bz_up = B_up[:, 2]
# Bz_down = B_down[:, 2]
# Bz_total = B_total[:, 2]

# # === Plot ===
# plt.figure(figsize=(7,5))
# plt.plot(z*1000, Bz_up*1e4, label='Upper Coil', color='tab:blue', linestyle='--')
# plt.plot(z*1000, Bz_down*1e4, label='Lower Coil', color='tab:orange', linestyle='--')
# plt.plot(z*1000, Bz_total*1e4, label='Total Field', color='black', linewidth=2)
# plt.axhline(0, color='gray', linestyle=':')
# plt.xlabel('z [mm]')
# plt.ylabel('Bz [G]')  # Tesla → Gauss for readability
# plt.title(f'Axial Magnetic Field — 3D MOT Coils ({config} config)')
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
