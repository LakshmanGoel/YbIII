import magpylib as mag
import numpy as np
import matplotlib.pyplot as plt

# ZEEMAN SLOWER COIL PARAMETERS
# All dimensions in mm, current in Amps

coil_params = [
    # z_center, length, inner_radius, outer_radius, turns, current
    (19, 38, 45.009/2, 111/2, 132, 10.0),
    (75, 74, 45.009/2, 84/2, 192, 10.0),
    (136, 48, 45.009/2, 78/2, 160, 10.0),
    (181, 42, 45.009/2, 73/2, 136, 10.0),
    (220, 35, 45.009/2, 67.2/2, 112, 10.0),
    (258, 42, 45.009/2, 62/2, 108, 10.0),
]

# FUNCTION TO BUILD A SOLENOID FROM LOOPS
def make_solenoid(z0, length, r_inner, r_outer, turns, current):
    loops = []

    n_layers = 4  # radial layers
    n_axial = turns // n_layers

    radii = np.linspace(r_inner, r_outer, n_layers)
    z_positions = np.linspace(
        z0 - length/2,
        z0 + length/2,
        n_axial
    )

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
slower = mag.Collection()

for params in coil_params:
    solenoid = make_solenoid(*params)
    slower.add(solenoid)

# import matplotlib.cm as cm

# colors = cm.viridis(np.linspace(0, 1, len(coil_params)))

# slower = mag.Collection()

# for params, color in zip(coil_params, colors):
#     solenoid = make_solenoid(*params)
#     solenoid.style.color = tuple(color[:3])
#     slower.add(solenoid)

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