import numpy as np
import matplotlib.pyplot as plt

from atomsmltr.environment import PlaneWaveLaserBeam
from atomsmltr.environment import ConstantForce
from atomsmltr.atoms import Ytterbium
from atomsmltr.simulation import Configuration, RK4

# - setup atom
atom = Ytterbium()
main = atom.trans["main"] # get transition, to help setting up lasers
intercombination = atom.trans["intercombination"]

# - setup laser
laser_1 = PlaneWaveLaserBeam()
laser_1.direction = (0, -np.sin(45), np.cos(45))
laser_1.set_power_from_I(intercombination.Isat) # set power to reach Isat
laser_1.tag = "las1"

# - config
config = Configuration()
config.atom = atom
config += laser_1
config.add_atomlight_coupling("las1", "intercombination", 0*intercombination.Gamma) # Arguments: laser = "las1", transition = intercombination", detuning = - or +200 * intercombination.Gamma
#config.add_atomlight_coupling("las2", "green", 1 * intercombination.Gamma)

# let's add gravity, pointing along +y
m = Ytterbium().mass  # kg
g = 9.81  # m/s^2
grav_force = (0,  m * g, 0)
gravity = ConstantForce(field_value=grav_force, tag="gravity")
config += gravity

# - simulation
sim = RK4(config=config)
t = np.linspace(0, 0.25, 6000) # timesteps for integration
u0 = (0, 0, 0, 0, 0, 0) # atom starts with vz=100m/s
res = sim.integrate(u0, t)

# plot
fix, axes = plt.subplots(3, 2, tight_layout=True)

axes[0, 0].plot(res.t * 1e3, res.y[0])
axes[0, 0].set_ylabel("x␣(m)")
axes[0, 1].plot(res.t * 1e3, res.y[3])
axes[0, 1].set_ylabel("vx␣(m/s)")

axes[1, 0].plot(res.t * 1e3, res.y[1])
axes[1, 0].set_ylabel("y␣(m)")
axes[1, 1].plot(res.t * 1e3, res.y[4])
axes[1, 1].set_ylabel("vy␣(m/s)")

axes[2, 0].plot(res.t * 1e3, res.y[2])
axes[2, 0].set_ylabel("z␣(m)")
axes[2, 1].plot(res.t * 1e3, res.y[5])
axes[2, 1].set_ylabel("vz␣(m/s)")

for ax in axes.flat:
    ax.set(xlabel="t␣(ms)")

plt.show()
