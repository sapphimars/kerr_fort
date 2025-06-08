import numpy as np
import matplotlib.pyplot as plt

# Read trajectory data
data = np.loadtxt("./trajectory.dat")

print("First few rows of trajectory data:")

tau = data[:, 0]  # proper time
t = data[:, 1]  # coordinate time
r = data[:, 2]  # radius
theta = data[:, 3]  # polar angle
phi = data[:, 4]  # azimuthal angle


# Plot r vs time to see the spiral
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.plot(tau, r)
plt.xlabel(f"Proper time $\\tau$")
plt.ylabel("Radius r")
plt.title("Radial trajectory")
plt.axhline(y=2, color="r", linestyle="--", label="Event horizon (r=2M)")
plt.axhline(y=6, color="g", linestyle="--", label="ISCO (r=6M)")
plt.legend()

plt.subplot(1, 3, 2)
# Convert to Cartesian for orbit plot
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Orbit in x-y plane")
plt.axis("equal")

plt.subplot(1, 3, 3)
plt.plot(tau, theta)
plt.xlabel(f"Proper time $\\tau$")
plt.ylabel(f"$\\theta$ (radians)")
plt.title("Polar angle evolution")

plt.tight_layout()
plt.show()
