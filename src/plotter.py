import numpy as np
import matplotlib.pyplot as plt

# Load data
with open("image.dat", "r") as f:
    lines = f.readlines()

data = []
for line in lines[1:]:
    if line.strip():
        row = [int(x) for x in line.split()]
        data.extend(row)

image = np.array(data).reshape((256, 256))

plt.figure(figsize=(10, 10))

plt.imshow(image, cmap="gist_heat", origin="lower", interpolation="bilinear")
plt.title("Kerr Black Hole with Accretion Disk (a=0.9)", fontsize=16)
plt.colorbar(label="Intensity")

plt.xlabel("Pixel X")
plt.ylabel("Pixel Y")

plt.tight_layout()
plt.show()

print(f"Image statistics:")
print(f"  Min: {image.min()}, Max: {image.max()}")
print(
    f"  Pixels hitting disk: {np.sum(image > 0)} ({100*np.sum(image > 0)/image.size:.1f}%)"
)
print(f"  Black hole shadow: {np.sum(image == 0)} pixels")
