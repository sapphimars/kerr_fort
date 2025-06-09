import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Load data
with open("image.dat", "r") as f:
    lines = f.readlines()

data = []
for line in lines[1:]:
    if line.strip():
        row = [float(x) for x in line.split()]  # Note: float now
        data.extend(row)

image = np.array(data).reshape((768, 1024))

# Create enhanced visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))

# Main image with enhanced contrast
im1 = ax1.imshow(
    image, cmap="hot", origin="lower", vmin=0, vmax=np.percentile(image[image > 0], 95)
)
ax1.set_title("Kerr Black Hole (a=0.9)")
ax1.set_xlabel("Pixel X")
ax1.set_ylabel("Pixel Y")
plt.colorbar(im1, ax=ax1, label="Intensity")

# Contour plot to show structure
im2 = ax2.contourf(image, levels=20, cmap="plasma")
ax2.set_title("Contour View - Structure Analysis")
ax2.set_xlabel("Pixel X")
ax2.set_ylabel("Pixel Y")
plt.colorbar(im2, ax=ax2, label="Intensity")

# Add annotations for expected features
ax1.text(
    0.02,
    0.98,
    "Expected:\n• Central shadow\n• Photon sphere\n• Primary disk\n• Secondary image",
    transform=ax1.transAxes,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
)

# plt.tight_layout()
plt.show()

print(f"Image statistics:")
print(f"  Min: {image.min():.3f}, Max: {image.max():.3f}")
print(f"  Photon sphere candidates: {np.sum((image > 1.2) & (image < 3.0))} pixels")
print(f"  Primary disk: {np.sum(image > 3.0)} pixels")
