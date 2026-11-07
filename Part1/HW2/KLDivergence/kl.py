import numpy as np
import matplotlib.pyplot as plt

mu1 = np.linspace(-5, 5, 200)
mu2 = np.linspace(-5, 5, 200)
M1, M2 = np.meshgrid(mu1, mu2)

sigmas = [0.5, 1.0, 2.0]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, sig in enumerate(sigmas):
    # KL(P || Q) with equal variances simplifies to (mu1 - mu2)^2 / (2 * sigma^2)
    KL = (M1 - M2)**2 / (2 * sig**2)
    Z_clipped = np.clip(KL, 0, 10) # Cap at 10 for consistent colormapping
    
    ax = axes[i]
    levels = np.linspace(0, 10, 40)
    cp = ax.contourf(M1, M2, Z_clipped, levels=levels, cmap='inferno')
    
    # Line where means are equal
    ax.plot([-5, 5], [-5, 5], 'w--', alpha=0.5, label=r'$\mu_1 = \mu_2$')
    
    ax.set_title(rf'$\sigma = {sig}$', fontsize=14)
    ax.set_xlabel(r'$\mu_1$ (Mean of P)', fontsize=12)
    if i == 0:
        ax.set_ylabel(r'$\mu_2$ (Mean of Q)', fontsize=12)
    if i == 2:
        ax.legend(loc='upper left')

cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
cbar = fig.colorbar(cp, cax=cbar_ax)
cbar.set_label('KL Divergence', fontsize=12)

plt.suptitle(r'Kullback-Leibler Divergence $D_{KL}(P \parallel Q)$ with Equal Variances ($\sigma_1 = \sigma_2 = \sigma$)', fontsize=16, y=1.05)
plt.subplots_adjust(right=0.9, wspace=0.15)
plt.savefig('kl2.png')
