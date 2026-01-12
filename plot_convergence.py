import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# --- Configuration ---
Path = './outputs/'
iterations_prefix = 'iteration_'

# 1. Get and sort files
if not os.path.exists(Path):
    os.makedirs(Path)

files = [f for f in os.listdir(Path) if f.startswith(iterations_prefix)]
files = sorted(files, key=lambda x: float(x.split('_')[1].split('.')[0]) if x.split('_')[1].split('.')[0].lstrip('-').isdigit() else float('inf'))

num_files = len(files)
colors = cm.coolwarm(np.linspace(0, 1, num_files))

# 2. Setup Figure with custom width ratios
# Ratio explanation: [Panel A, Panel B, Colorbar Slot, Panel C]
# Setting the 3rd element to 0.05 keeps B and C apart for the bar, 
# while a and b stay closer because there is no slot between them.
fig, axs_all = plt.subplots(1, 4, figsize=(8, 8/3), 
                            gridspec_kw={'width_ratios': [1, 1, 0.05,  1]},
                            constrained_layout=True)

# Assign readable names to our axes
ax_a, ax_b, ax_cbar, ax_c = axs_all
axs = [ax_a, ax_b, ax_c] # List for easy looping

# 3. Plotting with Conditional Transparency
for i, file in enumerate(files):
    try:
        radii, gr, pot, _ = np.loadtxt(os.path.join(Path, file), unpack=True)
        
        if i in [1, num_files - 1]: 
            alpha_val, lw, z, label_val = 1.0, 1.4, 10, None
        elif file == 'iteration_-1.dat':
            colors[i] = cm.seismic(0)  
            alpha_val, lw, z, label_val = 1.0, 1.4, 11, r'Potential $r^{-3}$'
        else:
            alpha_val, lw, z, label_val = 0.6, 0.7, 1, None

        if file != 'iteration_-1.dat':
            ax_a.plot(radii, gr, color=colors[i], alpha=alpha_val, linewidth=1.4, zorder=z)
        
        ax_b.plot(radii, pot, color=colors[i], alpha=alpha_val, linewidth=lw, zorder=z, label=label_val, linestyle='-' if file != 'iteration_-1.dat' else '--')
    except Exception as e:
        print(f"Skipping {file}: {e}")

# 4. Colorbar in the reserved middle slot
# 4. Colorbar in the reserved middle slot
sm = plt.cm.ScalarMappable(cmap=cm.coolwarm, norm=plt.Normalize(vmin=0, vmax=num_files-1))
cbar = fig.colorbar(sm, cax=ax_cbar, orientation='vertical')

# Remove the previous side-label
cbar.set_label('') 

#ctick on the left
ax_cbar.set_xlabel(r'$t$', labelpad=4, ha='left')

# 5. Convergence Data (Panel C)
try:
    it, _, err, it_diff, _, _ = np.loadtxt(os.path.join(Path, 'convergence_data.dat'), unpack=True)
    ax_c.semilogy(it, err, label=r'MSE$(g_t, g_{\text{ref}})$', color=colors[25], linewidth=2)
    ax_c.semilogy(it, it_diff, label=r'MSE$(g_t, g_{t-1})$', color=colors[num_files*2//3], linewidth=2     )
    ax_c.legend(fontsize='small', frameon=True)
    ax_c.set_xlim(0, max(it))
except:
    ax_c.set_title('Convergence Data Not Found', fontsize=10)

# 6. Formatting & Panel Labels
labels = ['(a)', '(b)', '(c)']
separations = [-0.18, -0.15, -0.23]  # Custom separations for panel labels
for i, ax in enumerate(axs):
    ax.tick_params(direction='in', which='both', top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    
    ax.text(separations[i], 1.05, labels[i], transform=ax.transAxes, 
            fontsize=11, fontweight='bold', va='top', ha='right')
    

# Panel A Limits/Labels
ax_a.set_xlabel(r'$r\, /\, \sigma$')
ax_a.set_ylabel(r'$g_t(r)$')

# Panel B Limits/Labels
ax_b.set_xlabel(r'$r\, /\, \sigma$')
ax_b.set_ylabel(r'$\beta u_t(r)$')


ax_b.legend(loc='upper right')

# Panel C Labels
ax_c.set_xlabel(r'Iteration $t$')
ax_c.set_ylabel('Convergence Metrics')

# Save and Show
plt.savefig(os.path.join(Path, '00Convergence_Highlight.pdf'), dpi=300)
plt.show()