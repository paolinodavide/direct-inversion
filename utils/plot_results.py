import numpy as np
import matplotlib.pyplot as plt
import sys
import glob


def load_data(file, skiprows=1):
    """Load data from a file, skipping the first row."""
    try:
        return np.loadtxt(file, skiprows=skiprows)
    except OSError:
        print(f"File {file} not found. Skipping.")
        return None


def plot_data(x, y, title, xlabel, ylabel, label=None, style='-', markersize=3):
    """Plot data with the given parameters."""
    plt.plot(x, y, style, label=label, markersize=markersize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(linestyle='--')
    if label:
        plt.legend()


def process_files(filetype, column_index, plot_title, ylabel):
    """Process and plot data from files matching the given filetype."""
    files = glob.glob(f"{filetype}_*.dat")
    for file in files:
        data = load_data(file)
        if data is not None:
            positions = data[:, 0]
            y_data = data[:, column_index]
            plot_data(positions, y_data, f"{plot_title} vs Positions", r"$r / \sigma$", ylabel, label=file)
    plt.savefig(f"{filetype}.svg", format="svg")
    plt.show()


def plot_error(filetype):
    """Plot error data from a specific file."""
    data = load_data(f"{filetype}.dat")
    if data is not None:
        first_column = data[:, 0]
        last_column = data[:, -1]

        plt.subplot(1, 2, 1)
        plt.semilogy(range(len(first_column)), first_column, 'o-', label=r"$\Delta(g)^2$")
        plt.title("Delta(g)^2 vs Iteration Number")
        plt.xlabel("Iteration")
        plt.ylabel("Error")
        plt.grid(linestyle='--')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.semilogy(range(len(last_column)), last_column, 'o-', label=r"$\Delta(\text{iter})^2$")
        plt.title("Delta(iter)^2 vs Iteration Number")
        plt.xlabel("Iteration")
        plt.ylabel("Err")
        plt.grid(linestyle='--')
        plt.legend()

        plt.tight_layout()
        plt.savefig(f"{filetype}.svg", format="svg")
        plt.show()


def plot_final():
    """Plot the final data for 'gr', 'pot', and 'forces'."""
    plt.figure(figsize=(16, 9))
    for i, (filetype, ylabel, title) in enumerate([('gr', 'g(r)', 'g(r)'),
                                                   ('pot', 'Potential', 'Potential'),
                                                   ('pot', 'Force', 'Force')], start=1):
        plt.subplot(1, 3, i)
        files = glob.glob(f"{filetype}_*.dat")
        files = [file for file in files if file.split('_')[-1].split('.')[0].isdigit()]
        files = sorted(files, key=lambda x: int(x.split('_')[-1].split('.')[0]), reverse=True)[:1]
        target_file = f"{filetype}_target.dat"
        
        # Converged results
        data = load_data(files[-1])
        if data is not None:
            positions = data[:, 0]
            y_data = data[:, 1] if ylabel != 'Force' else data[:, 2]
            plot_data(positions, y_data, f"{title} vs Positions", r"$r/\sigma$", ylabel, label=files[-1], style='o')

        # Target data
        target_data = load_data(target_file)
        if target_data is not None:
            target_positions = target_data[:, 0]
            target_y_data = target_data[:, 1] if ylabel != 'Force' else target_data[:, 2]
            plot_data(target_positions, target_y_data, f"{title} vs Positions", r"$r/\sigma$", ylabel, label=target_file, style='--')
            if ylabel != 'g(r)':
                MSE = np.linalg.norm(target_y_data - y_data)**2 / len(target_y_data)
                SNR = 10 * np.log10(np.sum(target_y_data**2) / np.sum((target_y_data - y_data)**2))
                plt.text(0.05, 0.95, f"MSE: {MSE:.4e}\nSNR: {SNR:.2f} dB", 
                         transform=plt.gca().transAxes, fontsize=10, 
                         verticalalignment='top', bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout()
    plt.show()


def main():
    try:
        filetype = sys.argv[1]
    except IndexError:
        print("Enter the name of the file to parse")
        return

    column_index = 1
    plot_type = filetype
    if filetype == 'forces':
        filetype = 'pot'
        column_index = 2

    if filetype not in ['err', 'final', 'all']:
        process_files(filetype, column_index, plot_type, plot_type)
    elif filetype == 'err':
        plot_error(filetype)
    elif filetype == 'final':
        plot_final()
    elif filetype == 'all':
        process_files('gr', 1, 'g(r)', 'g(r)')
        process_files('pot', 1, 'Potential', 'Potential')
        process_files('pot', 2, 'Forces', 'Force')
        plot_error('err')
        plot_final()


if __name__ == "__main__":
    main()


