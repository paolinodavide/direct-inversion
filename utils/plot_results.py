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
            plot_data(positions, y_data, f"{plot_title} vs Positions", "Position", ylabel, label=file)
    plt.savefig(f"{filetype}.svg", format="svg")
    plt.show()


def plot_error(filetype):
    """Plot error data from a specific file."""
    data = load_data(f"{filetype}.dat")
    if data is not None:
        first_column = data[:, 0]
        last_column = data[:, -1]

        plt.subplot(1, 2, 1)
        plt.scatter(range(len(first_column)), first_column, label="First Column")
        plt.title(f"{filetype} - First Column vs Line Index")
        plt.xlabel("Line Index")
        plt.ylabel("First Column")
        plt.grid(linestyle='--')
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.scatter(range(len(last_column)), last_column, label="Last Column")
        plt.title(f"{filetype} - Last Column vs Line Index")
        plt.xlabel("Line Index")
        plt.ylabel("Last Column")
        plt.grid(linestyle='--')
        plt.legend()

        plt.tight_layout()
        plt.savefig(f"{filetype}.svg", format="svg")
        plt.show()


def plot_final():
    """Plot the final data for 'gr', 'pot', and 'forces'."""
    for i, (filetype, ylabel, title) in enumerate([('gr', 'g(r)', 'g(r)'),
                                                   ('pot', 'Potential', 'Potential'),
                                                   ('pot', 'Force', 'Force')], start=1):
        plt.subplot(1, 3, i)
        files = glob.glob(f"{filetype}_*.dat")
        files = [file for file in files if file.split('_')[-1].split('.')[0].isdigit()]
        files = sorted(files, key=lambda x: int(x.split('_')[-1].split('.')[0]), reverse=True)[:1]
        target_file = f"{filetype}_target.dat"
        
        data = load_data(files[-1])
        if data is not None:
            positions = data[:, 0]
            y_data = data[:, 1] if ylabel != 'Force' else data[:, 2]
            plot_data(positions, y_data, f"{title} vs Positions", "Position", ylabel, label=files[-1], style='o')

        target_data = load_data(target_file)
        if target_data is not None:
            target_positions = target_data[:, 0]
            target_y_data = target_data[:, 1] if ylabel != 'Force' else target_data[:, 2]
            plot_data(target_positions, target_y_data, f"{title} vs Positions", "Position", ylabel, label=target_file, style='--')
    plt.tight_layout()
    plt.show()


def main():
    try:
        filetype = sys.argv[1]
    except IndexError:
        print("Enter the name of the file to parse")
        return

    column_index = 1
    if filetype == 'forces':
        filetype = 'pot'
        column_index = 2

    if filetype not in ['err', 'final']:
        process_files(filetype, column_index, filetype, filetype)
    elif filetype == 'err':
        plot_error(filetype)
    elif filetype == 'final':
        plot_final()


if __name__ == "__main__":
    main()


