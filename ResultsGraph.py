import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_displacements_x(numerical_file_path):
    data = pd.read_csv(numerical_file_path, delim_whitespace=True)
    data_filtered = data[data["x1"] == 0]

    data_filtered = data_filtered.drop_duplicates(subset="x0")
    sorted_by_x0 = data_filtered.sort_values(by="x0").reset_index(drop=True)

    P = 0.01
    mu = 7.92
    nu = 0.3
    a = 0.01

    x0 = sorted_by_x0['x0'].tolist()
    u0 = sorted_by_x0['u0'].tolist()

    u0analit = []
    for x in x0:
        u0analit.append(
            (P / (4 * mu)) * (
                    ((1 - nu) / (1 + nu)) * x + (a ** 2) / x +
                    (x + (4 * a ** 2) / ((1 + nu) * x) - (a ** 4) / (x ** 3))
            )
        )
    plt.figure(figsize=(10, 6))
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.plot(x0, u0, color='blue', linewidth=2, label="Numerical solution")
    plt.plot(x0, u0analit, color='green', linewidth=2, label="Analytical solution")
    plt.xlabel("x0")
    plt.ylabel("Displacement")
    plt.legend()
    plt.title("Comparison of Displacement")
    plt.savefig("displacement_down.png")
    plt.show()


def plot_displacements_y(numerical_file_path):
   
    data = pd.read_csv(numerical_file_path, delim_whitespace=True)

    data_filtered = data[data["x0"] == 0]

    data_filtered = data_filtered.drop_duplicates(subset="x1")
    sorted_by_x1 = data_filtered.sort_values(by="x1").reset_index(drop=True)

    P = 0.01
    mu = 7.92
    nu = 0.3
    a = 0.01

    x1 = sorted_by_x1['x1'].tolist()
    u1 = sorted_by_x1['u1'].tolist()

    u1analit = []
    for x in x1:
        u1analit.append(
            (-P / (4 * mu)) * (x + 2 * ((1 - nu) / (1 + nu)) * (a ** 2) / x + (a ** 4) / (x ** 3))
        )
       
    plt.figure(figsize=(10, 6))
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.plot(x1, u1, color='blue', linewidth=2, label="Numerical solution")
    plt.plot(x1, u1analit, color='green', linewidth=2, label="Analytical solution")
    plt.xlabel("x1")
    plt.ylabel("Displacement")
    plt.legend()
    plt.title("Comparison of Displacement")
    plt.savefig("displacement_left.png")
    plt.show()


def sort_numerical_data(file_path):
    data = np.loadtxt(file_path)
    sorted_data = data[data[:, 1].argsort()] 
    sorted_file_path = 'sorted_' + file_path
    np.savetxt(sorted_file_path, sorted_data, fmt='%.6e')
    return sorted_file_path


def kirsch_analytical_solution_y(y_values, radius=0.1, sigma=0.01):
    stresses_y = []
    for y in y_values:
        r = y
        if r > radius:
            stress_y = -sigma * ((radius ** 2 / r ** 2) - (3 * radius ** 4) / (2 * r ** 4))
            stresses_y.append(stress_y)
        else:
            stresses_y.append(0)
    return np.array(stresses_y)


def kirsch_analytical_solution_x(x_values, radius=0.1, sigma=0.01):
    stresses_x = []
    for x in x_values:
        r = x
        if r > radius:
            stress_x = sigma * (1 - (radius ** 2 / r ** 2) + (3 * radius ** 4) / (2 * r ** 4))
            stresses_x.append(stress_x)
        else:
            stresses_x.append(0)
    return np.array(stresses_x)


def plot_stress_left(numerical_file_path):
    sorted_file_path = sort_numerical_data(numerical_file_path)
    data = np.loadtxt(sorted_file_path)
    y_values_numerical = data[:, 0]
    stresses_numerical = data[:, 1]

    unique_indices = np.unique(y_values_numerical, return_index=True)[1]
    y_values_numerical = y_values_numerical[unique_indices]
    stresses_numerical = stresses_numerical[unique_indices]
    stresses_analytical_y = kirsch_analytical_solution_y(y_values_numerical)

    plt.figure(figsize=(10, 6))
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.plot(y_values_numerical, stresses_numerical, color='blue', linewidth=2, label='Numerical solution')
    plt.plot(y_values_numerical, stresses_analytical_y, color='green', linewidth=2, label='Analytical solution')
    plt.xlabel("y")
    plt.ylabel("Stress $\sigma_y$")
    plt.legend()
    plt.title("Comparison of Stress $\sigma_y$ on Left Boundary")
    plt.savefig("stress_left.png")
    plt.show()


def plot_stress_down(numerical_file_path):
    sorted_file_path = sort_numerical_data(numerical_file_path)
    data = np.loadtxt(sorted_file_path)

    x_values_numerical = data[:, 0]
    stresses_numerical = data[:, 1]

    mask = x_values_numerical >= 0.17
    x_values_numerical = x_values_numerical[mask]
    stresses_numerical = stresses_numerical[mask]

    unique_indices = np.unique(x_values_numerical, return_index=True)[1]
    x_values_numerical = x_values_numerical[unique_indices]
    stresses_numerical = stresses_numerical[unique_indices]

    stresses_analytical_x = kirsch_analytical_solution_x(x_values_numerical)

    plt.figure(figsize=(10, 6))
    plt.minorticks_on()
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.plot(x_values_numerical, stresses_numerical, color='blue', linewidth=2, label='Numerical solution')
    plt.plot(x_values_numerical, stresses_analytical_x, color='green', linewidth=2, label='Analytical solution')
    plt.xlabel("x")
    plt.ylabel("Stress $\sigma_x$")
    plt.legend()
    plt.title("Comparison of Stress $\sigma_x$ on Bottom Boundary")

    plt.savefig("stress_down.png")
    plt.show()

plot_stress_left("output_left.txt")
plot_stress_down("output_down.txt")
plot_displacements_x("output_displacements.txt")
plot_displacements_y("output_displacements.txt")
