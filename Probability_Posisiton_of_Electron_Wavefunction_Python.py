#SUMUR POTENSIAL

import numpy as np
import matplotlib.pyplot as plt

# Metode Simpson 1/3
def simpson_1_3(func, a, b, n):
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's 1/3 rule.")

    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = func(x)

    integral = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2])
    integral *= h / 3
    return integral

# Fungsi gelombang sumur potensial
def psi(x, n, L):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)  # Normalized wave function

def probability_density(x, n, L):
    return psi(x, n, L) ** 2

# Parameter
L = 1.0  # Lebar sumur
n = 11    # Orde

# Daerah survei
a = L/2
b = (3/4) * L

# jumlahan iterasi
n_intervals = 100

# Perhitungan probabilitas
P = simpson_1_3(lambda x: probability_density(x, n, L), a, b, n_intervals)

print(f"Probability of finding the particle in the interval [{a}, {b}]: {P:.6f}")

# Plot kurva
x_vals = np.linspace(0, L, 500)
y_wave = psi(x_vals, n, L)
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_wave, label=f'$\psi(x)$ for n={n}', color='blue')
plt.fill_between(x_vals, y_wave, where=(x_vals >= a) & (x_vals <= b), color='orange', alpha=0.5, label=f'Probability region in [{a}, {b}]')
plt.title("Wave Function for Infinite Potential Well")
plt.xlabel("x")
plt.ylabel("$\psi(x)$")
plt.legend()
plt.grid()
plt.show()

#OSILATOR HARMONIK

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite, factorial

'''
Fungsi Gelombang pada osilator harmonik memerlukan fungsi hermite
'''

# Fungsi gelombang osilator harmonik
def psi_n(x, n, m, omega, hbar):
    hermite_poly = hermite(n)
    normalization = 1.0 / np.sqrt(2**n * factorial(n)) * (m * omega / (np.pi * hbar))**0.25
    gaussian = np.exp(-m * omega * x**2 / (2 * hbar))
    return normalization * gaussian * hermite_poly(np.sqrt(m * omega / hbar) * x)

def probability_density(x, n, m, omega, hbar):
    return psi_n(x, n, m, omega, hbar)**2

# Metode Simpson 1/3
def simpson_1_3(func, a, b, n):
    if n % 2 != 0:
        raise ValueError("n must be even for Simpson's 1/3 rule.")

    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = func(x)

    integral = y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2])
    integral *= h / 3
    return integral

# Parameter gelombang
m = 9.11e-31    # massa elektron
omega = 1       # frekuensi sudut (anggap)
hbar = 1.05e-34 # konstanta planck reducted
n = 3           # orde

# Daerah survey
a = -1e-2
b = 1.5e-2

# Jumlah
n_intervals = 100

# Membuat fungsi probabilitas
func = lambda x: probability_density(x, n, m, omega, hbar)

# Menghitung probabilitas
P = simpson_1_3(func, a, b, n_intervals)

print(f"Probability of finding the particle in the interval [{a}, {b}]: {P:.6f}")

# Plot kurva
x_vals = np.linspace(-1e-1, 1e-1, 500)
y_wave = psi_n(x_vals, n, m, omega, hbar)
y_density = probability_density(x_vals, n, m, omega, hbar)

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(x_vals, y_wave, label=f'$\psi_{n}(x)$', color='blue')
plt.title("Wave Function for Quantum Harmonic Oscillator")
plt.xlabel("x")
plt.ylabel("$\psi_{n}(x)$")
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(x_vals, y_density, label=f'$|\psi_{n}(x)|^2$', color='green')
plt.fill_between(x_vals, y_density, where=(x_vals >= a) & (x_vals <= b), color='orange', alpha=0.5, label=f'Probability region in [{a}, {b}]')
plt.title("Probability Density for Quantum Harmonic Oscillator")
plt.xlabel("x")
plt.ylabel("$|\psi_{n}(x)|^2$")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
