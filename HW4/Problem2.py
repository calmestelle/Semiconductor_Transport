import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbols
MA, MB = sp.symbols('M_A M_B', real=True)
gamma1, gamma2 = sp.symbols('gamma_1 gamma_2', real=True)
omega, a, k = sp.symbols('omega a k', real=True)
x = sp.symbols('x', real=True)  # x = omega^2
i = sp.I

# Parameter values
params = {
    gamma1: 1.0,
    gamma2: 1.0,
    MA: 10.0,
    MB: 1.0,
    a: 1.0
}
a_val = params[a]

# Define symbolic matrix
M = sp.Matrix([
    [-MA * omega**2 + gamma2 + gamma1, -gamma1, -gamma2 * sp.exp(-i * k * a)],
    [-gamma1, -MB * omega**2 + 2 * gamma1, -gamma1],
    [-gamma2 * sp.exp(i * k * a), -gamma1, -MA * omega**2 + gamma1 + gamma2]
])

# Compute and simplify determinant
detM = sp.simplify(M.det()).subs(omega**2, x)
detM_sub = detM.subs(params)

# Extract coefficients as functions of k
coeffs_sym = sp.Poly(detM_sub, x).all_coeffs()
coeffs_sym = [sp.simplify(c) for c in coeffs_sym]
coeff_funcs = [sp.lambdify(k, c, modules='numpy') for c in coeffs_sym]

# Solve cubic in x = ω² for each k
k_vals = np.linspace(0, np.pi / a_val, 200)
omega_branches = [[], [], []]

for k_val in k_vals:
    coeff_vals = [f(k_val) for f in coeff_funcs]
    roots = np.roots(coeff_vals)
    real_omegas = [np.sqrt(r.real) for r in roots if np.isreal(r) and r.real >= 0]
    real_omegas.sort()
    for i in range(3):
        omega_branches[i].append(real_omegas[i] if i < len(real_omegas) else np.nan)

omega_branches = np.array(omega_branches)

# Plotting
plt.figure(figsize=(8, 5))
for i in range(3):
    label = f'Optical Branch {i}' if i > 0 else 'Acoustic'
    plt.plot(k_vals, omega_branches[i], label=label)

plt.xlabel(r'$k\ (\pi/a)$')
plt.ylabel(r'$\omega$ (rad/s)')
plt.title(r'Phonon Dispersion: $\omega(k)$ for 3-atom unit cell')
plt.grid(True)
plt.legend()

# Add parameter text box
param_text = (f"$\gamma_1/\gamma_2 = {params[gamma1]/params[gamma2]:.1f}$\n"
              f"$M_A = {params[MA]}$\n"
              f"$M_B = {params[MB]}$")
plt.text(0.05, 0.95, param_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig("phonon_dispersion.pdf", format='pdf', bbox_inches='tight')
plt.show()
