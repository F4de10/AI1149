import numpy as np
import scipy.stats

def chi_square_test(epsilon, P, variance_factor):
    omega = np.dot(np.dot(epsilon.T, P), epsilon)
    chi = omega / (variance_factor**2)
    p_value = 1 - scipy.stats.chi2.cdf(chi, len(epsilon))
    return p_value

# Test case 1
epsilon = np.array([[-1.5], [-1.5], [3.0], [3.0], [3.0]])
P = np.diag([4, 4, 1, 1, 2])
Q_ee = 1 / 28 * np.array([[5, 5, 2, 2, -4], [5, 5, 2, 2, -4], [2, 2, 12, 12, 4], [2, 2, 12, 12, 4], [-4, -4, 4, 4, 6]])
variance_factor = 9

p_value = chi_square_test(epsilon, P, Q_ee, variance_factor)
print("Test case 1: p-value =", p_value)

# Test case 2
epsilon = np.array([[0.5], [0.5], [-1.0], [-1.0], [-1.0]])
P = np.diag([1, 1, 1, 1, 1])
Q_ee = 1 / 10 * np.array([[2, 2, 1, 1, -2], [2, 2, 1, 1, -2], [1, 1, 5, 5, 2], [1, 1, 5, 5, 2], [-2, -2, 2, 2, 3]])
variance_factor = 5

p_value = chi_square_test(epsilon, P, Q_ee, variance_factor)
print("Test case 2: p-value =", p_value)