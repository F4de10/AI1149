from exercise_methods import *
from scipy.stats import chi2, norm, t

print("\n", "problem 5.1", "\n")
A = np.array([[1, 0], [0, 1], [-1, -1]])
P = np.diag([1, 1, 1])
Q_ee = 1 / 3 * np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])

ic(Q_ee, P)
calculate_redundancy(A, P, Q_ee)
print("\n")

A = np.array([[1, 0], [0, 1], [-1, -1]])
P = np.diag([2, 2, 1])
Q_ee = 1 / 8 * np.array([[1, 1, 2], [1, 1, 2], [2, 2, 4]])
ic(Q_ee, P)
calculate_redundancy(A, P, Q_ee)

print("\n", "problem 5.2", "\n")

observations = np.array([[1.680], [-0.148], [0.557]])
P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
A = np.array([[1], [1], [1]])
C = np.array([[-10.421], [-12.261], [-11.559]])

Q_ee = adjustment_by_elements(observations, A, P, C, "meter", "false")
ic(Q_ee, P)
R, r = calculate_redundancy(A, P, Q_ee)

print("\n", "problem 5.3", "\n")

observations = np.array([[100.006], [99.999], [200.002]])


P = np.diag([1, 1, 1])
A = np.array([[1, 0], [0, 1], [1, 1]])
ic(A)
C = np.array([[0], [0], [0]])
ic(C)
Q_ee = adjustment_by_elements(observations, A, P, C, "meter", "false")
ic(Q_ee, P)
calculate_redundancy(A, P, Q_ee)
print("\n")

P = np.diag([4, 4, 1])

Q_ee = adjustment_by_elements(observations, A, P, C, "meter", "false")
ic(Q_ee, P)
calculate_redundancy(A, P, Q_ee)

print("\n", "problem 5.4", "\n")
P = np.diag([2, 2, 1, 1, 2])
L = np.array([[1.002], [2.004], [-12.001], [8.998], [3.012]])
A = np.array([[-1, 1, 0], [0, -1, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 1]])
epsilon = np.array([[-1.5], [-1.5], [3.0], [3.0], [3.0]])
ic(epsilon)
Q_ee = (
    1
    / 28
    * np.array(
        [
            [5, 5, 2, 2, -4],
            [5, 5, 2, 2, -4],
            [2, 2, 12, 12, 4],
            [2, 2, 12, 12, 4],
            [-4, -4, 4, 4, 6],
        ]
    )
)
variance_factor = 9
estimated_variance = 22.5


print("\n", "Overall test:", "\n")
n = len(epsilon)
m = len(A[0])
ic(n - m)
ic(chi2.ppf(1 - 0.05, n - m))
critical_value = variance_factor * chi2.ppf(1 - 0.05, n - m)
omega, p_value = chi_square_test(epsilon, P, variance_factor)
ic(p_value)
ic(critical_value)

print("\n", "Individual u-test:", "\n")
u = np.zeros((len(epsilon), 1))
critical_value = np.zeros((len(epsilon), 1))
for i in range(len(epsilon)):
    u[i][0] = epsilon[i][0] / (np.sqrt(variance_factor) * np.sqrt(Q_ee[i][i]))
    critical_value = norm.ppf(1 - 0.025)
ic(u)
ic(critical_value)

print("\n", "Individual t-test:", "\n")
n = len(epsilon)
m = len(A[0])
w = np.zeros((len(epsilon), 1))
critical_value = np.zeros((len(epsilon), 1))
for i in range(len(epsilon)):
    w[i][0] = epsilon[i][0] / (np.sqrt(estimated_variance) * np.sqrt(Q_ee[i][i]))
    critical_value = t.ppf(1 - 0.025, n - m - 1)
ic(w)
ic(critical_value)

print("\n", "problem 5.5", "\n")

print("\n", "first case:", "\n")
A = np.array([[1, 0], [0, 1], [-1, -1]])
P = np.diag([1, 1, 1])
Q_ee = 1 / 3 * np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])

ic(Q_ee, P)
R, r = calculate_redundancy(A, P, Q_ee, "false")
print("\n")
ic(r)

variance_factor = 22.5
variance = np.zeros((len(r), 1))
for i in range(len(P)):
    variance[i] = np.sqrt(variance_factor) / np.sqrt(P[i][i])
ic(variance)

theta = np.zeros((len(variance), 1))
phi = np.zeros((len(variance), 1))
for i in range(len(variance)):
    theta[i] = 2.8 * (variance[i] / np.sqrt(r[i]))
    phi[i] = (1 - r[i]) * theta[i]
ic(theta)
ic(phi)

print("\n", "second case:", "\n")

A = np.array([[1, 0], [0, 1], [-1, -1]])
P = np.diag([2, 2, 1])
Q_ee = 1 / 8 * np.array([[1, 1, 2], [1, 1, 2], [2, 2, 4]])
ic(Q_ee, P)
R, r = calculate_redundancy(A, P, Q_ee, "false")
ic(r)
print("\n")

variance_factor = 22.5
variance = np.zeros((len(r), 1))
for i in range(len(P)):
    variance[i] = np.sqrt(variance_factor) / np.sqrt(P[i][i])
ic(variance)

theta = np.zeros((len(variance), 1))
phi = np.zeros((len(variance), 1))
for i in range(len(variance)):
    theta[i] = 2.8 * (variance[i] / np.sqrt(r[i]))
    phi[i] = (1 - r[i]) * theta[i]
ic(theta)
ic(phi)
