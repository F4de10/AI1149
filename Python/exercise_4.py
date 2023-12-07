from exercise_methods import *

print("problem 4", "\n")
l1 = 60 + 0 / 60 + 00.2 / 3600
l2 = 70 + 53 / 60 + 38.2 / 3600
l3 = 49 + 6 / 60 + 25.0 / 3600
l4 = 29 + 59 / 60 + 57.7 / 3600
l5 = 3464.121

observations = np.array([[l1], [l2], [l3], [l4], [l5]])
# 2 parameters: P(x, y)

variance = np.array([[3], [3], [3], [3], [1.5]])
variance_factor = 3

P = np.zeros((len(variance), len(variance)))
for i in range(len(variance)):
    P[i][i] = variance_factor**2 / variance[i][0] ** 2
print(P)

given_coordinates = np.array(
    [[2732.651, 1000.314], [1000.600, 2000.314], [1000.600, 6500.314]]
)

x_A, y_A = given_coordinates[0, 0], given_coordinates[0, 1]
x_B, y_B = given_coordinates[1, 0], given_coordinates[1, 1]
x_C, y_C = given_coordinates[2, 0], given_coordinates[2, 1]

azimuth_CB = azimuth(x_C, y_C, x_B, y_B)

print("azimuths_CB: ", azimuth_CB)

azimuth_CP = azimuth_CB + l4
print("azimuths_CP: ", azimuth_CP)

x_P_0 = x_C + l5 * np.cos(np.radians(azimuth_CP))
y_P_0 = y_C + l5 * np.sin(np.radians(azimuth_CP))
print("x_P_0: ", x_P_0)
print("y_P_0: ", y_P_0)

AP_a = azimuth(x_A, y_A, x_P_0, y_P_0)
AB_a = azimuth(x_A, y_A, x_B, y_B)

BA_a = azimuth(x_B, y_B, x_A, y_A)
BC_a = azimuth(x_B, y_B, x_C, y_C)
BP_a = azimuth(x_B, y_B, x_P_0, y_P_0)

CB_a = azimuth(x_C, y_C, x_B, y_B)
CP_a = azimuth(x_C, y_C, x_P_0, y_P_0)

PA_a = azimuth(x_P_0, y_P_0, x_A, y_A)
PB_a = azimuth(x_P_0, y_P_0, x_B, y_B)
PC_a = azimuth(x_P_0, y_P_0, x_C, y_C)

AP_s = distance(x_A, y_A, x_P_0, y_P_0)
BP_s = distance(x_B, y_B, x_P_0, y_P_0)
CP_s = distance(x_C, y_C, x_P_0, y_P_0)
PA_s = distance(x_P_0, y_P_0, x_A, y_A)
PB_s = distance(x_P_0, y_P_0, x_B, y_B)
PC_s = distance(x_P_0, y_P_0, x_C, y_C)

L = np.zeros((len(observations), 1))
A = np.zeros((len(observations), 2))
x_delta = np.zeros((3, 1))
y_delta = np.zeros((3, 1))
rho_bis = 3600 * 180 / np.pi
corr = rho_bis / 100


L[0] = (observations[0] - (AB_a - AP_a)) * 3600
L[1] = (observations[1] - ((BP_a + 360) - BA_a)) * 3600
L[2] = (observations[2] - (BC_a - BP_a)) * 3600
L[3] = (observations[3] - (CP_a - CB_a)) * 3600
L[4] = (observations[4] - CP_s) * 100

A[0][0] = corr * (np.sin(np.radians(AP_a)) / AP_s)
A[0][1] = corr * (-np.cos(np.radians(AP_a)) / AP_s)

A[1][0] = corr * (-np.sin(np.radians(BP_a)) / BP_s)
A[1][1] = corr * (np.cos(np.radians(BP_a)) / BP_s)

A[2][0] = corr * (np.sin(np.radians(BP_a)) / BP_s)
A[2][1] = corr * (-np.cos(np.radians(BP_a)) / BP_s)

A[3][0] = corr * (-np.sin(np.radians(CP_a)) / CP_s)
A[3][1] = corr * (np.cos(np.radians(CP_a)) / CP_s)

A[4][0] = (x_P_0 - x_C) / CP_s
A[4][1] = (y_P_0 - y_C) / CP_s

print("L: ", L)
print("A: ", A)

n = len(observations)
print("\n", "n: ", n)
m = len(A[0])
print("\n", "m: ", m)

L_prim = observations
print("\n", "L_prim: ", "\n", L_prim)
ATPA = np.dot(np.dot(A.T, P), A)
print("\n", "ATPA: ", "\n", ATPA)
ATPA_inv = np.linalg.inv(ATPA)
print("\n", "ATPA_inv: ", "\n", ATPA_inv)
ATPA_inv_ATP = np.dot(np.dot(ATPA_inv, A.T), P)
print("\n", "ATPA_inv_ATP: ", "\n", ATPA_inv_ATP)
dX_hat = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(A.T, P), A)), A.T), P), L)
print("\n", "dX_hat: ", "\n", dX_hat, "(cm)", "\n")
X_hat = np.array([[x_P_0], [y_P_0]]) + dX_hat / 100
print("\n", "X_hat: ", "\n", X_hat)

epsilon = L - np.dot(A, dX_hat)
print("\n", "epsilon: ", "\n", epsilon, "\n")

L_hat = np.dot(A, dX_hat)
print("\n", "L_hat: ", "\n", L_hat)

L_hat_prim = L_prim - epsilon
print("\n", "L_hat_prim: ", "\n", L_hat_prim)

variance = (np.dot(np.dot(epsilon.T, P), epsilon)) / (n - m)
print("\n", "variance: ", variance)

Q_xx = np.linalg.inv(np.dot(np.dot(A.T, P), A))
print("\n", "Q_xx: ", "\n", Q_xx)

C_xx = variance * Q_xx
print("\n", "C_xx: ", "\n", C_xx)

variance_x_P = np.sqrt(variance) * np.sqrt(Q_xx[0][0])
variance_y_P = np.sqrt(variance) * np.sqrt(Q_xx[1][1])
variance_P = np.sqrt(variance_x_P**2 + variance_y_P**2)
print("\n", "variance_P: ", variance_P, "(cm)")

Q_LL = np.dot(np.dot(A, Q_xx), A.T)
print("\n", "Q_LL: ", "\n", Q_LL)

C_LL = variance * Q_LL
print("\n", "C_LL: ", "\n", C_LL)

Q_epsilon_epsilon = np.linalg.inv(P) - Q_LL
print("\n", "Q_epsilon_epsilon: ", "\n", Q_epsilon_epsilon)

C_epsilon_epsilon = variance * Q_epsilon_epsilon
print("\n", "C_epsilon_epsilon: ", "\n", C_epsilon_epsilon)
