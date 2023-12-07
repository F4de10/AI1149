from exercise_methods import *

print("problem 4", "\n")
l1 = 60 + 0 / 60 + 00.2 / 3600
l2 = 70 + 53 / 60 + 38.2 / 3600
l3 = 49 + 6 / 60 + 25.0 / 3600
l4 = 29 + 59 / 60 + 57.7 / 3600
l5 = 3464.121

observations = np.array([[l1], [l2], [l3], [l4], [l5]])

print("\n", "a)", "\n")
print(f"Independent unknowns: P_x and P_y (2 unknowns)")

print("\n", "b)", "\n")
variance = np.array([[3], [3], [3], [3], [1.5]])
variance_factor = 3

P = np.zeros((len(variance), len(variance)))
for i in range(len(variance)):
    P[i][i] = variance_factor**2 / variance[i][0] ** 2
ic(P)

print("\n", "c)", "\n")
given_coordinates = np.array(
    [[2732.651, 1000.314], [1000.600, 2000.314], [1000.600, 6500.314]]
)

x_A, y_A = given_coordinates[0, 0], given_coordinates[0, 1]
x_B, y_B = given_coordinates[1, 0], given_coordinates[1, 1]
x_C, y_C = given_coordinates[2, 0], given_coordinates[2, 1]

azimuth_CB = azimuth(x_C, y_C, x_B, y_B)
azimuth_CP = azimuth_CB + l4

x_P_0 = x_C + l5 * np.cos(np.radians(azimuth_CP))
y_P_0 = y_C + l5 * np.sin(np.radians(azimuth_CP))
print(f"x_P_0: {x_P_0}, y_P_0: {y_P_0}")

print("\n", "d)", "\n")
APB = [x_A, y_A, x_P_0, y_P_0, x_B, y_B]
BAP = [x_B, y_B, x_A, y_A, x_P_0, y_P_0]
BPC = [x_B, y_B, x_P_0, y_P_0, x_C, y_C]
CBP = [x_C, y_C, x_B, y_B, x_P_0, y_P_0]
CP = [x_C, y_C, x_P_0, y_P_0]

angles = [[APB, "j"], [BAP, "k"], [BPC, "j"], [CBP, "k"]]
distances = [[CP, "j"]]

A, C = linj_figur(angles, distances, 3600, 100)
L_prim = observations
L = L_prim - C
ic(A)
ic(L_prim)
ic(C)
corrections = np.array([[3600], [3600], [3600], [3600], [100]])

L_prim = np.multiply(corrections, L_prim)
C = np.multiply(corrections, C)


print("\n", "e)", "\n")

dX, Q_epsilon = adjustment_by_elements(L_prim, A, P, C, "custom")
ic(dX)
X = x_P_0 + dX[0][0] / 100
Y = y_P_0 + dX[1][0] / 100
ic(X)
ic(Y)


