from exercise_methods import *

l1 = 60 + 0 / 60 + 5 / 3600
l2 = 60 + 0 / 60 + 5 / 3600
l3 = 59 + 59 / 60 + 58 / 3600
l4 = 100.008
l5 = 99.997

observations = np.array([[l1], [l2], [l3], [l4], [l5]])
weights = np.array([[0.25], [0.25], [0.25], [1], [1]])
P = np.zeros((len(weights), len(weights)))
for i in range(len(weights)):
    P[i][i] = weights[i][0]
ic(P)

x_A, y_A = 6500000, 1500000
x_B, y_B = 6500060, 1500080


azimuth_AB = azimuth(x_A, y_A, x_B, y_B)
azimuth_AP = azimuth_AB - l1

x_p_0 = x_A + l4 * np.cos(np.radians(azimuth_AP))
y_p_0 = y_A + l4 * np.sin(np.radians(azimuth_AP))

APB = [x_A, y_A, x_p_0, y_p_0, x_B, y_B]
BAP = [x_B, y_B, x_A, y_A, x_p_0, y_p_0]
PBA = [x_p_0, y_p_0, x_B, y_B, x_A, y_A]
AP = [x_A, y_A, x_p_0, y_p_0]
PB = [x_p_0, y_p_0, x_B, y_B]

angles = [[APB, "j"], [BAP, "k"], [PBA, "i"]]
distances = [[AP, "j"], [PB, "i"]]

L, A = linj_figur(observations, angles, distances)
ic(L)
ic(A)
