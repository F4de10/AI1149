from exercise_methods import *

print("problem 3.1", "\n")
observations = np.array([[60 + 3 / 3600], [60 + 3 / 3600], [59 + 59 / 60 + 51 / 3600]])
P = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
A = np.array([[1, 0], [0, 1], [-1, -1]])
C = np.array([[0], [0], [180]])

adjustment_by_elements(observations, A, P, C, "degree")

print("\n", "problem 3.2", "\n")
observations = np.array([[1.680], [-0.148], [0.557]])
P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
A = np.array([[1], [1], [1]])
C = np.array([[-10.421], [-12.261], [-11.559]])
adjustment_by_elements(observations, A, P, C, "degree")

print("\n", "problem 3.3", "\n")
observations = np.array([[100.006], [99.999], [200.002]])

# Initialize P as a diagonal matrix
P = np.zeros((len(observations), len(observations)))

# Populate the diagonal elements
for i in range(len(observations)):
    P[i][i] = np.round(40000 / (observations[i][0] ** 2))
ic(P)

A = np.array([[1, 0], [0, 1], [1, 1]])
C = np.array([[0], [0], [0]])
adjustment_by_elements(observations, A, P, C, "meter")

print("\n", "problem 3.4", "\n")
observations = np.array([[1.054], [-0.7], [-0.351]])
P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
A = np.array([[-1, 1], [0, -1], [1, 0]])
C = np.array([[0], [1], [-1]])
X = adjustment_by_elements(observations, A, P, C, "meter")
ic(X)

print("\n", "problem 3.5", "\n")
observations = np.array([[1.001], [1.000], [1.004], [-1.996]])
P = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
A = np.array([[1, 0], [-1, 0], [0, 1], [1, -1]])
C = np.array([[0], [2], [-2], [0]])
adjustment_by_elements(observations, A, P, C, "meter")
