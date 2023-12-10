from matplotlib import pyplot as plt
import seaborn as sns
from exercise_methods import *
from scipy.stats import chi2

print("problem 3.6", "\n")
print("a)", "\n")

x = np.array([0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 3.8]).reshape((-1, 1))
print(f"x: {x}")
y = np.array([3.40, 3.60, 3.98, 4.22, 4.61, 5.12, 5.40, 6.21, 6.48, 7.10])
print(f"y: {y}")
n = len(x)

sns.regplot(x=x, y=y, ci=False, line_kws={"color": "red"})
plt.show()

"""model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print_statement(f"coefficient of determination: {r_sq}")
print_statement(f"𝛼: {model.intercept_}")
print_statement(f"β: {model.coef}")"""

print("b)", "\n")
q_11 = n
print(f"q_11: {q_11}")
q_12 = sum(x)
print(f"q_12: {q_12}")
q_21 = q_12
print(f"q_21: {q_21}")
q_22 = sum(x**2)
print(f"q_22: {q_22}")
w_1 = sum(y)
print(f"w_1: {w_1}")
w_2 = np.dot(y, x)
print(f"w_2: {w_2}")

f = 1 / ((q_11 * q_22) - (q_12**2))
a = f * (q_22 * w_1 - q_12 * w_2)
b = f * (-q_12 * w_1 + q_11 * w_2)

print("𝛼:", a, "'\", "'β:"', b, "\n")

print("c)", "\n")

df = n - 2
print(f"degrees of freedom: {df}")
print("\n")
epsilon_i = sum((y - (a + (b * x)).reshape((1, -1))))
print(f"epsilon_i: {epsilon_i}")
print("\n")
sigma = np.sqrt((1 / df) * sum(sum((y - (a + (b * x)).reshape((1, -1))) ** 2)))
print(f"sigma_0: ±{sigma}")

sigma_a = sigma * np.sqrt(q_22 / ((q_11 * q_22) - (q_12**2)))
print(f"sigma_a: {sigma_a}")
sigma_b = sigma * np.sqrt(q_11 / ((q_11 * q_22) - (q_12**2)))
print(f"sigma_b: {sigma_b}")

print("\n", "d)", "\n")

known_error = 0.15
chi2_value = df * ((sigma**2) / (known_error**2))
print(f"chi2: {chi2_value}")

# chi2_0.05 (8) = 15.507
# 5.68 < 15.5 --> do not reject null hypothesis
