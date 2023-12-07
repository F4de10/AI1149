import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from icecream import ic

with open("RSL.txt", "r") as file:
    lines = file.readlines()
    dates = []
    values = []
    for line in lines:
        data = line.split()
        dates.append(int(data[0]))
        if len(data) < 2:
            break
        values.append(int(data[1]))

# Perform linear regression
regression = LinearRegression()
regression.fit([[date] for date in dates], values)
predicted_values = regression.predict([[date] for date in dates])

# Get the line equation
slope = regression.coef_[0]
intercept = regression.intercept_
line_equation = f"y = {slope:.2f}x + {intercept:.2f}"
landhöjning = f"landhöjningen är {abs(slope):.2f} mm/år"

# Plotting the series and regression line
plt.plot(
    dates,
    values,
    color="blue",
    linestyle="solid",
    linewidth=1,
    marker="o",
    markerfacecolor="blue",
    markersize=3,
    label="Actual",
)
plt.plot(
    dates,
    predicted_values,
    color="red",
    linestyle="dashed",
    linewidth=1,
    label="Predicted",
)
plt.xlabel("Dates")
plt.ylabel("Values")
plt.title("Series Plot with Linear Regression")
plt.legend()

# Add line equation to the plot
plt.text(
    dates[-1],
    values[0],
    line_equation,
    fontsize=12,
    verticalalignment="top",
    horizontalalignment="right",
    color="black",
)

plt.text(
    dates[-1] + 5,
    values[0] - 50,
    landhöjning,
    fontsize=12,
    verticalalignment="top",
    horizontalalignment="right",
    color="black",
)

plt.show()

# Create a residual time series by removing the linear trend
residuals = [value - (slope * date + intercept) for date, value in zip(dates, values)]
plt.plot(dates, residuals)
plt.xlabel("Dates")
plt.ylabel("residuals")
plt.title("Residuals time series")
plt.show()
ic(residuals)

# Perform DFT on the residual time series
N = len(residuals)
n = np.arange(N / 2)
ic(n)
t = np.arange(N)
ic(t)
a = [
    (2 / N) * np.sum(residuals * np.cos((2 * np.pi * (t + 1) * (n + 1)) / N)) for n in n
]
ic(a)
b = [
    (2 / N) * np.sum(residuals * np.sin((2 * np.pi * (t + 1) * (n + 1)) / N)) for n in n
]
ic(b)

# Calculate the frequencies
frequencies = n / N
ic(frequencies)

# Calculate the amplitude
amplitudes = [np.sqrt(a**2 + b**2) for a, b in zip(a, b)]
ic(amplitudes)


# Plot the frequency-amplitude plot
plt.plot(frequencies, amplitudes)
plt.xlabel("Frequency")
plt.ylabel("Amplitude")
plt.title("Frequency-Amplitude Plot")
plt.show()
