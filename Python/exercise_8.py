import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from icecream import ic

months = []
co2_values = []

with open("CO2.txt", "r") as file:
    for line in file:
        line = line.strip()
        months.append(int(line.split()[0]))
        co2_values.append(float(line.split()[1]))  # Convert CO2 values to float

# Compute simple moving average over a 12-month interval (M = 12)
M = 12
moving_average = pd.DataFrame(co2_values).rolling(window=12).mean()
ic(moving_average[10:-4])

# Compute centered moving average over a 12-month interval (M = 12)
M = 12
centered_moving_average = (
    pd.DataFrame(co2_values).rolling(window=12, center=True).mean()
)
ic(centered_moving_average[4:-4])

# Compute exponential moving average with alpha = 0.3
exp_moving_average_03 = pd.DataFrame(co2_values).ewm(alpha=0.3).mean()
exp_moving_average_07 = pd.DataFrame(co2_values).ewm(alpha=0.7).mean()
ic(exp_moving_average_03)
ic(exp_moving_average_07)

# Compute exponential moving average with alpha = 0.3
# exp_moving_average_03 = []
# exp_moving_average_07 = []

# for i in range(len(co2_values)):
#     if i == 0:
#         exp_moving_average_03.append(co2_values[i])
#         exp_moving_average_07.append(co2_values[i])
#     else:
#         exp_moving_average_03.append(
#             0.3 * co2_values[i] + 0.7 * exp_moving_average_03[i - 1]
#         )
#         exp_moving_average_07.append(
#             0.7 * co2_values[i] + 0.3 * exp_moving_average_07[i - 1]
#         )

# Compute autocorrelation coefficients
autocorrelation_coefficients = []
for lag in range(len(co2_values)):
    if lag == 0:
        autocorrelation_coefficients.append(1.0)
    else:
        numerator = sum(
            (co2_values[i] - sum(co2_values) / len(co2_values))
            * (co2_values[i - lag] - sum(co2_values) / len(co2_values))
            for i in range(lag, len(co2_values))
        )
        denominator = sum(
            (co2_values[i] - sum(co2_values) / len(co2_values)) ** 2
            for i in range(len(co2_values))
        )
        autocorrelation_coefficients.append(numerator / denominator)
ic(autocorrelation_coefficients)
# Plot all figures in the same figure
plt.figure(figsize=(13, 9))

# Plot original CO2 values
plt.subplot(2, 3, 1)
plt.plot(months, co2_values)
plt.xlabel("Months")
plt.ylabel("CO2 Value")
plt.title("CO2 Values Over Time")

# Plot computed moving average
plt.subplot(2, 3, 2)
plt.plot(months, moving_average)
plt.xlabel("Months")
plt.ylabel("Moving Average")
plt.title("Simple Moving Average (M = 12)")

# Plot computed centered moving average
plt.subplot(2, 3, 3)
plt.plot(months, centered_moving_average)
plt.xlabel("Months")
plt.ylabel("Centered Moving Average")
plt.title("Centered Moving Average (M = 12)")

# Plot computed exponential moving averages
plt.subplot(2, 3, 4)
plt.plot(months, exp_moving_average_03, label="Alpha = 0.3")
plt.plot(months, exp_moving_average_07, label="Alpha = 0.7")
plt.xlabel("Months")
plt.ylabel("Exponential Moving Average")
plt.title("Exponential Moving Averages")
plt.legend()

# Combined plot
plt.subplot(2, 3, 5)
plt.plot(months, co2_values, label="Original CO2 Values")
plt.plot(months, moving_average, label="Moving Average")
plt.plot(months, centered_moving_average, label="Centered Moving Average")
plt.plot(months, exp_moving_average_03, label="Alpha = 0.3")
plt.plot(months, exp_moving_average_07, label="Alpha = 0.7")
plt.xlabel("Months")
plt.ylabel("CO2 Value")
plt.title("Combined plot")
plt.legend()

# Plot autocorrelation coefficients
plt.subplot(2, 3, 6)
plt.stem(range(len(autocorrelation_coefficients)), autocorrelation_coefficients)
plt.xlabel("Lag")
plt.ylabel("Autocorrelation Coefficient")
plt.title("Autocorrelation Coefficients")

plt.tight_layout()
plt.show()
