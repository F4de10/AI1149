import numpy as np
from icecream import ic
import scipy.stats


def adjustment_by_elements(observations, A, P, C, unit="none", print="true"):
    """
    Perform adjustment by elements.

    Args:
        observations (2D array): The observations.
        A (1D array): The design matrix.
        P (2D array): The weight matrix.
        C (2D array): The constant matrix.
        unit (str): The unit of measurement. Either "degree" or "meter".
        print (str): Whether to enable printing of intermediate results. Either "true" or "false".

    Returns:
        tuple: A tuple containing the adjusted observations (X) and the covariance matrix of the residuals (Q_epsilon_epsilon).
    """

    if unit == "degree":
        correction = 3600
    elif unit == "meter":
        correction = 1000
    elif unit == "centimeter":
        correction = 100
    else:
        correction = 1

    if print == "false":
        ic.disable()
    else:
        ic.enable()

    n = len(observations)
    ic(n)
    m = len(A[0])
    ic(m)

    L_prim = observations
    L = L_prim - C
    ic(L)
    ATPA = np.dot(np.dot(A.T, P), A)
    ic(ATPA)
    ATPA_inv = np.linalg.inv(ATPA)
    ic(ATPA_inv)
    ATPA_inv_ATP = np.dot(np.dot(ATPA_inv, A.T), P)
    ic(ATPA_inv_ATP)
    dX = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(A.T, P), A)), A.T), P), L)
    ic(dX)

    L_hat = np.dot(A, dX)
    ic(L_hat)

    if unit == "custom":
        C = C / [[3600], [3600], [3600], [3600], [100]]
        L_hat_prim = (L_hat / [[3600], [3600], [3600], [3600], [100]]) + C
    else:
        L_hat_prim = L_hat / correction + C
    ic(L_hat_prim)

    epsilon = correction * (L - np.dot(A, dX))
    ic(epsilon)

    variance = (np.dot(np.dot(epsilon.T, P), epsilon)) / (n - m)
    ic(variance)
    unit_weight_standard_error = np.sqrt(variance)
    ic(unit_weight_standard_error)

    Q_xx = np.linalg.inv(np.dot(np.dot(A.T, P), A))
    ic(Q_xx)
    C_xx = variance * Q_xx
    ic(C_xx)

    """variance_x_P = np.sqrt(variance) * np.sqrt(Q_xx[0][0])
    variance_y_P = np.sqrt(variance) * np.sqrt(Q_xx[1][1])
    variance_P = np.sqrt(variance_x_P**2 + variance_y_P**2)
    ic(variance_P)"""

    Q_LL = np.dot(np.dot(A, Q_xx), A.T)
    ic(Q_LL)
    C_LL = variance * Q_LL
    ic(C_LL)

    Q_epsilon_epsilon = np.linalg.inv(P) - Q_LL
    ic(Q_epsilon_epsilon)
    C_epsilon_epsilon = variance * Q_epsilon_epsilon
    ic(C_epsilon_epsilon)

    ic.enable()

    return dX, Q_epsilon_epsilon


def distance(x1, y1, x2, y2):
    """
    Calculate the distance between two points.

    Args:
        x1 (float): The x-coordinate of the first point.
        y1 (float): The y-coordinate of the first point.
        x2 (float): The x-coordinate of the second point.
        y2 (float): The y-coordinate of the second point.

    Returns:
        float: The distance between the two points.
    """
    distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return distance


def azimuth(x1, y1, x2, y2):
    """
    Calculate the azimuth between two points.

    Args:
        x1 (float): The x-coordinate of the first point.
        y1 (float): The y-coordinate of the first point.
        x2 (float): The x-coordinate of the second point.
        y2 (float): The y-coordinate of the second point.

    Returns:
        float: The azimuth between the two points.
    """
    azimuth = np.arctan2(y2 - y1, x2 - x1) * 180 / np.pi
    if azimuth < 0:
        azimuth += 360
    return azimuth


def linj_vinkel(angle, point):
    """
    Calculate the line angle and distance based on given parameters.

    Args:
        angle (list): List of angle coordinates [x_i, y_i, x_j, y_j, x_k, y_k].
        point (str): The point to calculate the line angle and distance for. Can be "i", "j", "k", or any other value.

    Returns:
        tuple: A tuple containing the line angle and distance based on the given point.
            If point is "i", returns (C, a, b).
            If point is "j", returns (C, c, d).
            If point is "k", returns (C, e, f).
            If point is any other value, returns (C, a, b, c, d, e, f).

    Raises:
        ValueError: If the point parameter is not one of the valid options ("i", "j", "k").

    """
    [x_i, y_i, x_j, y_j, x_k, y_k] = angle
    alpha_ik = azimuth(x_i, y_i, x_k, y_k)
    alpha_ij = azimuth(x_i, y_i, x_j, y_j)

    beta_ijk = alpha_ik - alpha_ij
    if beta_ijk < 0:
        beta_ijk += 360

    s_ij = distance(x_i, y_i, x_j, y_j)
    s_ik = distance(x_i, y_i, x_k, y_k)
    a = np.sin(np.radians(alpha_ik)) / s_ik - np.sin(np.radians(alpha_ij)) / s_ij
    b = -np.cos(np.radians(alpha_ik)) / s_ik + np.cos(np.radians(alpha_ij)) / s_ij
    c = np.sin(np.radians(alpha_ij)) / s_ij
    d = -np.cos(np.radians(alpha_ij)) / s_ij
    e = -np.sin(np.radians(alpha_ik)) / s_ik
    f = np.cos(np.radians(alpha_ik)) / s_ik

    C = beta_ijk
    if point == "i":
        return C, a, b
    elif point == "j":
        return C, c, d
    elif point == "k":
        return C, e, f
    else:
        raise ValueError(
            "Invalid point parameter. Must be one of the valid options: 'i', 'j', 'k'."
        )


def linj_avstånd(distance_points, point):
    """
    Calculate the distance between a line segment and a point.

    Parameters:

    distance_points (list): List containing the coordinates of the line segment's endpoints [x_i, y_i, x_j, y_j].
    point (str): Indicates the point to calculate the distance from. Can be "i", "j", or any other value.

    Returns:
    tuple: Tuple containing the distance between the line segment and the point, as well as the coefficients of the line equation.

    """
    [x_i, y_i, x_j, y_j] = distance_points
    s_ij = distance(x_i, y_i, x_j, y_j)
    a = -(x_j - x_i) / s_ij
    b = -(y_j - y_i) / s_ij
    d = -a
    e = -b
    C = s_ij
    if point == "i":
        return C, a, b
    elif point == "j":
        return C, d, e
    else:
        return C, a, b, d, e


def linj_figur(angles, distances, correction_angles=3600, correction_distances=1000):
    """
    Calculate the linear figure based on observations, angles, and distances.

    Args:
        observations (list): List of observation points.
        angles (list): List of angle measurements.
        distances (list): List of distance measurements.
        correction_angles (float, optional): Correction factor for angles. Defaults to 3600 (seconds).
        correction_distances (float, optional): Correction factor for distances. Defaults to 1000 (mm).

    Returns:
        tuple: A tuple containing the calculated linear figure L and the adjustment matrix A.
    """
    A = np.zeros((len(angles + distances), 2))
    C = np.zeros((len(angles + distances), 1))

    rho = correction_angles * 180 / np.pi
    corr = rho / correction_distances

    for i in range(len(angles)):
        A[i][0] = corr * linj_vinkel(angles[i][0], angles[i][1])[1]
        A[i][1] = corr * linj_vinkel(angles[i][0], angles[i][1])[2]
        C[i][0] = linj_vinkel(angles[i][0], angles[i][1])[0]

    for i in range(len(angles), len(angles) + len(distances)):
        j = i - len(angles)
        A[i][0] = linj_avstånd(distances[j][0], distances[j][1])[1]
        A[i][1] = linj_avstånd(distances[j][0], distances[j][1])[2]
        C[i][0] = linj_avstånd(distances[j][0], distances[j][1])[0]
    return A, C


def calculate_redundancy(A, P, Q_ee, print_check="true"):
    """
    Calculate the redundancy matrix and diagonal redundancy vector.

    Parameters:
    P (numpy.array): The matrix P.
    Q_ee (numpy.array): The matrix Q_ee.

    Returns:
    tuple: A tuple containing the redundancy matrix R and the diagonal redundancy vector r.
    """
    R = Q_ee * P
    r = np.diag(R)

    if print_check == "false":
        ic.disable()
    else:
        ic.enable()
    ic(R)
    ic(r)

    n = len(A)
    m = len(A[0])
    ic(n, m)

    ic(n - m)
    ic(round(sum(r), 9))
    ic.enable()

    if print_check == "true":
        if n - m == round(sum(r), 9):
            print("The redundancy matrix is correct. n-m = sum(r)")
        else:
            print("The redundancy matrix is incorrect. n-m != sum(r)")
    else:
        pass

    return R, r


def chi_square_test(epsilon, P, variance_factor):
    omega = np.dot(np.dot(epsilon.T, P), epsilon)
    ic(omega)
    chi = omega / variance_factor
    ic(chi)
    p_value = scipy.stats.chi2.cdf(chi, len(epsilon))
    return omega, p_value


if __name__ == "__main__":
    a = 5
