from math import factorial
import scipy.signal
import numpy as np

from Polynomial import Polynomial


def combination(m: int, n: int):
    return factorial(m) / (factorial(n) * factorial(m - n))


# taylor series 1 / (1 - a * z)^ n
def taylor_series_kth_term(k: int, a, n: int):
    if k < 0:
        return 0
    return k, a ** k * combination(k + n - 1, k)


def taylor_series_k_first_term(k: int, a, n: int):
    series = []
    for i in range(k):
        series.append(taylor_series_kth_term(i, a, n))
    return series


# laurent series 1 / (1 - a * z) ^ n
def laurent_series_kth_term(k: int, a, n: int):
    return -k - n, taylor_series_kth_term(k, a, n)[1]


def laurent_series_k_first_term(k: int, a, n: int):
    series = []
    for i in range(k):
        series.append(laurent_series_kth_term(i, a, n))
    return series


def extract_poles(numerator: Polynomial, denominator: Polynomial, z0: complex):
    coefficient, poles, poly = scipy.signal.residue(numerator.coefficient, denominator.coefficient)
    poles_terms = {}
    for (c, p) in zip(coefficient, poles):
        if p in poles_terms:
            poles_terms[p - z0].append(c)
        else:
            poles_terms[p - z0] = [c]

    return poles_terms, poly


# print the first m term series around point z
def pole_series(pole: complex, coefficient: list, m: int, z: complex):
    series = {}
    if abs(pole) > abs(z) and pole != 0:
        # taylor series
        for power, c in enumerate(coefficient):
            power += 1
            taylor = taylor_series_k_first_term(m, 1 / pole, power)
            for p, val in taylor:
                if p in series:
                    series[p] += (-pole) ** (-power) * c * val
                else:
                    series[p] = (-pole) ** (-power) * c * val

    elif abs(pole) < abs(z) and pole != 0:
        # laurent series
        for power, c in enumerate(coefficient):
            power += 1
            laurent = laurent_series_k_first_term(m, pole, power)
            for p, val in laurent:
                if p in series:
                    series[p] += (c * val)
                else:
                    series[p] = c * val
    else:
        for i, c in enumerate(coefficient):
            series[-(i + 1)] = c
    required_key = sorted(series.keys(), reverse=True)[:m]
    series = {k: v for k, v in series.items() if k in required_key}
    return series


def get_series_in_convergence_area(pole_term: dict, m: int, z: complex) -> dict:
    series = {}
    for pole in pole_term:
        fractional_series = pole_series(pole, pole_term[pole], m, z)
        for k in fractional_series:
            if k in series:
                series[k] += fractional_series[k]
            else:
                series[k] = fractional_series[k]
    return series


# TODO handle 1/z
# TODO handle if numerator is bigger

def get_series_in_all_convergence_area(numerator: Polynomial, denominator: Polynomial, m: int, z0: complex):
    pole_terms, poly = extract_poles(numerator, denominator, z0)
    poles = []
    for p in pole_terms.keys():
        if abs(p) not in poles:
            poles.append(abs(p))

    if 0 not in pole_terms.keys():
        poles = [0] + poles
    poles = sorted(poles)
    for i in range(1, len(poles) + 1):
        if i == len(poles):
            z = poles[i - 1] + 1
            if poles[i - 1] != 0:
                print("Series in the convergence area", round(poles[i - 1], 2), "< |z| ")
            else:
                print("Series in the convergence area 0 < |z|")

        else:
            z = (poles[i - 1] + poles[i]) / 2
            print("Series in the convergence area ", round(poles[i - 1], 2), " < |z| < ", round(poles[i], 2))
        series = get_series_in_convergence_area(pole_terms, m, z)
        poly = poly[::-1]
        for deg in range(len(poly)):
            if deg in series:
                series[deg] += poly[deg]
            else:
                series[deg] = poly[deg]
        show_series(series, 3, z0)


def show_series(poly: dict, precision: int, z0: complex):
    for pow, coeff in sorted(poly.items()):
        if coeff != 0:
            x, y = round(coeff.real, precision), round(coeff.imag, precision)
            if y == 0:
                print(f'{x:+}', end="")
            else:
                print(f'{x:+}{y:+}i', end="")
            if pow != 0:
                print("*(x{0:+}{1:+}i)**{2} ".format(round(z0.real, precision), round(z0.imag, precision), int(pow)),
                      end="")
    print()


get_series_in_all_convergence_area(Polynomial([1]), Polynomial([1, 6, 12, 8, 0]),10, 0)
