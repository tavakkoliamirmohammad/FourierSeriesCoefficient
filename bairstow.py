from Polynomial import Polynomial
import math


class Bairstow:
    def __init__(self, r: float, s: float, poly: Polynomial):
        self.r = r
        self.s = s
        self.poly = poly
        self.b = []
        self.c = []

    def calculate_b(self, k: int):
        n = len(self.poly)
        if k == n + 1 or k == n + 2:
            return 0
        return self.poly[k] + self.r * self.calculate_b(k + 1) + self.s * self.calculate_b(k + 2)

    def calculate_c(self, k: int):
        n = len(self.poly)
        if k == n + 1 or k == n + 2:
            return 0
        return self.calculate_b(k) + self.r * self.calculate_c(k + 1) + self.s * self.calculate_c(k + 2)

    def calculate_d(self):
        return self.c[0] * self.c[2] - math.pow(self.c[1], 2)

    def calculate_d1(self):
        return -self.b[0] * self.c[2] + self.c[1] * self.b[1]

    def calculate_d2(self):
        return -self.c[0] * self.b[1] + self.c[1] * self.b[0]

    def do_iteration(self):
        self.update_values()
        d = self.calculate_d()
        d1 = self.calculate_d1()
        d2 = self.calculate_d2()

        self.r = self.r + d1 / d
        self.s = self.s + d2 / d

    def update_values(self):
        self.b = []
        self.c = []
        self.b.append(self.calculate_b(0))
        self.b.append(self.calculate_b(1))
        self.c.append(self.calculate_c(1))
        self.c.append(self.calculate_c(2))
        self.c.append(self.calculate_c(3))
