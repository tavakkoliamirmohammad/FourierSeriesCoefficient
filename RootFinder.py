from Polynomial import Polynomial
from bairstow import Bairstow


class RootFinder:
    def __init__(self, epsilon: float, poly: Polynomial, r: float, s: float):
        self.epsilon = epsilon
        self.poly = poly
        self.r = r
        self.s = s
        self.roots = []

    def do_iteration(self):
        bairstow = Bairstow(self.r, self.s, self.poly)
        while True:
            r = bairstow.r
            s = bairstow.s
            bairstow.do_iteration()
            if abs(bairstow.r - r) < self.epsilon or abs(bairstow.s - s) < self.epsilon:
                break
        poly_roots = Polynomial.find_quadratic_root(1, -bairstow.r, -bairstow.s)
        self.roots += poly_roots
        result_poly = Polynomial([1, -bairstow.r, -bairstow.s])
        self.poly = Polynomial.extended_synthetic_division(self.poly, result_poly)

    def is_done(self):
        if len(self.poly) == 0:
            return False
        return True

    def find(self):
        while True:
            if not self.is_done():
                break
            if len(self.poly) < 2:
                self.roots.append(complex(-self.poly[0]/self.poly[1]))
                break
            self.do_iteration()

        return self.roots
