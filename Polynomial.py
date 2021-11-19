import math


class Polynomial:
    def __init__(self, coefficient: list):
        self.coefficient = coefficient

    def __len__(self):
        return len(self.coefficient) - 1

    def __getitem__(self, item):
        return self.coefficient[len(self) - item]

    @staticmethod
    def find_quadratic_root(a: float, b: float, c: float):
        delta = b ** 2 - 4 * a * c
        roots = []
        if delta < 0:
            r1 = complex(-b / (2 * a), math.sqrt(abs(delta)) / (2 * a))
            r2 = complex(-b / (2 * a), -math.sqrt(abs(delta)) / (2 * a))
            roots.append(r1)
            roots.append(r2)
        elif delta == 0:
            r = complex(-b / (2 * a))
            roots.append(r)
        else:
            r1 = complex(-b / (2 * a) + math.sqrt(abs(delta)) / (2 * a))
            r2 = complex(-b / (2 * a) - math.sqrt(abs(delta)) / (2 * a))
            roots.append(r1)
            roots.append(r2)
        return roots

    
    @staticmethod
    def extended_synthetic_division(poly_dividend, poly_divisor):
        dividend = poly_dividend.coefficient
        divisor = poly_divisor.coefficient
        out = list(dividend)  # Copy the dividend
        normalizer = divisor[0]
        for i in range(len(dividend) - (len(divisor) - 1)):
            out[i] /= normalizer
            coef = out[i]
            if coef != 0:
                for j in range(1, len(divisor)):
                    out[i + j] += -divisor[j] * coef
        separator = -(len(divisor) - 1)
        # Returns a Quotient and a Remainder
        return Polynomial(out[:separator]), Polynomial(out[separator:])
    
    @staticmethod
    def multiply(poly1, poly2):
        coefs1 = poly1.coefficient
        coefs2 = poly2.coefficient
        product = [0] * (len(coefs1) + len(coefs2) - 1)
        for i in range(len(coefs1)):
            for j in range(len(coefs2)):
                product[i + j] += coefs1[i] * coefs2[j]
        result = Polynomial(product)
        return result 

    # Evaluating a polynomial using the Horner's method
    @staticmethod
    def evaluate(poly, x):
        coefs = poly.coefficient
        result = coefs[0]
        for i in range(1, len(coefs)):
            result = result*x + coefs[i]
        return result
    
    @staticmethod
    def add(poly1, poly2):
        coefs1 = poly1.coefficient
        coefs2 = poly2.coefficient
        if (len(coefs1) < len(coefs2)): # If the degree of the second poly is more, then choose it to intialize the result array
            tmp = coefs1
            coefs1 = coefs2
            coefs2 = tmp
        coefs1.reverse() # The coefficient lists are reversed, so that they are in an acsending order,
        coefs2.reverse() # For the algorithm to work right.
        result = coefs1
        for i in range(0, len(coefs2)): 
            result[i] += coefs2[i]
        result.reverse() # Change the order of coefficients to descending order again
        polyResult = Polynomial(result)
        return polyResult
    
    @staticmethod
    def sub(poly1, poly2):
        operand1IsLess = False # Variable to see if the degree of the first operand is less, so that we would negate the result gotten from subtracting.
        coefs1 = poly1.coefficient
        coefs2 = poly2.coefficient
        if (len(coefs1) < len(coefs2)): # If the first operand has a less degree, swap the coefficients, so that the subtraction is done with the first operand
            tmp = coefs1                # Having more degree than the second one
            coefs1 = coefs2
            coefs2 = tmp
            operand1IsLess = True
        coefs1.reverse()
        coefs2.reverse()
        result = coefs1
        for i in range(0, len(coefs2)):
            result[i] -= coefs2[i]
        result.reverse()
        if (operand1IsLess): # Multiply the result by -1 if the degree of the first operand is less than the second one
            result = [-x for x in result]
        polyResult = Polynomial(result)
        return polyResult