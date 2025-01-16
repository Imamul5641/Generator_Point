import sys
import sympy

class GeneratorPoints:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.generator_points = []
        if sympy.isprime(p):
            self.p = p
        else:
            print(p, "is not a prime number")
            self.p = sympy.nextprime(p)
            print("The next prime number is", self.p)

    def modular_sqrt(self, n, p):
        """ Return the square root of n modulo p, if exists, else return None. """
        if pow(n, (p - 1) // 2, p) != 1:
            return None
        elif p % 4 == 3:
            return self.euler_formula(n, p)
        elif p % 4 == 1:
            return self.tonelli_shanks(n, p)
        return None

    def euler_formula(self, n, p):
        return pow(n, (p + 1) // 4, p)

    def tonelli_shanks(self, n, p):
        """Use Tonelli–Shanks algorithm for p ≡ 1 (mod 4)."""
        if pow(n, (p - 1) // 2, p) != 1:
            return None
        q = p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1
        m = s
        c = pow(z, q, p)
        t = pow(n, q, p)
        r = pow(n, (q + 1) // 2, p)
        while t != 1:
            t2i = t
            for i in range(1, m):
                t2i = pow(t2i, 2, p)
                if t2i == 1:
                    break
            b = pow(c, 2 ** (m - i - 1), p)
            m = i
            c = pow(b, 2, p)
            t = (t * c) % p
            r = (r * b) % p
        return r

    def mod_inverse(self, a, p):
        """Calculate modular inverse of a under modulo p."""
        return pow(a, p - 2, p)

    def elliptic_curve_points(self):
        points = []
        for x in range(self.p):
            y2 = (x**3 + self.a * x + self.b) % self.p
            y = self.modular_sqrt(y2, self.p)
            if y is not None:
                points.append((x, y))
                if y != self.p - y:
                    points.append((x, self.p - y))
            # Check for points where y = 0
            if y2 == 0:
                points.append((x, 0))
        self.curve_order = len(points) + 1  # +1 for the point at infinity
        return points

    def add_point(self, xp, yp, xq, yq):
        if xp == 0 and yp == 0:
            return (xq, yq)
        if xq == 0 and yq == 0:
            return (xp, yp)
        if xp == xq and yp == (-yq % self.p):
            return (0, 0)

        if xp == xq and yp == yq:
            # Doubling case
            l = (3 * xp * xp + self.a) * self.mod_inverse(2 * yp, self.p) % self.p
        else:
            # Addition case
            l = (yq - yp) * self.mod_inverse(xq - xp, self.p) % self.p

        xr = (l * l - xp - xq) % self.p
        yr = (l * (xp - xr) - yp) % self.p
        return (xr, yr)

    def point_order(self, xp, yp):
        xq, yq = xp, yp
        xr, yr = self.add_point(xq, yq, xq, yq)  # First doubling
        order = 2

        while (xr != xp or yr != yp):
            xq, yq = xr, yr
            xr, yr = self.add_point(xp, yp, xq, yq)
            order += 1

            # Safety check to avoid infinite loops
            if xr == 0 and yr == 0:
                return order

        return order

    def find_points(self):
        points = self.elliptic_curve_points()
        for point in points:
            x, y = point
            if self.point_order(x, y) == self.curve_order:
                self.generator_points.append((x, y))
        print("Generator points:", self.generator_points)


def main():
    if len(sys.argv) != 4:
        print("Usage: python elliptic_curve.py <a> <b> <p>")
        sys.exit(1)

    try:
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        p = int(sys.argv[3])
    except ValueError:
        print("Error: a, b, and p must be integers.")
        sys.exit(1)

    gp = GeneratorPoints(a, b, p)
    gp.find_points()


if __name__ == "__main__":
    main()

