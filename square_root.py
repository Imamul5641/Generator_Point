import argparse
import matplotlib.pyplot as plt
import sympy

class Square_root:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
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
            return None  # n is not a quadratic residue modulo p
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

    def plot_elliptic_curve(self, points):
        x_vals = [point[0] for point in points]
        y_vals = [point[1] for point in points]

        plt.figure(figsize=(8, 8))
        plt.scatter(x_vals, y_vals, color='blue')
        plt.title(f'Elliptic Curve: y^2 = x^3 + {self.a}x + {self.b} (mod {self.p})')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)
        plt.show()

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
        return points


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Elliptic Curve Point Calculator")
    parser.add_argument("a", type=int, help="Coefficient 'a' in the elliptic curve equation")
    parser.add_argument("b", type=int, help="Coefficient 'b' in the elliptic curve equation")
    parser.add_argument("p", type=int, help="Prime modulus 'p'")

    args = parser.parse_args()

    sr = Square_root(args.a, args.b, args.p)
    points = sr.elliptic_curve_points()
    print("Points on the elliptic curve:")
    for point in points:
        print(point)
    print("Number of points including origin o(0,0):", len(points) + 1)
    sr.plot_elliptic_curve(points)

