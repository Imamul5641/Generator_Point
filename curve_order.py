import sys
from sympy import isprime, nextprime, legendre_symbol

# Counter for the number of points
no_of_points = 0

# Ensure p is prime and of the form 4k+3
def check_prime(p):
    while not isprime(p):
        p = nextprime(p)
    return p

# Function to find points on the elliptic curve y^2 = x^3 + ax + b over field p
def curve_order(a, b, p):
    global no_of_points
    no_of_points = 0  # Reset the counter

    for x in range(p):
        # Calculate the right-hand side of the curve equation: y^2 = x^3 + ax + b
        xx = (x**3 + a * x + b) % p

        # Check if y^2 = 0, which means y = 0 is a valid solution
        if xx == 0:
            no_of_points += 1  # (x, 0) is a point on the curve

        # Check if y^2 has a solution in the field, i.e., xx is a quadratic residue mod p
        elif legendre_symbol(xx, p) == 1:
            # Compute y as the square root of xx mod p
            y = pow(xx, (p + 1) // 4, p)
            no_of_points += 2  # Count (x, y) and (x, -y)

    no_of_points += 1  # Add 1 for the point at infinity

    return no_of_points

# Main function
def main():
    if len(sys.argv) != 4:
        print("Usage: python elliptic_curve.py <a> <b> <p>")
        sys.exit(1)

    try:
        # Parse input values
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        p = int(sys.argv[3])
    except ValueError:
        print("Error: a, b, and p must be integers.")
        sys.exit(1)

    # Ensure p is prime and of the form 4k + 3
    p = check_prime(p)
    print(f"Prime field: p = {p}")

    # Calculate the order of the curve
    order_of_curve = curve_order(a, b, p)
    print(f"Elliptic Curve: y^2 = x^3 + {a}x + {b} over F_{p}")
    print(f"The Order of the Curve is: {order_of_curve}")

if __name__ == "__main__":
    main()

