import sys
from sympy import mod_inverse

# Define the point addition function
def add_point(a, b, p, xp, yp, xq, yq):
    if xp == 0 and yp == 0:
        return (xq, yq)
    if xq == 0 and yq == 0:
        return (xp, yp)
    if xp == xq and yp == (-yq % p):
        return (0, 0)

    if xp == xq and yp == yq:
        # Doubling case
        l = (3 * xp * xp + a) * mod_inverse(2 * yp, p) % p
    else:
        # Addition case
        l = (yq - yp) * mod_inverse(xq - xp, p) % p

    xr = (l * l - xp - xq) % p
    yr = (l * (xp - xr) - yp) % p
    return (xr, yr)

# Define the order-finding function
def point_order(a, b, p, xp, yp):
    xq, yq = xp, yp
    xr, yr = add_point(a, b, p, xq, yq, xq, yq)  # First doubling
    order = 2

    while (xr != xp or yr != yp):
        xq, yq = xr, yr
        xr, yr = add_point(a, b, p, xp, yp, xq, yq)
        order += 1

        # Safety check to avoid infinite loops
        if xr == 0 and yr == 0:
            return order

    return order

# Main function for command-line execution
def main():
    if len(sys.argv) != 6:
        print("Usage: python elliptic_curve.py <a> <b> <p> <xp> <yp>")
        sys.exit(1)

    try:
        # Parse input values
        a = int(sys.argv[1])
        b = int(sys.argv[2])
        p = int(sys.argv[3])
        xp = int(sys.argv[4])
        yp = int(sys.argv[5])
    except ValueError:
        print("Error: All inputs must be integers.")
        sys.exit(1)

    # Calculate the order of the point
    order = point_order(a, b, p, xp, yp)
    print(f"Order of the point ({xp}, {yp}): {order}")

if __name__ == "__main__":
    main()

