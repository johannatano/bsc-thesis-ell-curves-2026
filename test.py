def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result


def factorize(n):
    # return list(factorint(n).items())  # → [(2, 2), (3, 1)]
    """Return list of (prime, exponent) pairs for n > 1."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors


print(euler_phi(10))  # Should print 4, since the numbers coprime to 12 are 1, 5, 7, 11
