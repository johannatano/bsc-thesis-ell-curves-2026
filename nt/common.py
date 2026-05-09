def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result


def legendre(a, p):
    if p == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1


def valuation(n, l):
    if n == 0:
        return float("inf")
    v = 0
    while n % l == 0:
        n //= l
        v += 1
    return v


def fmt_factored(n):
    if n == 0:
        return "0"
    if n > 10**20:
        return str("")
    factors = factorize(abs(n))
    if not factors:
        return str(n)
    s = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
    return f"-{s}" if n < 0 else s


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

def divisors(n):
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)
