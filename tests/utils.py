def isclose(a, b, rel_tol=1e-07, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
