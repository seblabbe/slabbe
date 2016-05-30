
def random_composition(n, length):
    r"""
    EXEMPLES::

        sage: random_composition(4,2)
        [1, 3]
        sage: random_composition(4,2)
        [2, 2]
        sage: random_composition(4,2)
        [3, 1]
        sage: random_composition(24,9)
        [1, 4, 2, 2, 3, 5, 4, 2, 1]

    Because this is very slow!!::

        sage: C = Compositions(24, length=9)
        sage: %time C.random_element()
        CPU times: user 43.3 s, sys: 31.9 ms, total: 43.3 s
        Wall time: 43.2 s
        [2, 2, 5, 2, 8, 1, 1, 2, 1]
    """
    L = [randint(0,n-length) for _ in range(length-1)]
    L.sort()
    L = [a+i+1 for i,a in enumerate(L)]
    L.insert(0,0)
    L.append(n)
    return [L[i+1]-L[i] for i in range(length)]
