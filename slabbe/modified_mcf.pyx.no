
cdef class Modified(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Modified, Brun, Reverse
        sage: algo = Modified(Brun(), {123:1.2})
        sage: algo
        Modified 3-dimensional continued fraction algorithm

    ::

        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()          # not tested
        sage: algo._test_definition()

    For example::

        sage: f = lambda g: Modified(Brun(), {123:g}).lyapunov_exponents(n_iterations=1000000)[-1]
        sage: plot(f, (0.95, 1.05), adaptive_recursion=0)   # not tested
        Graphics object consisting of 1 graphics primitive


    The motivation for this class is to study the constant in front of the not
    unimodular matrix for branch 4 in Reverse algorithm::

        sage: algo = Modified(Reverse(), {4:1})
        sage: algo.lyapunov_exponents(n_iterations=1000000)    # abs tol 0.1
        (0.407428273374779, -0.10394067692044444, 1.2551140500375464)
        sage: algo = Modified(Reverse(), {4:2.^(1/3)})       # unimodular case
        sage: algo.lyapunov_exponents(n_iterations=1000000)    # abs tol 0.1
        (0.3629971890178995, -0.14441058920752092, 1.397828395305838)

    Question: what is the 1-theta2/theta1 de l'algo Reverse?

    More systematically::

        sage: f = lambda g: Modified(Reverse(), {4:g}).lyapunov_exponents(n_iterations=1000000)[-1]
        sage: plot(f, (0.95, 1.05), adaptive_recursion=0)  # not tested
        Graphics object consisting of 1 graphics primitive
    """
    cdef MCFAlgorithm _algo
    cdef dict _gamma
    def __init__(self, int dim, MCFAlgorithm, dict gamma):
        self._algo = MCFAlgorithm
        self._gamma = gamma
        self.dim = dim

    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Modified, Brun
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: algo = Modified(Brun(), {132:1})
            sage: algo(P)
            ((0.3, 0.6, 0.20000000000000007), (0.2, 0.6, 0.3))

        ::

            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: algo = Modified(Brun(), {123:2})
            sage: algo(P)
            ((0.6, 1.2, 0.40000000000000013), (0.1, 0.3, 0.15))
        """
        cdef int branch
        branch = self._algo.call(P)
        if branch in self._gamma:
            g = self._gamma[branch]
            P.x[0] *= g
            P.x[1] *= g
            P.x[2] *= g
            P.a[0] /= g
            P.a[1] /= g
            P.a[2] /= g
        return branch

    def substitutions(self):
        r"""
        EXAMPLES::

        """
        return self._algo.substitutions()

    def dual_substitutions(self):
        r"""
        EXAMPLES::

        """
        return self._algo.dual_substitutions()

