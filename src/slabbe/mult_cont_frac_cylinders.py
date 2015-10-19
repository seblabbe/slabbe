r"""
Multidimensional Continued Fraction Algorithm's Cylinders

EXAMPLES:

The 1-cylinders of ARP transformation as matrices::

    sage: ARP = TransformationARP()
    sage: list(ARP.n_cylinders_iterator(1))
    [
    [1 0 0]  [1 0 0]  [1 1 0]  [1 1 0]  [1 1 1]  [1 1 1]
    [1 1 0]  [1 1 1]  [1 1 1]  [2 2 1]  [2 1 1]  [2 2 1]
    [3 2 1], [3 2 1], [3 2 1], [3 2 1], [3 2 1], [3 2 1]
    ]

Ces calculs illustrent le bounded distorsion de ratio=4 pour ARP
multiplicatif (2 avril 2014)::

    sage: T = TransformationARP_multi(2)
    sage: T.norm_ratio_max(1)
    3
    sage: T.norm_ratio_max(2)
    10/3
    sage: T.norm_ratio_max(3)
    25/7
    sage: n(_)
    3.57142857142857
    sage: T.norm_ratio_max(4)
    62/17
    sage: n(_)
    3.64705882352941

::

    sage: T = TransformationARP_multi(3)
    sage: T.norm_ratio_max(1)
    3
    sage: T.norm_ratio_max(2)
    7/2
    sage: T.norm_ratio_max(3)
    48/13
    sage: n(_)
    3.69230769230769
    sage: T.norm_ratio_max(4)
    161/43
    sage: n(_)
    3.74418604651163

The limit 4 is achievable (when m->oo, not much when the level increases) ::

    sage: sup = (vector((0,0,1)) * A2m^10 * P1 * vector((1,1,1))).expand()
    sage: inf = (vector((0,0,1)) * A2m^10 * P1 * vector((0,0,1))).expand()
    sage: sup
    4*m^10 + 3*m^9 + 36*m^8 + 24*m^7 + 112*m^6 + 63*m^5 + 140*m^4 + 60*m^3
    + 60*m^2 + 15*m + 3
    sage: inf
    m^10 + m^9 + 9*m^8 + 8*m^7 + 28*m^6 + 21*m^5 + 35*m^4 + 20*m^3 + 15*m^2
    + 5*m + 1
    sage: lim(sup/inf, m=oo)
    4

::

    sage: sup = (vector((0,0,1)) * A2m^3 * P1 * vector((1,1,1))).expand()
    sage: inf = (vector((0,0,1)) * A2m^3 * P1 * vector((0,0,1))).expand()
    10:47:45 ###
    sage: sup
    4*m^3 + 3*m^2 + 8*m + 2
    sage: inf
    m^3 + m^2 + 2*m + 1
    sage: lim(sup/inf, m=oo)
    4



"""

import itertools
from collections import Counter

sqrt3 = sqrt(3)
M2to3 = matrix(3,[-sqrt3,-1,sqrt3,-1,0,2],ring=RR)/3
M3to2 = matrix(2,[-sqrt3,sqrt3,0,-1,-1,2],ring=RR)/2

##################
# matrices
##################
A1 = matrix(3, [1,0,0,0,1,0,-1,-1,1]).inverse()
A2 = matrix(3, [1,0,0,-1,-1,1,0,1,0]).inverse()
A3 = matrix(3, [-1,-1,1,1,0,0,0,1,0]).inverse()
P1 = matrix(3, [0,-1,1,1,0,0,-1,1,0]).inverse()
P2 = matrix(3, [-1,1,0,0,-1,1,1,0,0]).inverse()
P3 = matrix(3, [0,-1,1,-1,1,0,1,0,0]).inverse()

m = var('m')

A1m = matrix(3, [1,0,0,0,1,0,m,m,1])
A2m = matrix(3, [1,0,0,0,0,1,m,1,m])
A3m = matrix(3, [0,1,0,0,0,1,1,m,m])

Delta = matrix(3, [1,0,1,0,1,1])

B1 = matrix(3, [1,0,0, 0,1,0, 0,-1,1]).inverse()
B2 = matrix(3, [1,0,0, 0,-1,1, 0,1,0]).inverse()
B3 = matrix(3, [0,-1,1, 1,0,0, 0,1,0]).inverse()

def H(j,k, normalize=False):
    r"""
    EXAMPLES:

    These are the six inside triangles::

        sage: possible = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
        sage: for (j,k) in possible: print (j,k);print P_mat(j,k) * H(j,k)
        (3, 1)
        [  1   1   1]
        [1/2   1   0]
        [1/2   1   1]
        (2, 1)
        [  1   1   1]
        [1/2   1   1]
        [1/2   0   1]
        (3, 2)
        [  1 1/2   0]
        [  1   1   1]
        [  1 1/2   1]
        (1, 2)
        [  1 1/2   1]
        [  1   1   1]
        [  0 1/2   1]
        (2, 3)
        [  1   0 1/2]
        [  1   1 1/2]
        [  1   1   1]
        (1, 3)
        [  1   1 1/2]
        [  0   1 1/2]
        [  1   1   1]

    These are the six outside triangles::

        sage: possible = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
        sage: for (j,k) in possible: print (j,k);print AR_mat(j) * H(j,k)
        (3, 1)
        [1/2   0   0]
        [1/2   1   0]
        [  1   1   1]
        (2, 1)
        [1/2   0   0]
        [  1   1   1]
        [1/2   0   1]
        (3, 2)
        [  1 1/2   0]
        [  0 1/2   0]
        [  1   1   1]
        (1, 2)
        [  1   1   1]
        [  0 1/2   0]
        [  0 1/2   1]
        (2, 3)
        [  1   0 1/2]
        [  1   1   1]
        [  0   0 1/2]
        (1, 3)
        [  1   1   1]
        [  0   1 1/2]
        [  0   0 1/2]
    """
    if normalize:
        a = 1/2
    else:
        a = 1
    if j==3 and k==1:
        return matrix([(a,a,0), (0,1,0), (0,0,1)]).transpose()
    elif j==2 and k==1:
        return matrix([(a,0,a), (0,1,0), (0,0,1)]).transpose()
    elif j==3 and k==2:
        return matrix([(1,0,0), (a,a,0), (0,0,1)]).transpose()
    elif j==1 and k==2:
        return matrix([(1,0,0), (0,a,a), (0,0,1)]).transpose()
    elif j==2 and k==3:
        return matrix([(1,0,0), (0,1,0), (a,0,a)]).transpose()
    elif j==1 and k==3:
        return matrix([(1,0,0), (0,1,0), (0,a,a)]).transpose()
    else:
        raise ValueError("invalid j(=%s) and k(=%s)" % (j,k))

##################
# helper functions
##################
def norm_ratio(m):
    r"""
    1 Avril 2014. L'ancien ratio n'était pas le bon. Je n'utilisais pas les
    bonnes normes.

    EXAMPLES::

        sage: M = matrix(3, (1,2,3,4,5,6,7,8,9))
        sage: M
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: norm_ratio(M)
        2.66666666667
        sage: (7+8+9) / 9.
        2.66666666666667
    """
    last_row = vector((0,0,1)) * m
    return max(last_row) / min(last_row)

def is_pisot(m):
    r"""
    """
    S = sorted((abs(e) for e in m.eigenvalues()), reverse=True)
    return S[0] > 1 and S[1] < 1

def perron_eigenvector(M):
    r"""
    EXAMPLES::

        sage: L = [-13, 14, -14, -26, 12, -12, -18, 13, -21]
        sage: MI = matrix(3, L)
        sage: perron_eigenvector(MI)
        (-19.55058680967682?, (1, 1.617076203480530?, 2.084975261314588?))

    ::

        sage: perron_eigenvector(stochastisize(MI))
        Using numpy instead of Sage for finding eigenvectors...
        ((1.0000000000000004+0j), (0.57735026919, 0.57735026919, 0.57735026919))

    """
    try:
        #raise NotImplementedError
        vec = M.eigenvectors_right()
    except NotImplementedError:
        #print "Using numpy instead of Sage for finding eigenvectors..."
        eig, vec = numpy.linalg.eig(M)
        index = abs(eig).argmax()
        rightv = vec.transpose()[index]
        if eig[index].imag == 0:
            eig_sage = RR(eig[index].real)
            vec_sage = vector(a.real for a in rightv)
        else:
            eig_sage = CC(eig[index])
            vec_sage = vector(CC, rightv)
        return eig_sage, vec_sage
    else:
        rep = max(vec, key=lambda a:abs(a[0].n()))
        a,[perronvec],mul = rep
        return a, perronvec

def semi_norm_v(M, v,  p=2, verbose=False):
    r"""
    Return the semi norm on the hyperplane orthogonal to v.

    EXAMPLES::

        sage: A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
        sage: semi_norm_v(A1, vector( (1,1,1)))
        0.9999999999890247
        sage: semi_norm_v(A1, vector( (1,1,1)), p=1)
        0.9999394820959548
        sage: semi_norm_v(A1, vector( (1,1,1)), p=oo)
        1.0

    """
    def func(z):
        vz = vector(z)
        return - (M*vz).norm(p) / vz.norm(p)
    cons = [lambda z: v * vector(z),
            lambda z: - v * vector(z)]
    x0 = range(len(v))
    x0[0] = v[1]
    x0[1] = -v[0]
    rep = minimize_constrained(func, cons, x0)
    if verbose:
        print rep, rep.norm(), rep*v
    return -func(rep)

def semi_norm_cone(M, cone,  p=2, verbose=False):
    r"""
    Return the semi norm on the hyperplane orthogonal to v where v lives in
    the cone.

    EXAMPLES:

    For Arnoux-Rauzy, only the 1-norm works::

        sage: A1 = matrix(3, [1,1,1, 0,1,0, 0,0,1])
        sage: cone = A1
        sage: semi_norm_cone(A1.transpose(), cone, p=1)
        0.9999999999999998
        sage: semi_norm_cone(A1.transpose(), cone, p=oo)
        1.9999757223144654
        sage: semi_norm_cone(A1.transpose(), cone, p=2)
        1.3065629648763757

    For Poincaré, all norms work::

        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: cone = P21 * H(2,1)
        sage: semi_norm_cone(P21.transpose(), cone, p=1)
        0.9999957276014074
        sage: semi_norm_cone(P21.transpose(), cone, p=oo)
        1.0
        sage: semi_norm_cone(P21.transpose(), cone, p=2)
        0.9999999999670175

    For Poincaré on the whole domain, it works for some norms::

        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: cone = P21
        sage: semi_norm_cone(P21.transpose(), cone, p=1)
        0.9999944176794844
        sage: semi_norm_cone(P21.transpose(), cone, p=2)
        1.6180339887021953
        sage: semi_norm_cone(P21.transpose(), cone, p=oo)
        1.0

    For a product, all norms work::

        sage: A1 = matrix(3, [1,1,1, 0,1,0, 0,0,1])
        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: M = A1 * P21
        sage: cone = A1 * P21 * H(2,1)
        sage: semi_norm_cone(M.transpose(), cone, p=1)
        0.999993244882415
        sage: semi_norm_cone(M.transpose(), cone, p=oo)
        0.9999935206958908
        sage: semi_norm_cone(M.transpose(), cone, p=2)
        0.7529377601317161
    """
    a,b,c = cone.columns()
    ab = vector(matrix((a,b)).right_kernel_matrix().row(0))
    ac = vector(matrix((a,c)).right_kernel_matrix().row(0))
    bc = vector(matrix((b,c)).right_kernel_matrix().row(0))
    middle = a+b+c
    cons = []
    if ab * middle < 0:
        cons.append(lambda z: -ab * vector(z))
    else:
        cons.append(lambda z: ab * vector(z))
    if ac * middle < 0:
        cons.append(lambda z: -ac * vector(z))
    else:
        cons.append(lambda z: ac * vector(z))
    if bc * middle < 0:
        cons.append(lambda z: -bc * vector(z))
    else:
        cons.append(lambda z: bc * vector(z))
    if not all(con(middle) > 0 for con in cons):
        raise ValueError("the middle should be in the cone")
    func = lambda v : - semi_norm_v(M,vector(v),p)
    x0 = middle
    rep = minimize_constrained(func, cons, x0)
    if not all((con(rep) >= 0 or abs(con(rep)) < 1e-7) for con in cons):
        raise ValueError("the answer (={}) should be in the cone".format(rep))
    if not all(r >= 0 or abs(r) < 1e-7 for r in rep):
        raise ValueError("the answer (={}) should be positive".format(rep))
    if verbose:
        print "optimal found at ", rep / rep.norm(p)
    return -func(rep)

######################
# class Transformation
######################
class Transformation:
    r"""
    (What is the good name for this class?)

    INPUT:

    - ``M`` -- list, tuple or dict; the matrices. Keys 0,...,n-1 are used
      for list and tuple.
    - ``domain`` -- dict or matrix; the domain
    - ``language`` -- regular language or None; if None, the language is
      the full shift.

    EXAMPLES::

        sage: B1 = matrix(3, [1,0,0, 0,1,0, 0,1,1])
        sage: B2 = matrix(3, [1,0,0, 0,0,1, 0,1,1])
        sage: B3 = matrix(3, [0,1,0, 0,0,1, 1,0,1])
        sage: M = {'B1':B1, 'B2':B2, 'B3':B3}
        sage: domain = matrix(3, [1,1,1,0,1,1,0,0,1])
        sage: T = Transformation(M, domain)
    """
    def __init__(self, M, domain, language=None):
        if isinstance(M, dict):
            self._M = M
        elif isinstance(M, (list, tuple)):
            self._M = dict(enumerate(M))
        else:
            raise ValueError("M must be a list, tuple or a dict")
        if isinstance(domain, dict):
            self._domain = domain
        else:
            self._domain = {letter:domain for letter in self._M.keys()}
        if language is None:
            self._language = Language(self._M.keys())
        else:
            self._language = language

    @cached_method
    def identity_matrix(self):
        return self._M.values()[0].parent().one()

    def word_to_matrix(self, w):
        r"""
        EXAMPLES::

            sage: ARPunsorted.word_to_matrix(Word())
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return prod((self._M[a] for a in w), z=self.identity_matrix())
    def word_last_letter_to_domain(self, w):
        r"""
        EXAMPLES::

            sage: ARPunsorted.word_last_letter_to_domain(Word())
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if w:
            return self._domain[w[-1]]
        else:
            return self.identity_matrix()

    def n_words_iterator(self, n):
        r"""
        EXAMPLES::
            
            sage: list(ARP.n_words_iterator(1))
            [word: 0, word: 1, word: 2, word: 3, word: 4, word: 5]
        """
        return self._language.iterate_by_length(n)

    def n_matrices_iterator(self, n):
        r"""
        EXAMPLES::

            sage: ARP = TransformationARP()
            sage: A,B = zip(*list(ARP.n_matrices_iterator(1)))
            sage: A
            (word: P2, word: P3, word: P1, word: A1, word: A3, word: A2)
            sage: B
            (
            [0 0 1]  [0 0 1]  [0 1 0]  [1 0 0]  [0 1 0]  [1 0 0]
            [1 0 1]  [0 1 1]  [0 1 1]  [0 1 0]  [0 0 1]  [0 0 1]
            [1 1 1], [1 1 1], [1 1 1], [1 1 1], [1 1 1], [1 1 1]
            )
        """
        for w in self.n_words_iterator(n):
            yield w, self.word_to_matrix(w)
    def n_cylinders_iterator(self, n):
        r"""
        EXAMPLES::

            sage: it = ARPunsorted.n_cylinders_iterator(1)
            sage: for w,cyl in it: print "{}\n{}".format(w,cyl)
            A1
            [1 1 1]
            [0 1 0]
            [0 0 1]
            A2
            [1 0 0]
            [1 1 1]
            [0 0 1]
            A3
            [1 0 0]
            [0 1 0]
            [1 1 1]
            P12
            [1 1 1]
            [1 2 1]
            [0 1 1]
            P13
            [1 1 1]
            [0 1 1]
            [1 1 2]
            P23
            [1 0 1]
            [1 1 1]
            [1 1 2]
            P21
            [2 1 1]
            [1 1 1]
            [1 0 1]
            P31
            [2 1 1]
            [1 1 0]
            [1 1 1]
            P32
            [1 1 0]
            [1 2 1]
            [1 1 1]
            
        """
        for w in self.n_words_iterator(n):
            yield w, self.word_to_matrix(w)*self.word_last_letter_to_domain(w)

    def plot_partition(self, n, label=True):
        r"""
        EXAMPLES::

            sage: ARPunsorted.plot_partition(3)

        """
        G = Graphics()
        for w,cyl in self.n_cylinders_iterator(n):
            columns = cyl.columns()
            G += polygon((M3to2*col/col.norm(1) for col in columns), fill=False) 
            if label:
                sum_cols = sum(columns)
                G += text("{}".format(w), M3to2*sum_cols/sum_cols.norm(1))
        return G

    def is_pisot(self, w, verbose=False):
        r"""
        """
        m = self.word_to_matrix(w)
        S = sorted((abs(e) for e in m.eigenvalues()), reverse=True)
        if S[0] > 1 and S[1] < 1:
            return True
        else:
            if verbose:
                print "not pisot:", S[0], S[1]
                print "indices of matrices:", w
                print m
            return False

    @cached_method
    def non_pisot_matrices(self,n, verbose=False):
        r"""
        Return the list of non pisot matrices (as list of indices of base
        matrices).

        EXAMPLES::

            sage: ARP.non_pisot_matrices(1)
            [word: A1, word: A2]
            sage: ARP.non_pisot_matrices(2)
            [word: A1,A1, word: A1,A2, word: A2,A1, word: A2,A2]
            sage: ARP.non_pisot_matrices(3)
            [word: A1,A1,A1,
             word: A1,A1,A2,
             word: A1,A2,A1,
             word: A1,A2,A2,
             word: A2,A1,A1,
             word: A2,A1,A2,
             word: A2,A2,A1,
             word: A2,A2,A2]
            sage: len(ARP.non_pisot_matrices(4))  # long time
            16

        ::

            sage: Brun.non_pisot_matrices(2)
            [word: B1,B1, word: B1,B2, word: B2,B1, word: B2,B2]
            sage: Brun.non_pisot_matrices(3)
            [word: B1,B1,B1,
             word: B1,B1,B2,
             word: B1,B2,B1,
             word: B1,B2,B2,
             word: B2,B1,B1,
             word: B2,B1,B2,
             word: B2,B2,B1,
             word: B2,B2,B2]
        """
        return [w for w in self.n_words_iterator(n) if not self.is_pisot(w, verbose=verbose)]

    def non_pisot_automaton(self, n):
        r"""
        EXAMPLES::

            sage: A = ARPunsorted.non_pisot_automaton(2)
            sage: A
            Automaton with 9 states
            sage: A.graph().plot(edge_labels=True)   # not tested
        """
        L = []
        for i in range(n):
            L.extend(self.non_pisot_matrices(i))
        alphabet = self._language._alphabet
        F = FiniteLanguage(alphabet, L)
        A = F.automaton()
        A = A.minimization().relabeled()
        #return A
        G = A.graph()
        to_remove = set(A5.states()) - set(A5.final_states())
        G.delete_vertices(to_remove)
        return G

    def eigenvalues_n_matrices(self,n, verbose=False):
        r"""
        Return the eigenvalues of the matrices of level n.

        EXAMPLES::

            sage: ARP.non_pisot_matrices(1)
            [(0,), (1,)]
            sage: ARP.non_pisot_matrices(2)
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: ARP.non_pisot_matrices(3)       # long time
            [(0, 0, 0),
             (0, 0, 1),
             (0, 1, 0),
             (0, 1, 1),
             (1, 0, 0),
             (1, 0, 1),
             (1, 1, 0),
             (1, 1, 1)]
            sage: len(ARP.non_pisot_matrices(4))  # long time
            16

        ::

            sage: Brun.non_pisot_matrices(2)
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: Brun.non_pisot_matrices(3)
            [(0, 0, 0),
             (0, 0, 1),
             (0, 1, 0),
             (0, 1, 1),
             (1, 0, 0),
             (1, 0, 1),
             (1, 1, 0),
             (1, 1, 1)]
        """
        R = []
        for w in self.n_words_iterator(n):
            m = self.word_to_matrix(w)
            S = m.eigenvalues()
            R.append((m,S))
            if verbose:
                print "indices of matrices:", w
                print m
                print "eigenvalues:", S
        return R

    def left_right_eigenvectors_n_matrices(self,n, verbose=False):
        r"""
        Return the left and right eigenvectors of the matrices of level n.

        EXAMPLES::

            sage: ARPunsorted.left_right_eigenvectors_n_matrices(1)
            [(word: P12, (0, 1, 0), (0, 0, 1)),
             (word: P13, (0, 0, 1), (0, 1, 0)),
             (word: P23, (0, 0, 1), (1, 0, 0)),
             (word: P21, (1, 0, 0), (0, 0, 1)),
             (word: P31, (1, 0, 0), (0, 1, 0)),
             (word: P32, (0, 1, 0), (1, 0, 0))]

        """
        R = []
        for w in self.n_words_iterator(n):
            m = prod(self._M[i] for i in w)
            try:
                a,v_right = perron_eigenvector(m)
                b,v_left = perron_eigenvector(m.transpose())
            except ValueError:
                print "problem with :\n",m
            else:
                R.append((w, v_right,v_left))
                if verbose:
                    print "indices of matrices:", w
                    print m
                    print "eigenvectors:", v_right, v_left
        return R
    def plot_eigenvectors_n_matrices(self, n, side='right', color_index=0, draw_line=False):
        r"""
        INPUT:

        - color_index -- 0 for first letter, -1 for last letter

        EXAMPLES::

            sage: ARP.plot_left_right_eigenvectors_n_matrices(2)
        """
        R = self.left_right_eigenvectors_n_matrices(n)
        L = [(w, M3to2*(a/sum(a)), M3to2*(b/sum(b))) for (w,a,b) in R]
        G = Graphics()
        alphabet = self._language._alphabet
        color_ = dict( (letter, hue(i/float(len(alphabet)))) for i,letter in
                enumerate(alphabet))
        for letter in alphabet:
            L_filtered = [(w,p1,p2) for (w,p1,p2) in L if w[color_index] == letter]
            words,rights,lefts = zip(*L_filtered)
            if side == 'right':
                G += point(rights, color=color_[letter], legend_label=letter)
            elif side == 'left':
                G += point(lefts,  color=color_[letter], legend_label=letter)
            else:
                raise ValueError("side(=%s) should be left or right" % side)

        if draw_line:
            for (a,b) in L:
                G += line([a,b], color='black', linestyle=":")
        G += line([M3to2*vector(a) for a in [(1,0,0), (0,1,0), (0,0,1), (1,0,0)]]) 
        title = "%s eigenvectors, colored by letter w[%s] of cylinder w" % (side, color_index)
        G += text(title, (0.5, 1.05), axis_coords=True)
        G.axes(False)
        return G

    def plot_pisot_conjugates(self, n):
        r"""
        EXAMPLES::

            sage: Brun.plot_pisot_conjugates(5)

        Image envoyee a Timo (6 mai 2014)::

            sage: sum(Brun.plot_pisot_conjugates(i) for i in [1..6])  #not tested
        """
        Lreal = []
        Limag = []
        for w,m in self.n_matrices_iterator(n):
            a,b,c = sorted(m.eigenvalues(), key=abs)
            if a.imag() == 0 and b.imag() == 0:
                Lreal.append((a,b))
            else:
                Limag.append((a.real(),a.imag()))
                Limag.append((b.real(),b.imag()))
        return points(Lreal) + points(Limag, color='red')

    @cached_method
    def plot_pisot_conjugates_brun(self, n):
        r"""
        EXAMPLES::

            sage: Brun.plot_pisot_conjugates(5)

        Image envoyee a Timo (6 mai 2014)::

            sage: sum(Brun.plot_pisot_conjugates(i) for i in [1..6])  #not tested
        """
        T = [('a','a',0), ('a','b',1), ('a','d',2)]
        T += [('b','b',0), ('b','a',1)]
        T += [('c','c',0), ('c','b',2)]
        T += [('d','c',1), ('d','c',2)]
        A = Automaton(T, initial_states=['a','b','c','d'],
                      final_states=['a','b','c','d'])
        B = A.determinisation()
        LrealIN = []
        LrealNO = []
        Limag = []
        for w in self.n_words_iterator(n):
            m = prod(self._M[i] for i in w)
            a,b,c = sorted(m.eigenvalues(), key=abs)
            if a.imag() == 0 and b.imag() == 0:
                if B(w*100)[0]:
                    LrealIN.append((a,b))
                else:
                    LrealNO.append((a,b))
            else:
                Limag.append((a.real(),a.imag()))
                Limag.append((b.real(),b.imag()))
        return (points(LrealNO, color='green') 
            + points(LrealIN, color='blue') 
            + points(Limag, color='red'))



    def position_of_last_column(self, n):
        r"""
        EXAMPLES:

        On dirait que la troisieme colonne n'est jamais la plus petite.
        Soit la plus grande (pos=2) soit la milieu (pos=1).

            sage: T = TransformationARP_multi(4)
            sage: T.position_of_last_column(1)
            Counter({1: 13, 2: 10})
            sage: T.position_of_last_column(2)
            Counter({2: 342, 1: 187})
            sage: T.position_of_last_column(3)
            Counter({2: 7866, 1: 4301})

        ::

            sage: T = TransformationARP_multi(6)
            sage: T.position_of_last_column(1)
            Counter({1: 19, 2: 14})
            sage: T.position_of_last_column(2)
            Counter({2: 702, 1: 387})
            sage: T.position_of_last_column(3)
            Counter({2: 23166, 1: 12771})
        """
        c = Counter()
        for m in self.m_matrices(n):
            norms = tuple(col.norm(1) for col in m.columns())
            index = sorted(norms).index(norms[-1])
            c[index] += 1
        return c

    def norm_ratio_iterator(self,n):
        r"""
        returns an iterator of the ratio max/min of the norm of the columns
        of the n-cylinders.
        """
        for w,m in self.n_cylinders_iterator(n):
            yield norm_ratio(m)
    def semi_norm_study(self, n, p=2):
        r"""
        EXAMPLES:

        For the 1-norm, all matrices contracts the hyperplane::
            
            sage: ARPunsorted.semi_norm_study(1, p=1)
            A1 1.0 False
            A2 1.0 False
            A3 1.0 False
            P12 0.999990663839 False
            P13 0.999987801801 False
            P23 0.999983738803 False
            P21 0.999995727601 False
            P31 0.999994171082 False
            P32 0.999996491473 False

        For the 2-norm, AR matrices do not contract::

            sage: ARPunsorted.semi_norm_study(1, p=2)
            A1 1.30656296488 False
            A2 1.30656296486 False
            A3 1.30656296475 False
            P12 0.99999999996 False
            P13 0.999999999967 False
            P23 0.999999999997 False
            P21 0.999999999967 False
            P31 0.999999999769 False
            P32 0.999999999839 False

        When, the 1-norm is < 1, the product is pisot::

            sage: ARPunsorted.semi_norm_study(2, p=1)
            A1,A1 1.0 False
            A1,A2 1.0 False
            A1,A3 1.0 False
            A1,P12 0.999998922557 False
            A1,P13 0.999997464905 False
            A1,P23 0.999999150973 True
            A1,P21 0.999993244882 False
            A1,P31 0.999994030522 False
            A1,P32 0.999998046513 True
            A2,A1 1.0 False
            A2,A2 1.0 False
            A2,A3 1.0 False
            A2,P12 0.99999375291 False
            A2,P13 0.999995591588 True
            ...
            P31,A3 0.999988326888 False
            P31,P12 0.749998931902 True
            P31,P23 0.799999157344 True
            P31,P32 0.749993104833 True
            P32,A1 0.999997170005 True
            P32,A3 0.99999420509 False
            P32,P13 0.666665046248 True
            P32,P21 0.666665629351 True
            P32,P31 0.666664488371 True
        """
        for w in self.n_words_iterator(n):
            m = self.word_to_matrix(w)
            cone = m*self.word_last_letter_to_domain(w)
            print w, semi_norm_cone(m.transpose(), cone, p=p), self.is_pisot(w)

    def norm_ratio_max(self, n):
        r"""
        EXAMPLES:

        Non borné::

            sage: T = TransformationARP()
            sage: T.norm_ratio_max(1)
            6.0
            sage: T.norm_ratio_max(2)
            9.0
            sage: T.norm_ratio_max(3)
            12.0
            sage: T.norm_ratio_max(4)
            15.0
            sage: T.norm_ratio_max(5)
            18.0

        """
        return max(self.norm_ratio_iterator(n))


    def cylinder_with_max_ratio(self,n):
        r"""
        """
        it = self.n_words_iterator(n)
        key = lambda w:norm_ratio(prod(self._M[a] for a in w)*self._domain)
        return max(it, key=key)

    def triangle_partition(self,n):
        L = [Triangle(*m.columns()) for w,m in self.n_cylinders_iterator(n)]
        return EnsembleDeTriangles(L)





####################
# ARP Language
####################
def ARP_automaton():
    r"""
    EXAMPLES::

        sage: A = ARP_automaton(None)
        sage: A.process(['A1', 'P12', 'A1', 'P13'])
        (True, 'H13')
        sage: A(['A1', 'P12', 'A1', 'P13'])
        True
    """
    import sage.combinat.finite_state_machine
    sage.combinat.finite_state_machine.FSMOldProcessOutput = False
    def H(j,k):
        return 'H%s%s' % (j,k)
    def P(j,k):
        return 'P%s%s' % (j,k)
    def A(k):
        return 'A%s' % (k)
    jk = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
    D = 'Delta'
    states = [H(j,k) for (j,k) in jk] + [D]
    autom = Automaton(initial_states=[D], final_states=states)
    for (j,k) in jk:
        i = 6 - j - k
        v = H(j,k)
        autom.add_transition(v, H(i,j), P(i,j))
        autom.add_transition(v, H(k,i), P(k,i))
        autom.add_transition(v, H(j,i), P(j,i))
        autom.add_transition(v, D, A(i))
        autom.add_transition(v, v, A(j))
        autom.add_transition(D, v, P(j,k))
    for k in [1,2,3]:
        autom.add_transition(D, D, A(k))
    return autom

class Language(object):
    def __init__(self, alphabet):
        self._alphabet = alphabet

    def iterate_by_length(self, length):
        W = Words(self._alphabet)
        it = W.iterate_by_length(length)
        return it

class RegularLanguage(Language):
    r"""
    EXAMPLES::

        sage: alphabet = ['A1', 'A2', 'A3', 'P31', 'P21', 'P32', 'P12', 'P23', 'P13']
        sage: automaton = ARP_automaton()
        sage: L = RegularLanguage(alphabet, automaton)
        sage: len(list(L.iterate_by_length(0)))
        1
        sage: len(list(L.iterate_by_length(1)))
        9
        sage: len(list(L.iterate_by_length(2)))
        57
        sage: len(list(L.iterate_by_length(3)))
        345

    """
    def __init__(self, alphabet, automaton):
        Language.__init__(self, alphabet)
        self._automaton = automaton

    def iterate_by_length(self, length):
        W = Words(self._alphabet)
        it = W.iterate_by_length(length)
        return itertools.ifilter(self._automaton, it)
class FiniteLanguage(Language):
    r"""
    EXAMPLES::

        sage: L = ['a', 'aa', 'aaa']
        sage: F = FiniteLanguage(alphabet=['a'], words=L)
    """
    def __init__(self, alphabet, words):
        Language.__init__(self, alphabet)
        self._words = words
    def automaton(self):
        r"""
        EXAMPLES::

            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.automaton()
            Automaton with 7 states
        """
        transitions = []
        final_states = []
        end = 1
        for w in self._words:
            start = 0
            for a in w:
                transitions.append((start, end, a))
                start, end = end, end+1
            final_states.append(start)
        return Automaton(transitions, initial_states=[0], final_states=final_states)
    def minimal_automaton(self):
        r"""
        .. NOTE:: 
        
            One of the state is not final. You may want to remove it...

        EXAMPLES::

            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.minimal_automaton()
            Automaton with 5 states
        """
        A = self.automaton()
        A = A.minimization().relabeled()
        return A
    def number_of_states(self):
        r"""
        EXAMPLES::

            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.number_of_states()
            5
        """
        return len(self.minimal_automaton().states())

####################
# quick construction
####################
def TransformationARPunsorted():
    A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
    A2 = matrix(3, [1,0,0, -1,1,-1, 0,0,1]).inverse()
    A3 = matrix(3, [1,0,0, 0,1,0, -1,-1,1]).inverse()
    P12 = matrix(3, [1, 0, 1, 1, 1, 1, 0, 0, 1])
    P13 = matrix(3, [1, 1, 0, 0, 1, 0, 1, 1, 1])
    P23 = matrix(3, [1, 0, 0, 1, 1, 0, 1, 1, 1])
    P21 = matrix(3, [1, 1, 1, 0, 1, 1, 0, 0, 1])
    P31 = matrix(3, [1, 1, 1, 0, 1, 0, 0, 1, 1])
    P32 = matrix(3, [1, 0, 0, 1, 1, 1, 1, 0, 1])
    M = (A1, A2, A3, P12, P13, P23, P21, P31, P32)
    alphabet = ['A1', 'A2', 'A3', 'P12', 'P13', 'P23', 'P21', 'P31', 'P32']
    M = dict(zip(alphabet, M))
    automaton = ARP_automaton()
    L = RegularLanguage(alphabet, automaton)
    pairs = [(1,2), (1,3), (2,3), (2,1), (3,2), (3,1)]
    domain = {"P{}{}".format(j,k):H(j,k) for (j,k) in pairs}
    domain.update({"A{}".format(k):identity_matrix(3) for k in [1,2,3]})
    return Transformation(M, domain, language=L)

def TransformationArnouxRauzy():
    A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
    A2 = matrix(3, [1,0,0, -1,1,-1, 0,0,1]).inverse()
    A3 = matrix(3, [1,0,0, 0,1,0, -1,-1,1]).inverse()
    M = (A1, A2, A3)
    alphabet = ['A1', 'A2', 'A3']
    M = dict(zip(alphabet, M))
    domain = matrix(3, [1,0,0,0,1,0,0,0,1])
    return Transformation(M, domain)

def TransformationBrun():
    B1 = matrix(3, [1,0,0, 0,1,0, 0,-1,1]).inverse()
    B2 = matrix(3, [1,0,0, 0,-1,1, 0,1,0]).inverse()
    B3 = matrix(3, [0,-1,1, 1,0,0, 0,1,0]).inverse()
    M = (B1, B2, B3)
    alphabet = ['B1', 'B2', 'B3']
    M = dict(zip(alphabet, M))
    domain = matrix(3, [1,1,1,0,1,1,0,0,1])
    return Transformation(M, domain)

def TransformationARP():
    A1 = matrix(3, [1,0,0,0,1,0,-1,-1,1]).inverse()
    A2 = matrix(3, [1,0,0,-1,-1,1,0,1,0]).inverse()
    A3 = matrix(3, [-1,-1,1,1,0,0,0,1,0]).inverse()
    P1 = matrix(3, [0,-1,1,1,0,0,-1,1,0]).inverse()
    P2 = matrix(3, [-1,1,0,0,-1,1,1,0,0]).inverse()
    P3 = matrix(3, [0,-1,1,-1,1,0,1,0,0]).inverse()
    M = (A1, A2, A3, P1, P2, P3)
    alphabet = ['A1', 'A2', 'A3', 'P1', 'P2', 'P3']
    M = dict(zip(alphabet, M))
    domain = matrix(3, [1,0,0,1,1,0,1,1,1])
    return Transformation(M, domain)

def TransformationARP_multi(order):
    domain = matrix(3, [1,0,0,1,1,0,1,1,1])
    M = [P1,P2,P3]
    t12 = matrix(3, [1,0,0,0,0,1,0,1,0])
    t132 = matrix(3, [0,1,0,0,0,1,1,0,0])
    for i in range(1, order+1):
        A1_i = A1**i
        M.append(A1_i * P1)
        M.append(A1_i * P2)
        M.append(A1_i * P3)
        M.append(A1_i * t12)
        M.append(A1_i * t132)
    return Transformation(M, domain)

AR = TransformationArnouxRauzy()
Brun = TransformationBrun()
ARP = TransformationARP()
ARPunsorted = TransformationARPunsorted()
ARPm = TransformationARP_multi(4)

####################
# Polyhedron Partition
####################
def arp_polyhedron(d=3):
    r"""
    Return the d-dimensional 1-cylinders of the ARP algorithm.

    EXAMPLES::

        sage: A,P = arp_polyhedron(3)
        sage: A.vertices_list()
        [[0, 0, 0], [1/2, 1/2, 0], [1/2, 1/4, 1/4], [1, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0], [1/2, 1/2, 0], [1/2, 1/4, 1/4], [1/3, 1/3, 1/3]]

    ::

        sage: A,P = arp_polyhedron(4)
        sage: A.vertices_list()
        [[0, 0, 0, 0],
         [1/2, 1/2, 0, 0],
         [1/2, 1/6, 1/6, 1/6],
         [1/2, 1/4, 1/4, 0],
         [1, 0, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0, 0],
         [1/2, 1/2, 0, 0],
         [1/2, 1/4, 1/4, 0],
         [1/2, 1/6, 1/6, 1/6],
         [1/4, 1/4, 1/4, 1/4],
         [1/3, 1/3, 1/3, 0]]

    ::

        sage: A,P = arp_polyhedron(5)
        sage: A.vertices_list()
        [[0, 0, 0, 0, 0],
         [1/2, 1/2, 0, 0, 0],
         [1/2, 1/8, 1/8, 1/8, 1/8],
         [1/2, 1/6, 1/6, 1/6, 0],
         [1/2, 1/4, 1/4, 0, 0],
         [1, 0, 0, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0, 0, 0],
         [1/2, 1/2, 0, 0, 0],
         [1/2, 1/6, 1/6, 1/6, 0],
         [1/2, 1/8, 1/8, 1/8, 1/8],
         [1/2, 1/4, 1/4, 0, 0],
         [1/3, 1/3, 1/3, 0, 0],
         [1/5, 1/5, 1/5, 1/5, 1/5],
         [1/4, 1/4, 1/4, 1/4, 0]]
    """
    positive = [ [0]*i + [1] + [0]*(d-i) for i in range(1, d+1)]
    atmostone = [[1] + [-1]*d]
    ieq_arnoux = [[0]+[1]+[-1]*(d-1)]
    ieq_arnoux_not = [[0]+[-1]+[1]*(d-1)]
    ieq_sorted = [ [0]*i + [1,-1] + [0]*(d-i-1) for i in range(1,d)]
    A = Polyhedron(ieqs=positive + atmostone + ieq_sorted + ieq_arnoux)
    P = Polyhedron(ieqs=positive + atmostone + ieq_sorted + ieq_arnoux_not)
    L = Polyhedron(ieqs=positive + atmostone + ieq_sorted)
    return A,P,L

def cassaigne_polyhedron(d=3):
    r"""
    Return the d-dimensional 1-cylinders of the Cassaigne algorithm.

    (of the dual!)

    EXAMPLES::

        sage: L,La,Lb = cassaigne_polyhedron(3)
        sage: L.vertices_list()
        [[0, 0, 0], [0, 1/2, 1/2], [1/3, 1/3, 1/3], [1/2, 1/2, 0]]
        sage: La.vertices_list()
        [[0, 0, 0], [0, 1/2, 1/2], [1/3, 1/3, 1/3], [1/4, 1/2, 1/4]]
        sage: Lb.vertices_list()
        [[0, 0, 0], [1/3, 1/3, 1/3], [1/2, 1/2, 0], [1/4, 1/2, 1/4]]

    ::

        sage: L,La,Lb = cassaigne_polyhedron(4)
        sage: L.vertices_list()
        [[0, 0, 0, 0],
         [0, 1/3, 1/3, 1/3],
         [1/3, 1/3, 1/3, 0],
         [1/4, 1/4, 1/4, 1/4],
         [1/5, 2/5, 1/5, 1/5],
         [1/5, 1/5, 2/5, 1/5]]

    ::

        sage: L,La,Lb = cassaigne_polyhedron(5)
        sage: L.vertices_list()
        [[0, 0, 0, 0, 0],
         [0, 1/4, 1/4, 1/4, 1/4],
         [1/4, 1/4, 1/4, 1/4, 0],
         [1/6, 1/6, 1/3, 1/6, 1/6],
         [1/5, 1/5, 1/5, 1/5, 1/5],
         [1/6, 1/3, 1/6, 1/6, 1/6],
         [1/7, 2/7, 2/7, 1/7, 1/7],
         [1/7, 2/7, 1/7, 2/7, 1/7],
         [1/7, 1/7, 2/7, 2/7, 1/7],
         [1/6, 1/6, 1/6, 1/3, 1/6]]
    """
    # [-1,7,3,4] represents the inequality 7x_1+3x_2+4x_3>= 1.
    positive = [ [0]*i + [1] + [0]*(d-i) for i in range(1, d+1)]
    atmostone = [[1] + [-1]*d]
    ai_lt_a1d = [[0]+[1]+[0]*(i-2)+[-1]+[0]*(d-i-1)+[1] for i in range(2,d)]
    ai_gt_a1 = [[0]+[-1]+[0]*(i-2)+[1]+[0]*(d-i-1)+[0] for i in range(2,d)]
    ai_gt_ad = [[0]+[0]+[0]*(i-2)+[1]+[0]*(d-i-1)+[-1] for i in range(2,d)]
    a1_gt_ad = [[0]+[1]+[0]*(d-2)+[-1]]
    a1_lt_ad = [[0]+[-1]+[0]*(d-2)+[1]]
    L = Polyhedron(ieqs=positive + atmostone + ai_lt_a1d + ai_gt_a1 + ai_gt_ad)
    La = Polyhedron(ieqs=positive+atmostone+ai_lt_a1d+ai_gt_a1+ai_gt_ad+a1_lt_ad)
    Lb = Polyhedron(ieqs=positive+atmostone+ai_lt_a1d+ai_gt_a1+ai_gt_ad+a1_gt_ad)
    return L, La, Lb

