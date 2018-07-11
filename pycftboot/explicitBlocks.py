import scipy.misc as scp

def pochhammer(x, m):
	if m==0:
		return 1
	else:
		return reduce(lambda y, z: y*z, [x - i for i in range(m)])


class ExplicitConvolvedBlockTable:
    """
    Explicitly gives the convolved Conformal blocks

    Parameters
    ----------
    block_table: A `ConformalBlockTable` from which to produce the convolved blocks.

    Attributes
    ----------
    dim:         The spatial dimension, inherited from `block_table`.
    k_max:       Numer controlling the accuracy of the rational approximation,
                 inherited from `block_table`.
    l_max:       The highest spin kept in the convolved block table. This is at most
                 the `l_max` of `block_table`.
    m_max:       Number controlling how many `a` derivatives there are where the
                 standard co-ordinates are expressed as `(a + sqrt(b)) / 2` and
                 `(a - sqrt(b)) / 2`. This is at most the `m_max` of `block_table`.
    n_max:       The number of `b` derivatives there are where the standard
                 co-ordinates are expressed as `(a + sqrt(b)) / 2` and
                 `(a - sqrt(b)) / 2`. This is at most the `n_max` of `block_table`.
    delta_12:    The difference between the external scaling dimensions of operator
                 1 and operator 2, inherited from `block_table`.
    delta_32:    The difference b .5etween the external scaling dimensions of operator
                 3 and operator 4, inherited from `block_table`.
    table:       A list of `PolynomialVector`s. A block's position in the table is
                 equal to its spin if `odd_spins` is `True`. Otherwise it is equal
                 to half of the spin.
    m_order:     A list stating how many `a` derivatives are being described by the
                 corresponding entry in a `PolynomialVector` in `table`. Different
                 from the `m_order` of `block_table` because some derivatives vanish
                 by symmetry.
    n_order:     A list stating how many `b` derivatives are being described by the
                 corresponding entry in a `PolynomialVector` in `table`.
    v0:		 A list of exact values at delta = 0
    """
	

    def __init__(self, block_table):
        # Copying everything but the unconvolved table is fine from a memory standpoint
        self.dim = block_table.dim
        self.k_max = block_table.k_max
        self.l_max = block_table.l_max
        self.m_max = block_table.m_max
        self.n_max = block_table.n_max
        self.delta_12 = block_table.delta_12
        self.delta_34 = block_table.delta_34
        self.v0 = []
        self.m_order = []
        self.n_order = []
        self.table = []
	

        for i in range(len(block_table.m_order)):
                if block_table.m_order[i] % 2 == 1:
                        self.m_order.append(block_table.m_order[i])
                        self.n_order.append(block_table.n_order[i])

	
        indexList = zip(block_table.m_order, block_table.n_order)
        indexListConv = zip(self.m_order, self.n_order)
        poleProd = [reduce(lambda x, y: x * y, [(delta - x) for x in tab.poles]) for tab in block_table.table]
        #poleProd = [1 for tab in block_table.table]
        gamma = pow(3 - 2 * sp.sqrt(2), delta)
        #gamma = 1

        self.v0 = [-2 * pow(-1, n) * pow(4, -delta_ext) * pochhammer(delta_ext, n) * \
                pochhammer(2 * (delta_ext - n), m) for (m, n) in zip(self.m_order, self.n_order)]



        for ii in range(len(block_table.table)):

                vector = [sum([scp.comb(m, i) * scp.comb(n, j) * pow(4, -delta_ext) *  \
                        pochhammer(delta_ext, j) * pochhammer(2 * delta_ext - 2 * j, i) * pow(-1, j + i) * \
                        block_table.table[ii].vector[indexList.index((m-i,n-j))] / \
                        poleProd[ii] for i in range(m+1) for j in range(n+1)]) for (m,n) in indexListConv]

                vector2 = [sum([scp.comb(m, i) * scp.comb(n, j) * pow(4, -delta_ext) *  \
                        pochhammer(delta_ext, j) * pochhammer(2 * delta_ext - 2 * j, i) * pow(-1, m - i + j) * \
                        block_table.table[ii].vector[indexList.index((m-i,n-j))] / \
                        poleProd[ii] for i in range(m+1) for j in range(n+1)]) for (m,n) in indexListConv]

                self.table.append([gamma * (x - y) for (x, y) in zip(vector, vector2)])




class ExplicitBlockTable:
    """
    Explicitly gives the Conformal blocks

    Parameters
    ----------
    block_table: A `ConformalBlockTable` from which to produce the explicit blocks.

    Attributes
    ----------
    dim:         The spatial dimension, inherited from `block_table`.
    k_max:       Numer controlling the accuracy of the rational approximation,
                 inherited from `block_table`.
    l_max:       The highest spin kept in the convolved block table. This is at most
                 the `l_max` of `block_table`.
    m_max:       Number controlling how many `a` derivatives there are where the
                 standard co-ordinates are expressed as `(a + sqrt(b)) / 2` and
                 `(a - sqrt(b)) / 2`. This is at most the `m_max` of `block_table`.
    n_max:       The number of `b` derivatives there are where the standard
                 co-ordinates are expressed as `(a + sqrt(b)) / 2` and
                 `(a - sqrt(b)) / 2`. This is at most the `n_max` of `block_table`.
    delta_12:    The difference between the external scaling dimensions of operator
                 1 and operator 2, inherited from `block_table`.
    delta_32:    The difference b .5etween the external scaling dimensions of operator
                 3 and operator 4, inherited from `block_table`.
    table:       A list of `PolynomialVector`s. A block's position in the table is
                 equal to its spin if `odd_spins` is `True`. Otherwise it is equal
                 to half of the spin.
    m_order:     A list stating how many `a` derivatives are being described by the
                 corresponding entry in a `PolynomialVector` in `table`. Different
                 from the `m_order` of `block_table` because some derivatives vanish
                 by symmetry.
    n_order:     A list stating how many `b` derivatives are being described by the
                 corresponding entry in a `PolynomialVector` in `table`.
    """
	

    def __init__(self, block_table):
        # Copying everything but the unconvolved table is fine from a memory standpoint
        self.dim = block_table.dim
        self.k_max = block_table.k_max
        self.l_max = block_table.l_max
        self.m_max = block_table.m_max
        self.n_max = block_table.n_max
        self.delta_12 = block_table.delta_12
        self.delta_34 = block_table.delta_34

        self.m_order = block_table.m_order
        self.n_order = block_table.n_order
        self.table = []
	

        for block_table_instance in block_table.table:

                poleProd = reduce(lambda x, y: x * y, [(delta - x) for x in block_table_instance.poles])

                self.table.append([block / poleProd for block in block_table_instance.vector])




def TwoDBlock(z, Del, l, delExt):
        b1 = (Del + l) / 2
        b2 = (Del - l) / 2
        return 2 * (pow(z, b1 + 2 * delExt) * mpmath.hyp2f1(b1, b1, 2 * b1, z) * pow(z, b2) * mpmath.hyp2f1(b2, b2, 2 * b2, z) - \
                pow(1 - z, b1 + 2 * delExt) * mpmath.hyp2f1(b1, b1, 2 * b1, 1 - z) * pow(1 - z, b2) * mpmath.hyp2f1(b2, b2, 2 * b2, 1 - z))







