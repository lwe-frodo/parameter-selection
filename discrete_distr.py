"""Support functions for manipulating discrete probability distributions.

A discrete probability distribution is represented as a dictionary mapping
values in its support to probability masses. Ex: {-1: .5, +1: .5}

A distribution may be defined over Z_q, which affects how the sums and products
are evaluated.

Based on the paper:
    Joppe Bos, Craig Costello, Leo Ducas, Ilya Mironov, Michael Naehrig, Valeria
    Nikolaenko, Ananth Raghunathan, Douglas Stebila.  Frodo: Take off the ring!
    Practical, quantum-secure key exchange from LWE.  In ACM Conference on Computer
    and Communications Security (CCS) 2016, ACM, October, 2016.
    DOI: http://dx.doi.org/10.1145/2976749.2978425
    Eprint http://eprint.iacr.org/2016/659

Copyright (c) 2016 Joppe Bos, Leo Ducas, Ilya Mironov, Valeria Nikolaenko,
                   Ananth Raghunathan, Douglas Stebila

Released under the MIT License; see LICENSE.txt for details.
"""

from scipy.stats import norm, binom
from math import fmod

NEGLIGIBLE = 1E-78  # < 2^{-256}

### Utility functions ###


def valid_distr(d):
    """Sanity check on a distribution.

    Args:
       d: A distribution represented as a dictionary.

    Returns:
      True iff d is a valid distribution.
    """
    return abs(sum(d.itervalues()) - 1.) < 1E-9 and all(x >= 0
                                                        for x in d.itervalues())


def valid_symmetric_distr(d):
    """Sanity check on a distribution, and check it is symmetric

    Args:
       d: A distribution represented as a dictionary.

    Returns:
      True iff d is a valid symmetric distribution
    """

    return valid_distr(d) and all(
        abs(d[x] - d[-x]) < 1E-9 for x in d.iterkeys())


def filter_negl(d):
    """Filters out negligible values from a distribution.

    Args:
       d: A distribution represented as a dictionary.

    Returns:
      Distribution represented as a dictionary.
    """
    return {v: p for (v, p) in d.iteritems() if p >= NEGLIGIBLE}


def std_modulo(d, q):
    """Computes standard deviation of a distribution defined over Z_q.

    Args:
      d: A distribution represented as a dictionary.
      q: A modulus.

    Returns:
      The distribution's standard deviation.

    """
    return sum([p * min(x, q - x) ** 2 for (x, p) in d.iteritems()]) ** .5


### Operations on distributions ###

def convolution(a, b, q=None):
    """Computes a convolution of two probability distributions.

    Args:
      a, b: Two distributions represented as dictionaries.
      q: If q is defined, computation is done modulo q. If q is None,
        distributions are taken over integers.

    Returns:
      A distribution represented as a dictionary.
    """
    if not (valid_distr(a) and valid_distr(b)):
        raise ValueError("One of the input distributions is not valid")

    c = {}
    for v, p_v in a.iteritems():
        for u, p_u in b.iteritems():
            w = v + u
            if q is not None:
                w %= q
            if w in c:
                c[w] += p_v * p_u
            else:
                c[w] = p_v * p_u
    return filter_negl(c)


def pdf_product(a, b, q=None):
    """Computes the product of two probability distributions.

    Args:
      a, b: Two distributions represented as dictionaries.
      q: If q is defined, computation is done modulo q. If q is None,
        distributions are taken over integers.

    Returns:
      A distribution represented as a dictionary.
    """
    if not (valid_distr(a) and valid_distr(b)):
        raise ValueError("One of the input distribution is not valid")

    c = {}
    for v, p_v in a.iteritems():
        for u, p_u in b.iteritems():
            w = v * u
            if q is not None:
                w %= q
            if w in c:
                c[w] += p_v * p_u
            else:
                c[w] = p_v * p_u
    return filter_negl(c)


def nfoldconvolution(n, a, q):
    """Computes n-fold convolution of a distribution in log n time.

    Args:
      n: A number of iterations.
      a: A distribution represented as a dictionary.
      q: If q is defined, computation is done modulo q. If q is None,
        distributions are taken over integers.

    Returns:
      A distribution represented as a dictionary.
    """
    r = {0: 1.0}
    n_bin = bin(n)[2:]  # binary representation of n
    for ch in n_bin:
        r = convolution(r, r, q)
        if ch == '1':
            r = convolution(r, a, q)
    return r

# Standard distributions


def dgauss(sigma):
    """Computes pdf of the rounded continuous Gaussian of std sigma.

    Args:
      sigma: Standard deviation.

    Returns:
      A distribution represented as a dictionary.
    """
    d = {}
    k = 0
    while True:
        # more numerically stable than cdf (which is ~1 for large k)
        p = norm.sf(k - .5, scale=sigma) - norm.sf(k + .5, scale=sigma)
        if p < NEGLIGIBLE:
            break
        d[k] = p
        d[-k] = p
        k += 1
    return d


def sym_binomial(n):
    """Computes pdf of the symmetric binomial distribution.

    Args:
      n: A number of coins. n is even.

    Returns:
      A distribution on [-n/2 ..n/2] represented as a dictionary.
    """
    if not n % 2 == 0:
        raise ValueError("n should be even")
    d = {}
    for k in xrange(n + 1):
        d[k - n / 2] = binom.pmf(k, n, .5)
    return d


def nonnegative_half(d):
    """Given a symmetric distribution d, outputs its nonnegative half.

    Assuming that d is symmetric distribution, prepares a distribution d2 so that
    d = sign * d2, where sign is random {-1, +1}.

    Args:
      d: A distribution represented as a dictionary.

    Returns:
      A distribution represented as a dictionary.
    """
    if not valid_symmetric_distr(d):
        raise ValueError(
            "The input distribution is not a valid symmetric distribution")

    d1 = {v: x for (v, x) in d.iteritems() if v >= 0}
    if 0 in d1:
        d1[0] = d1[0] / 2
    s = sum(d1.itervalues())  # must be 1 - d1[0]
    return {v: x / s for (v, x) in d1.iteritems()}  # normalize


def distr_to_str(d):
    """Prints a distribution to string.

    Args:
       d: A distribution represented as a dictionary.

    Returns:
      Distribution as a comma-separated string in a canonical order.
    """
    r = []
    for key in sorted(d):
        r.append("({}: {})".format(key, d[key]))
    return ", ".join(r)


def bits_needed_to_sample(d):
    """Counts the number of bits required to sample from the distribution.

    Args:
       d: A distribution represented as a dictionary.

    Returns:
      Number of bits required to sample from a distribution. None if undefined.
    """
    b = 1
    while b <= 64:
        if all(fmod(x, 2 ** -b) == 0. for x in d.itervalues()):
            return b
        b += 1

    return None


def bits_to_C_type(bits):
    if bits <= 8:
        return "uint8_t"
    elif bits <= 16:
        return "uint16_t"
    elif bits <= 32:
        return "uint32_t"
    elif bits <= 64:
        return "uint64_t"
    else:
        assert False, "Cannot represent thresholds using integers"
        return ""


def distribution_to_C(distr, suffix):
    """Formats CDF as a string for use in C code.

    Args:
      d: A discrete distribution as a dictionary.
      suffix: A common suffix for all constants.
    """

    n = max(distr.iterkeys()) + 1  # range is [0..n)
    ret = "const size_t CDF_LENGTH_{} = {};\n".format(suffix, n)

    b = bits_needed_to_sample(distr)

    ret += "const uint8_t CDF_BITS_{} = {}; // random bits required to sample from the"\
        " distribution (+1 for sign)\n".format(
            suffix,
            b)

    cdf = [0] * n
    s = 0  # cumulative sum

    for i in xrange(n):
        if i in distr:
            s += distr[i] * (1 << b)
        cdf[i] = s

    ret += "const {} CDF_{}[{}] = {{".format(bits_to_C_type(b + 1), suffix, n)
    ret += ", ".join(["{}".format(int(x) - 1) for x in cdf])
    ret += "}}; // out of [0, {}]".format(2**b - 1)
    return ret
