"""Support functions for computing Renyi divergence.

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

from math import log
import warnings


def renyi(p1, p2, a):
    """Computes the Renyi divergence between two distributions.

    Args:
      p1, p2: Discrete distributions represented as dictionaries.
      a: Order of the Renyi divergence, a <> 1, a can be infinity

    Returns:
      Renyi divergence D_a(p1 || p2). Can be infinity.
    """
    if any([v not in p2 for v in p1 if p1[v] > 0]):
        return float('inf')

    if a == float('inf'):
        return max([p1[v] / p2[v] for v in p1 if p1[v] > 0])

    with warnings.catch_warnings():
        # Suppressing the overflow warning. In case of an overflow the divergence is
        # going to be infinity anyway.
        warnings.simplefilter("ignore", RuntimeWarning)
        # more numerically stable than p1[v]**a / p1[v]**(a-1)
        S = sum([p1[v] * (p1[v] / p2[v]) ** (a - 1) for v in p1 if p1[v] > 0])
    return S ** (1. / (a - 1.))


def renyi_bound(ln_pr, ln_r, a):
    """Computes an upper bound on probability via Renyi divergence.

    Computes an upper bound on probability of a distribution under p1 given a
    bound on its probability under p2 and a bound on the Renyi divergence
    between p1 and p2.

    Args:
      ln_pr: Natural logarithm of the probability of an event under p2.
      ln_r: Natural logarithm of the Renyi divergence D_a(p1 || p2).
      a: The Renyi divergence order.

    Returns:
      Natural logarithm on an upper bound on the probability of the same event under p1.
    """
    if a == float('inf'):
        return ln_pr + ln_r
    else:
        return (ln_pr + ln_r) * ((a - 1.) / a)


def opt_renyi_bound(ln_pr, p1, p2, n):
    """Finds the optimal Renyi order and the corresponding bound.

    Given two distributions p1 and p2, and an event of probability exp(ln_pr)
    under p2^n (n-fold direct product), finds the optimal Renyi order so that
    the upper bound via Renyi divergence on the event's probability under p1^n
    is minimized.

    Args:
      ln_pr: Natural logarithm of the event's probability under p2^n.
      p1, p2: Two discrete distributions represented as dictionaries.
      n: Number of dimensions.

    Returns:
      A pair of floats - an optimal Renyi order and the natural logarithm on an
      upper bound on probability under p1^n.
    """

    minbound = 0
    best_a = 1.0
    a = 1.0
    while a < 500:
        a *= 1.1
        ln_r = log(renyi(p1, p2, a)) * n
        bound = renyi_bound(ln_pr, ln_r, a)  # bound on ln(prob under p1^n)
        if bound < minbound:
            minbound = bound
            best_a = a
    return best_a, minbound
