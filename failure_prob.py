"""Computing failure probabilities (heuristic and exact) of the key agreement protocol.

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


from math import sqrt, exp, log
from discrete_distr import pdf_product, std_modulo, nfoldconvolution, convolution, sym_binomial, dgauss


def heuristic_failure_prob_upper(q, n, sigma, w, reclen):
    """Computes an upper bound on the heuristic probability of failure.

    Args:
      q: Modulus.
      n: Vector length.
      sigma: Standard deviation of the noise distribution.
      w: Number of extracted bits per coordinate.
      reclen: Number of coordinates to be reconciled.

    Returns:
      The upper bound on the heuristic probability of failure.
    """
    std = sqrt(2 * n * sigma ** 4 + sigma ** 2)
    cuttail = q / (std * 2 ** (w + 2))
    r1r2 = (reclen + w - 1) / w
    return 2 * r1r2 * exp(-cuttail ** 2 / 2) / log(2)


def noise_failure_prob(noise, q, n, w, reclen):
    """Computes the failure probability of key agreement for given noise distribution.

    Args:
      noise: Noise distribution given as a dictionary.
      q: Modulus.
      n: Vector length.
      w: Number of extracted bits per coordinate.
      reclen: Number of coordinates to be reconciled.

    Returns:
      The probability of failure.
    """

    def pr_rec_failure(x):
        """ Probability of failing to recover from an error of magnitude x:
            0% if x <= B/4
            100% if x >= 3B/4
            and linear in-between
        """
        x = min(x, q - x)
        B = q / (2 ** w)
        if x < B / 4:
            return 0
        elif x > 3 * B / 4:
            return 1.
        else:
            return (x - B / 4.) / (2. * B / 4.)

    noise_sqr = pdf_product(noise, noise, q)

    v = nfoldconvolution(2 * n, noise_sqr, q)

    v = convolution(v, noise, q)  # v = 2n * (noise^2) + noise

    exact_pr = {x: p * pr_rec_failure(x) for (x, p) in v.iteritems()}

    failure_pr = reclen * sum(exact_pr.itervalues())

    return failure_pr


def print_failure_props(noise, qlog, n, w, agree_bits=256):
    reclen = (agree_bits + w - 1) / w
    s = std_modulo(noise, 2 ** qlog)

    print "p = 2^{}, n = {}, B = {}, |key| = {},".format(qlog, n, w, agree_bits),
    print "exact pr of failure = 2^{:.2f}".format(log(noise_failure_prob(noise, 2 ** qlog, n, w, reclen), 2)),
    print "(heuristic = 2^{:.2f})".format(log(heuristic_failure_prob_loose(2 ** qlog, n, s, w, reclen), 2))
    print


def main():
    d1 = {0: 44 / 128., 1: 61 / 128., 2: 20 / 128., 3: 3. / 128}
    sym_d1 = pdf_product(d1, {+1: .5, -1: .5})
    # Parameters with deliberately large probability of failure
    print "Distribution = D1:",
    print_failure_props(sym_d1, 10, 320, 1, 64)

    dg175 = dgauss(sqrt(1.75))
    print "Distribution = rounded Gaussian with sigma^2 = 1.75:",
    print_failure_props(dg175, 15, 752, 4, 256)

    d3 = {
        0: 603 / 2048.,
        1: 919 / 2048.,
        2: 406 / 2048.,
        3: 104 / 2048.,
        4: 15 / 2048.,
        5: 1 / 2048.}
    sym_d3 = pdf_product(d3, {+1: .5, -1: .5})
    print "Recommended parameters:",
    print_failure_props(sym_d3, 15, 752, 4, 256)


if __name__ == "__main__":
    main()
