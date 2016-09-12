"""Generation of LaTeX tables for the paper and CDFs for the C code.

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
from math import sqrt, log
from discrete_distr import dgauss, bits_needed_to_sample, nonnegative_half, pdf_product, distribution_to_C
from renyi import renyi, opt_renyi_bound
from approx_distr import approximate_dgauss
from pqsec import optimize_attack, svp_classical, svp_quantum, svp_plausible, primal_cost, dual_cost
from collections import namedtuple
from failure_prob import noise_failure_prob


def distribution_to_TeX(p):
    """Formats distribution for use in TeX file.

    Args:
      p: A dictionary with entries for 'distr', 'sigma', 'a', 'D'

    Returns:
      LaTeX string.
    """

    distr, sigma, a, distr_name = p['distr'], p['sigma'], p['a'], p['D']

    n = max(distr.iterkeys()) + 1  # range is [0..n)
    b = bits_needed_to_sample(distr)
    ret = "${}$ & {} & {:.2f} & ".format(distr_name, b + 1, sigma**2)

    for i in sorted(distr.iterkeys()):
        if i >= 0:
            p = int(round(distr[i] * 2**b))
            if i == 0:
                p *= 2
            ret += r"\!\!{}\!\!&".format(p)

    for i in xrange(max(distr.iterkeys()), 6):
        ret += "&"

    divergence = renyi(distr, nonnegative_half(dgauss(sigma)), a)
    ret += " {:.1f} & {:.7f}".format(a, divergence)
    return ret + r" \\"


def security_to_TeX(p, nbar, print_sec=True):
    """Formats security estimates for use in TeX file.

    Args:
      p: A dictionary with entries for 'name', 'n', 'q', 'distr', 'sigma'
      nbar: Number of columns in exchanged matrices.
      print_sec: If true, output will include security estimates

    Returns:
      LaTeX string.
    """
    name, n, qlog, d, sigma = p['name'], p['n'], p['q'], p['distr'], p['sigma']

    samples = 2 * n * nbar + nbar**2
    q = 2**qlog
    ret = ""

    ret += r"\multirow{2}{*}{" + name.capitalize() + "} "

    for cost in [primal_cost, dual_cost]:
        m_pc, b_pc, cost_pc = optimize_attack(
            q, n, samples, sigma, cost, svp_classical, verbose=False)
        m_pq, b_pq, cost_pq = optimize_attack(
            q, n, samples, sigma, cost, svp_quantum, verbose=False)
        m_pp, b_pp, cost_pp = optimize_attack(
            q, n, samples, sigma, cost, svp_plausible, verbose=False)

        if cost == primal_cost:
            ret += "& Primal & "
        else:
            ret += "& Dual & "

        ret += "{} & {} &".format(m_pc, b_pc)

        if print_sec:
            sym_d = pdf_product(d, {+1: .5, -1: .5})
            dg = dgauss(sigma)

            _, cost_pc_reduced = opt_renyi_bound(-cost_pc *
                                                 log(2), sym_d, dg, samples)
            _, cost_pq_reduced = opt_renyi_bound(-cost_pq *
                                                 log(2), sym_d, dg, samples)
            _, cost_pp_reduced = opt_renyi_bound(-cost_pp *
                                                 log(2), sym_d, dg, samples)

            ret += "{} & {} & {} & {} & {} & {} \\\\".format(
                int(cost_pc),
                int(cost_pq),
                int(cost_pp),
                int(-cost_pc_reduced / log(2)),
                int(-cost_pq_reduced / log(2)),
                int(-cost_pp_reduced / log(2)))  # always round down
        else:
            ret += "-- & -- & -- & -- & -- & -- \\\\"

        ret += "\n"
    return ret


def parameters_to_TeX(p, nbar):
    """Formats parameters for use in TeX file.

    Args:
      p: A dictionary with entries for 'name', 'n', 'q', 'D', 'B', 'distr'
      nbar: Number of columns in exchanged matrices.

    Returns:
      LaTeX string.
    """
    name, n, q, dname, B, distr = p['name'], p[
        'n'], p['q'], p['D'], p['B'], p['distr']

    s = name.capitalize()
    s += " & \!\! {} \!\!".format(n)
    s += " & \!\! $2^{{{}}}$ \!\!".format(q)
    s += " & \!\! ${}$ \!\!".format(dname)

    agreed_bits = B * nbar * nbar
    s += " & ${}\cdot {}^2 = {}$".format(B, nbar, agreed_bits)
    if agreed_bits < 100:  # pad to equal length
        s += "\hphantom{0}"

    sym_distr = pdf_product(distr, {+1: .5, -1: .5})
    failure_prob = noise_failure_prob(sym_distr, 2 ** q, n, B, agreed_bits)
    s += " & $2^{{{:.1f}}}$".format(log(failure_prob, 2))

    bandwidth = nbar * q * 2 * n + nbar * nbar  # bandwidth in bits
    s += " & {:.2f} KB".format(bandwidth / 8000.)

    s += r" \\"

    return s


def main():
    parameters = [{'name': 'challenge',
                   'D': 'D_1', 'sigma': sqrt(1.25), 'n': 352, 'q': 11, 'B': 1, 'bits': 8, 'base': 85},
                  {'name': 'classical',
                   'D': 'D_2', 'sigma': sqrt(1.00), 'n': 592, 'q': 12, 'B': 2, 'bits': 12, 'base': 137},
                  {'name': 'recommended',
                   'D': 'D_3', 'sigma': sqrt(1.75), 'n': 752, 'q': 15, 'B': 4, 'bits': 12, 'base': 138},
                  {'name': 'paranoid',
                   'D': 'D_4', 'sigma': sqrt(1.75), 'n': 864, 'q': 15, 'B': 4, 'bits': 16, 'base': 129},
                  ]

    for p in parameters:
        nbar, mbar = 8, 8
        samples = (nbar + mbar) * p['n'] + nbar * mbar
        _, p['distr'], p['a'] = approximate_dgauss(
            p['sigma'], samples, p['base'], None, p['bits'], quiet=True)

    print "### C Code ###"
    for p in parameters:
        suffix = p['D'].replace('_', '')
        print distribution_to_C(p['distr'], suffix)
        print

    print "### TABLE 1 ###"
    for p in parameters:
        print distribution_to_TeX(p)
    print

    print "### TABLE 2 ###"
    for p in parameters:
        print parameters_to_TeX(p, nbar=8)
    print

    print "### TABLE 3 ###"
    for p in parameters:
        print security_to_TeX(p, nbar=8, print_sec=p['name'] != 'challenge'),
        if p['name'] == 'paranoid':
            print r"\bottomrule"
        else:
            print r"\midrule"

if __name__ == "__main__":
    main()
