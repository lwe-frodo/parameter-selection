"""Searches parameter space for LWE-based key exchange protocols.

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

from math import log, sqrt, fmod, ceil, isinf
from pqsec import optimize_attack, svp_classical, svp_quantum, svp_plausible, primal_cost, dual_cost
from failure_prob import heuristic_failure_prob_upper, noise_failure_prob
from approx_distr import approximate_dgauss
from discrete_distr import pdf_product, dgauss, sym_binomial, distr_to_str, nonnegative_half, bits_needed_to_sample
from renyi import opt_renyi_bound


def estimate_cost(q, n, max_m, s, attacks):
    """ Finds attack costs (classical, quantum, plausible quantum) for an LWE instance.

    Args:
      q: LWE modulus.
      n: LWE dimension.
      max_m: maximum number of samples
      s: standard deviation of the error distribution
      attacks: A vector specifying types of attacks to consider (svp_classical,
          svp_quantum, svp_plausible)

    Returns:
      A triple of log_2 of the costs of three attacks: classical, quantum, and plausible quantum.
      Infinity if analysis is not run.
    """
    r = [float('inf')] * 3
    for cost in [primal_cost, dual_cost]:
        for i, attack in enumerate(attacks):
            _, _, c = optimize_attack(q, n, max_m, s, cost, attack)
            r[i] = min(r[i], c)
    return r

# Minimally acceptable probability of key agreement failure
_PROB_FAILURE_CUTOFF = 2 ** -32


def minimize_bandwidth(
        classical_lb,
        quantum_lb,
        plausible_lb,
        ubits,
        binomials_only=False,
        reduce_to_LWE=True,
        agree_bits=256):
    """Searches the parameter space to minimize bandwidth. Prints to stdout.

    Args:
        classical_lb: A lower bound on classical security (None if undefined).
        quantum_lb: A lower bound on quantum security (None if undefined).
        plausible_lb: A conservative lower bound on security (None if undefined).
        ubits: The bound on the number of uniform bits required for sampling.
        binomials_only: If True, considers only binomial distributions.
        reduce_to_LWE: If True, uses reduction to rounded-Gaussian LWE.
        agree_bits: Target key length.
    """

    def check_costs(cost_cl, cost_pq, cost_pp):
        """Checks whether the costs satisfy lower bounds.

        Returns:
            True if lower bounds are satisfied.
        """
        return ((classical_lb is None or cost_cl >= classical_lb)
                and (quantum_lb is None or cost_pq >= quantum_lb)
                and (plausible_lb is None or cost_pp >= plausible_lb))

    def find_opt_distr(sigma, samples, ubits, cost_cl, cost_pq, cost_pp):
        """Finds an optimal distribution approximating rounded continuous Gaussian.

        Args:
            sigma: The standard deviation of the target (rounded) Gaussian.
            samples: The total number of samples drawn by both parties combined.
            ubits: The bound on the number of uniform bits required for sampling.
            cost_cl, cost_pq, cost_pp: Estimated costs of the rounded Gaussian.

        Returns:
            Four-tuple consisting of the distribution and the cost triplet.
        """
        cost_cl_opt, d, _ = approximate_dgauss(
            sigma, samples, cost_cl, None, ubits, quiet=True)

        sym_d = pdf_product(d, {+1: .5, -1: .5})

        dg = dgauss(sigma)

        _, cost_pq_opt = opt_renyi_bound(-cost_pq * log(2), sym_d, dg, samples)
        _, cost_pp_opt = opt_renyi_bound(-cost_pp * log(2), sym_d, dg, samples)

        return [sym_d, cost_cl_opt, -cost_pq_opt /
                log(2), -cost_pp_opt / log(2)]

    def find_binomial_cost(sigma, samples, cost_cl, cost_pq, cost_pp):
        """Estimates the cost of replacing a rounded Gaussian with a binomial.

        Args:
            sigma: The standard deviation of the Gaussian.
            samples: The total number of samples drawn by Alice and Bob.
            cost_cl, cost_pq, cost_pp: Estimated costs of the rounded Gaussian.

        Returns:
            Four-tuple consisting of the distribution and the cost triplet.
        """

        dg = dgauss(sigma)
        # The binomial is defined as B(2*z, .5) - z.
        sb = sym_binomial(2 * sigma**2)

        _, cost_cl_binomial = opt_renyi_bound(
            -cost_cl * log(2), sb, dg, samples)
        _, cost_pq_binomial = opt_renyi_bound(
            -cost_pq * log(2), sb, dg, samples)
        _, cost_pp_binomial = opt_renyi_bound(
            -cost_pp * log(2), sb, dg, samples)

        return [
            sb,
            - cost_cl_binomial / log(2),
            - cost_pq_binomial / log(2),
            - cost_pp_binomial / log(2)]

    def print_intermediate_result(
            opt_d,
            qlog,
            n,
            w,
            sigma,
            bandwidth,
            heuristic_pr_failure,
            actual_pr_failure,
            costs_gauss,
            costs_opt):

        print distr_to_str(nonnegative_half(opt_d))
        print "q = 2**{}, n = {}, w = {}, sigma^2 = {:.2f}:".format(qlog, n, w, sigma**2),
        print "bandwidth = {:.2f} KB,".format(bandwidth / (8. * 1000)),

        print ("heuristic Pr of failure = {:.1f}, "
               "actual Pr of failure = {:.1f},".format(
                   log(heuristic_pr_failure, 2),
                   log(actual_pr_failure, 2))),

        formatted_costs = ", ".join("{:.1f}".format(c)
                                    for c in costs_gauss if not isinf(c))
        print "security = [{}],".format(formatted_costs),

        if reduce_to_LWE or not binomials_only:
            formatted_costs_opt = ", ".join(
                "{:.1f}".format(c) for c in costs_opt if not isinf(c))
            print "security after reduction = [{}],".format(formatted_costs_opt),

        if binomials_only:
            bits = int(round(4 * sigma**2)) + 1
        else:
            bits = bits_needed_to_sample(opt_d)
        print "distribution on [{}, {}], {} bits to sample".format(
            min(opt_d.iterkeys()),
            max(opt_d.iterkeys()),
            bits)

        print "Parameters for print_tables.py: 'sigma': sqrt({:.2f}), 'n': {}, 'q': {}, 'B': {}, 'bits': {}, 'base': {:.0f}".format(sigma**2, n, qlog, w, bits, min(costs_gauss))

    # Main body of minimize_bandwidth starts here.

    attacks = [svp_classical]
    if quantum_lb is not None:
        attacks.append(svp_quantum)
    if plausible_lb is not None:
        attacks.append(svp_plausible)

    qlog_list = range(10, 16)  # Enumerate moduli in the range 2^10..2^15.
    # Enumerate n in the range 256..900.
    n_list = [n for n in xrange(256, 900, 16)]
    # The number of bits to extract; can range between 1 and be as large as
    # the modulus.
    w_list = range(1, 1 + max(qlog_list))
    # Feasible values for the variance.
    v_list = (1., 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 2.)

    min_bandwidth_bits = 8 * 30000  # bandwidth in bits

    for qlog in qlog_list:
        print "q = 2**{}".format(qlog)
        for n in n_list:
            prev_k = None  # number of columns already considered
            for w in w_list:
                # number of coordinates needed to achieve agree_bits
                required_coordinates = (agree_bits + w - 1) / w
                # number of columns
                k = int(ceil(sqrt(required_coordinates)))
                if k == prev_k:
                    # Nothing to do: w is larger than before but the number
                    # of columns is the same.
                    continue
                prev_k = k

                max_m = n + k
                bandwidth = k * qlog * 2 * n + required_coordinates
                if bandwidth > min_bandwidth_bits:
                    continue

                for v in v_list:
                    if binomials_only and fmod(2 * v, 1.0) != 0:
                        continue

                    sigma = sqrt(v)

                    # Quick-and-dirty check on the probability of failure
                    heuristic_pr_failure = heuristic_failure_prob_upper(
                        2 ** qlog, n, sigma, w, agree_bits)
                    if heuristic_pr_failure > 2**-20:  # very loose bound
                        continue

                    samples = 2 * k * n + required_coordinates

                    costs_gauss = estimate_cost(
                        2 ** qlog, n, max_m, sigma, attacks)

                    # Checking upper bounds on costs, before trying to find an optimal
                    # approximation.
                    if not check_costs(*costs_gauss):
                        continue

                    if binomials_only:
                        if reduce_to_LWE:
                            [opt_d, cost_cl_opt, cost_pq_opt, cost_pp_opt] = find_binomial_cost(
                                int(round(2 * v)), samples, *costs_gauss)
                        else:
                            # no reduction
                            [cost_cl_opt, cost_pq_opt, cost_pp_opt] = costs_gauss
                            opt_d = sym_binomial(2 * int(round(2 * v)))
                    else:
                        [opt_d, cost_cl_opt, cost_pq_opt, cost_pp_opt] = find_opt_distr(
                            sigma, samples, ubits, *costs_gauss)

                    costs_opt = [cost_cl_opt, cost_pq_opt, cost_pp_opt]

                    if reduce_to_LWE and not check_costs(*costs_opt):
                        continue

                    actual_pr_failure = noise_failure_prob(
                        opt_d, 2 ** qlog, n, w, required_coordinates)
                    if actual_pr_failure > _PROB_FAILURE_CUTOFF:
                        continue

                    min_bandwidth_bits = bandwidth
                    print_intermediate_result(
                        opt_d,
                        qlog,
                        n,
                        w,
                        sigma,
                        bandwidth,
                        heuristic_pr_failure,
                        actual_pr_failure,
                        costs_gauss,
                        costs_opt)

    return


def main():
    binomials_only = False
    reduce_to_LWE = True

    print "### Challenge ###"
    minimize_bandwidth(64, None, None, 8, binomials_only, reduce_to_LWE, 64)
    print

    print "### Classical ###"
    minimize_bandwidth(128, None, None, 12, binomials_only, reduce_to_LWE, 128)
    print

    print "### Recommended ###"
    minimize_bandwidth(128, 128, None, 12, binomials_only, reduce_to_LWE, 256)
    print

    print "### Paranoid ###"
    minimize_bandwidth(128, 128, 128, 16, binomials_only, reduce_to_LWE, 256)
    print

if __name__ == "__main__":
    main()
