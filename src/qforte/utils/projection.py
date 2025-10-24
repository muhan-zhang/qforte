"""
Functions for constructing the projection operator for a Jordon-Wigner transformed Hamiltonian.
"""

import qforte as qf
import math
from numpy.polynomial.legendre import leggauss
from scipy.integrate import lebedev_rule

def construct_proj(
    is_sz_eig=None,
    target_s=None,
    target_ms=None,
    target_n=None,
    target_irrep=None, 
    _nqb=None,
    _ref=None,
    _point_group=None,
    _orb_irreps_to_int=None,
    ntrapz=None, 
    ngl=None,
    nlebedev=None,
):
    """
    Attributes
    ----------
    is_sz_eig: bool
        True if the ansatz is always an eigenstate of Sz.

    target_s: int or float
        Total spin quantum number (integer or half-integer).
        P(s; ms) = \int_\Omega d\Omega exp(-i \alpha Sz) exp(-i \beta Sy) exp(-i \gamma Sz)
        where \Omega = (\alpha, \beta, \gamma) is the Euler angle.
        If is_sz_eig is True, P(s; ms) can be simplified as \beta -only integration using 
        Gauss-legendre quadrature. Otherwise, it would be a Lebedev quadrature for (\alpha, \beta)
        combined with a trapezoidal quadrature for \gamma by default. (The user can set ngl value
        with a nlebedev=None to translate to a Trap + GL + Trap quadrature instead of Leb + Trap)

    target_ms: int or float
        Spin z (azimuthal projection) quantum number (integer or half integer).
        Projected within total spin projection by default.
        Trapezoidal quadrature by default if separated.

    target_n: int
        Particle number.
        Trapezoidal quadrature by default.
    
    target_irrep: int
        Irreducible representation mapped to an integer using Cotton ordering.
        Only Abelian groups are supported.
        Discrete summation.
    
    ntrapz: int
        Number of trapezoidal sampling points.
        By default, max(N_e, N_o - N_e) with N_e number of electrons and N_o number of spinorbitals.

    ngl: int
        Number of Gauss-Legendre sampling points.
        By default, 2.
    
    nlebedev: int
        Lebedev quadrature order.
        Only selected integers are allowed: {
            3, 5, 7, 9, 11, 13, 15, 17,
            19, 21, 23, 25, 27, 29, 31, 35,
            41, 47, 53, 59, 65, 71, 77, 83,
            89, 95, 101, 107, 113, 119, 125, 131
        }
        The corresponding degrees (number of sampling points) are: {
            6, 14, 26, 38, 50, 74, 86, 110,
            146, 170, 194, 230, 266, 302, 350, 434,
            590, 770, 974, 1202, 1454, 1730, 2030, 2354,
            2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
        }
    """

    grad_meas_coeff = 1
    n_cnot_proj = 0
    projectors = []
    projection = {}

    if (target_n is not None) or (not is_sz_eig):
        if ntrapz is None:
            ntrapz = int(max(sum(_ref), _nqb - sum(_ref))) # * 2 # See paper
        intvl = 2 * math.pi / ntrapz

    # NOTE: number projector
    if target_n is not None:
        projn = qf.QubitOperator()
        for idn in range(ntrapz):
            phi = intvl * idn
            wg = complex(math.cos(phi * (_nqb / 2 - target_n)), 
                         math.sin(phi * (_nqb / 2 - target_n))) / ntrapz
            Ug = qf.Circuit()
            for ig in range(0, _nqb, 2):
                Ug.add(qf.gate("Rz", ig, ig, phi))
                Ug.add(qf.gate("Rz", ig + 1, ig + 1, phi))
            projn.add(wg, Ug)
        projectors.append(projn)
        grad_meas_coeff *= ntrapz
        n_cnot_proj += 2 * _nqb
    
    # NOTE: spatial projector
    if target_irrep is not None:
        projirrep = qf.QubitOperator()
        irrep_ops, op_counts = qf.symop_system(
            _point_group[0],  # string specifying the point group
            _orb_irreps_to_int  # list of irreps for this system
        )
        op_sum = sum(op_counts)
        max_nczgates = -1
        for gm, op_count in zip(irrep_ops, op_counts):
            Ug = qf.Circuit()
            wg = op_count / op_sum
            nczgates = 0
            for ig in range(0, int(_nqb / 2)):
                if gm[ig] == -1:
                    Ug.add(qf.gate("Z", ig * 2))
                    Ug.add(qf.gate("Z", ig * 2 + 1))
                    nczgates += 2
            projirrep.add(wg, Ug)
            max_nczgates = max(max_nczgates, nczgates)
        projectors.append(projirrep)
        grad_meas_coeff *= len(irrep_ops)
        n_cnot_proj += 2 * max_nczgates

    # NOTE: the following code block is to generate spin projection operator (projector).
    #       However, the projector generated is only for specific evaluations with a caveat
    #       to project an arbitrary state to a non-Sz-symmetry-adapted space.
    #       This a priori simplifies calculations of energy expval and energy gradients

    if is_sz_eig is not None:

        projs2 = qf.QubitOperator()

        if is_sz_eig:
            method = "gauss-legendre"
            if ngl is None:
                ngl = 2 # default value
        else:
            if ngl is None and nlebedev is None:
                nlebedev = 3 # lebedev grid (default value)
                ntrapzsz = ntrapz
                method = "lebedev"
            elif ngl:
                ntrapzsz = ntrapz
                method = "gauss-legendre"
            else:
                ntrapzsz = ntrapz
                method = "lebedev"

        small_ds = []

        if not is_sz_eig:
            projsz = qf.QubitOperator()
            for ida in range(ntrapzsz):
                alpha = intvl * ida
                wg = complex(math.cos(alpha * target_ms), math.sin(alpha * target_ms)) / ntrapzsz
                Ug = qf.Circuit()
                # NOTE: exp(-i gamma Sz)
                for ia in range(0, _nqb, 2):
                    Ug.add(qf.gate("Rz", ia, ia, -alpha / 2))
                    Ug.add(qf.gate("Rz", ia + 1, ia + 1, alpha / 2))
                projsz.add(wg, Ug)

        if method == "gauss-legendre":
            gl_quad_points, gl_quad_weights = leggauss(ngl)
            projection.update(
                {
                    "gl_quad_points": gl_quad_points, 
                    "gl_quad_weights": gl_quad_weights,
                    "min_n_cnot_proj": 2 * _nqb
                }
            )
            betas = []
            for pt, wG in zip(
                projection.get("gl_quad_points"),
                projection.get("gl_quad_weights"),
            ):
                beta = math.pi - math.acos(pt)  # Gauss-Legendre quad

                small_d = wigner_d(beta, target_s, target_ms)

                # NOTE: weight of each Ug
                #       defined by Euler angle integration, wigner small d, quadrature
                wg = (target_s + 0.5) * small_d * wG 

                # NOTE: Unitary unit for summation
                Ug = qf.Circuit()
                for ib in range(0, _nqb, 2):
                    Ug.add(
                        qf.compact_excitation_circuit(
                            beta / 2.0, [ib + 1], [ib], qubit_excitations=False
                        )
                    )

                betas.append(beta)
                small_ds.append(small_d)
                projs2.add_term(wg, Ug)
                projection.update(
                    {
                        "betas": betas,
                        "small_ds": small_ds,
                    }
                )

        elif method == "lebedev":
            cartcoords, wl = lebedev_rule(nlebedev) # wl has been normalized to 4pi
            nlebsamp = len(wl)
            betas = []
            gammas = []
            for i in range(len(wl)):
                beta = math.acos(cartcoords[2, i])
                small_d = wigner_d(beta, target_s, target_ms)
                xyr = math.sqrt(cartcoords[0, i]**2 + cartcoords[1, i]**2)
                gamma = math.pi * (1 - cartcoords[2, i]) / 2 if xyr == 0 else math.copysign(
                    math.acos(cartcoords[0, i] / xyr),
                    cartcoords[1, i]
                ) + math.pi * (1 - math.copysign(1, cartcoords[1, i]))
                wg = (target_s + 0.5) * small_d * wl[i] / (2 * math.pi) \
                   * complex(math.cos(gamma * target_ms), math.sin(gamma * target_ms)) 

                Ug = qf.Circuit()
                # NOTE: exp(-i gamma Sz)
                for ig in range(0, _nqb, 2):
                    Ug.add(qf.gate("Rz", ig, ig, -gamma / 2))
                    Ug.add(qf.gate("Rz", ig + 1, ig + 1, gamma / 2))
                for ib in range(0, _nqb, 2):
                    Ug.add(
                        qf.compact_excitation_circuit(
                            beta / 2.0, [ib + 1], [ib], qubit_excitations=False
                        )
                    )
                betas.append(beta)
                gammas.append(gamma)
                small_ds.append(small_d)
                projs2.add_term(wg, Ug)
                projection.update(
                    {
                        "betas": betas,
                        "alphas": gammas,
                        "small_ds": small_ds,
                    }
                )

        if is_sz_eig:
            n_cnot_proj += int(
                int(_nqb / 2) * 8
            )  # 8 by controlled-Rz with one ancilla qubit
            grad_meas_coeff *= ngl
            projectors.append(projs2)
        else:
            n_cnot_proj += int(int(_nqb / 2) * 8) + 4 * _nqb
            if method == "gauss-legendre":
                grad_meas_coeff *= ngl
                grad_meas_coeff *= ntrapz**2
                projectors.append(projsz)
                projectors.append(projs2)
                projectors.append(projsz)
            elif method == "lebedev":
                grad_meas_coeff *= nlebsamp
                grad_meas_coeff *= ntrapz
                projectors.append(projs2)
                projectors.append(projsz)
    
    projection.update(
        {
            "projector": projectors,
            "n_cnot_proj": n_cnot_proj,
            "grad_meas_coeff": grad_meas_coeff,
        }
    )

    return projection

def wigner_d(beta, target_s, target_ms):
    """Function to calculate Wigner-(small)d value.
    """
    # NOTE: wigner small d calculation -- real number
    jmax = min(target_s + target_ms, target_s - target_ms)
    small_d = 0.0
    for j in range(int(jmax) + 1):
        small_d += (
            ((-1) ** j)
            * (math.cos(beta / 2.0) ** (2.0 * (target_s - j)))
            * (math.sin(beta / 2.0) ** (2.0 * j))
            / math.factorial(int(target_s + target_ms - j))
            / (math.factorial(j) ** 2)
            / math.factorial(int(target_s - target_ms - j))
        )
    small_d *= math.factorial(int(target_s + target_ms)) * math.factorial(
        int(target_s - target_ms)
    )
    return small_d