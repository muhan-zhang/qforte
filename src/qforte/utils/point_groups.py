"""
point_groups.py
====================================================
This module contains information about the character
tables of the point groups supported by QForte and
a function that finds the irrep of a given SQOp.
Currently, this module pertains only to molecular
objects with "build_type = 'psi4'".
"""

import numpy as np


def irreps_of_point_groups(point_group):
    """
    Function that returns the irreps of a given point group,
    using the Cotton ordering.
    """

    point_group = point_group.lower()
    point_group_to_irreps = {
        "c1": ["A"],
        "c2": ["A", "B"],
        "ci": ["Ag", "Au"],
        "cs": ["Ap", "App"],
        "d2": ["A", "B1", "B2", "B3"],
        "c2h": ["Ag", "Bg", "Au", "Bu"],
        "c2v": ["A1", "A2", "B1", "B2"],
        "d2h": ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"],
    }

    try:
        return point_group_to_irreps[point_group]
    except KeyError:
        raise ValueError(
            f'The given point group is not supported. Choose one of:\n{", ".join(point_group_to_irreps.keys())}'
        )


def char_table(point_group):
    """
    Function that prints the character table of a chosen point group.

    Parameters
    ----------
    point_group: list of two elements; point_group[0] is a string holding the name of the point group
                                       point_group[1] is a list that holds the irreps of the group in the Cotton ordering

    """

    if not isinstance(point_group, (list, tuple)):
        raise TypeError(
            """{0} is not a list.
                This function takes arguments of the form:\n
                ['c2v', ['A1', 'A2', 'B1', 'B2']]\n
                using the so-called Cotton ordering of the irreps.""".format(
                type(point_group)
            )
        )

    group = point_group[0].lower()
    groups = ["c1", "c2", "ci", "cs", "d2", "c2h", "c2v", "d2h"]

    if group not in groups:
        raise ValueError(
            "The given point group is not supported. Choose one of:\n{0}".format(groups)
        )

    irreps = point_group[1]

    print("==========>", group.capitalize(), "<==========")
    print("")
    print("      ", end="")

    for irrep in irreps:
        print(irrep.ljust(3), " ", end="")

    print("\n")

    for idx1, irrep1 in enumerate(irreps):
        print(irrep1.ljust(3), "  ", end="")
        for idx2, irrep2 in enumerate(irreps):
            print(irreps[idx1 ^ idx2].ljust(3), " ", end="")
        print()


def symmetry_operations(point_group):
    """
    Function that prints the character table of a chosen point group corresponding to symmetry operations.

    Parameters
    ----------
    point_group: string 
        holding the name of the point group
    """
    group = point_group.lower()
    groups = ["c1", "c2", "ci", "cs", "d2", "c2h", "c2v", "d2h"]

    if group in {"c1"}:
        sym_op = [[1]]
    elif group in {
        "c2",
        "ci",
        "cs"
    }:
        sym_op = [
            [1, 1],
            [1, -1]
        ]
    elif group in {
        "d2",
        "c2h",
        "c2v"
    }:
        sym_op = [
            [1, 1, 1, 1],
            [1, 1, -1, -1],
            [1, -1, 1, -1],
            [1, -1, -1, 1]
        ]
    elif group in {"d2h"}: # d2h
        sym_op = [
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, -1, -1, 1, 1, -1, -1],
            [1, -1, 1, -1, 1, -1, 1, -1],
            [1, -1, -1, 1, 1, -1, -1, 1],
            [1, 1, 1, 1, -1, -1, -1, -1],
            [1, 1, -1, -1, -1, -1, 1, 1],
            [1, -1, 1, -1, -1, 1, -1, 1],
            [1, -1, -1, 1, -1, 1, 1, -1],
        ]
    else:
        raise ValueError(
            "The given point group is not supported. Choose one of:\n{0}".format(groups)
        )

    return sym_op

def symop_system(point_group, orb_irreps_to_int):
    """
    Function that prints the character table of a chosen point group corresponding to symmetry operations.

    Parameters
    ----------
    point_group: string 
        holding the name of the point group

    orb_irreps_to_int: list[int]
    """
    sym_op = symmetry_operations(point_group)
    symop_sys = [[operation[irrep] for irrep in orb_irreps_to_int] for operation in sym_op]
    unique_ops, op_count = np.unique(symop_sys, axis=0, return_counts=True)
    return unique_ops, op_count