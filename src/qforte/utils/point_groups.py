"""
point_groups.py
====================================================
This module contains information about the character
tables of the point groups supported by QForte and
a function that finds the irrep of a given SQOp.
Currently, this module pertains only to molecular
objects with "build_type = 'psi4'".
"""

def char_table(point_group):
    """
    Function that prints the character table of a chosen point group.

    Parameters
    ----------
    point_group: list of two elements; point_group[0] is a string holding the name of the point group
                                       point_group[1] is a list that holds the irreps of the group in the Cotton ordering

    """

    if not isinstance(point_group, (list, tuple)):
        raise TypeError("""{0} is not a list.
                This function takes arguments of the form:\n
                ['c2v', ['A1', 'A2', 'B1', 'B2']]\n
                using the so-called Cotton ordering of the irreps.""".format(type(point_group)))

    group = point_group[0].lower()
    groups = ['c1', 'c2', 'ci', 'cs', 'd2', 'c2h', 'c2v', 'd2h']

    if group not in groups:
        raise ValueError('The given point group is not supported. Choose one of:\n{0}'.format(groups))

    irreps = point_group[1]

    print('==========>', group.capitalize(), '<==========')
    print('')
    print('      ', end = '')

    for irrep in irreps:
        print(irrep.ljust(3), ' ', end = '')

    print('\n')

    for idx1, irrep1 in enumerate(irreps):
        print(irrep1.ljust(3), '  ', end = '')
        for idx2, irrep2 in enumerate(irreps):
            print(irreps[idx1 ^ idx2].ljust(3), ' ', end = '')
        print()
