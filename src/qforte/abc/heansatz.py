"""
HEAnsatz base classes
====================================
The mixin classes inherited by any algorithm that uses a parameterized
ansatz. Member functions should be minimal and aim only to implement
the ansatz circut and potential supporting utility functions.
"""

import qforte as qf

class HEAnsatz:
    """A mixin class for implementing the hardware-efficient ansatz, to be inherited by a
    concrete class HEAnsatz+algorithm class.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class TiledUPS(HEAnsatz):
    """A class for implementing the Tiled Unitary Product State (tUPS) ansatz.

        Attributes
        ----------
        nlayer: int >= 1
            Number of layers

        option: [None, "oo", "pp"]
            None:  tUPS
            "oo": Orbital Optimization oo-tUPS
            "pp": Perfect Pairing pp-tUPS
    """
    def __init__(self, *args, nlayer, option=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._nlayer = nlayer
        self._option = option

    def fill_layers(self, reference):
        """This function populates an operator pool with SQOperator objects."""

        self._pool_obj = qf.SQOpPool()
        self._pool_obj.set_orb_spaces(reference)
        
        for i in range(self._nlayer):
            if self._option == "pp":
                self._pool_obj.fill_pool("pp")
            else:
                self._pool_obj.fill_pool("tUPS")

        if self._option == "oo":
            self._pool_obj.fill_pool("oo")

    # def ansatz_circuit(self, amplitudes=None):
    #     """This function returns the Circuit object built
    #     from the appropriate amplitudes.

    #     Parameters
    #     ----------
    #     amplitudes : list
    #         A list of parameters that define the variational degrees of freedom in
    #         the state preparation circuit Uvqc. This is needed for the scipy minimizer.
    #     """
    #     temp_pool = qf.SQOpPool()
    #     tamps = self._tamps if amplitudes is None else amplitudes

    #     for tamp, top in zip(tamps, self._tops):
    #         temp_pool.add(tamp, self._pool_obj[top][1])

    #     U = qf.Circuit()
    #     for tamp, sq_op in temp_pool:
    #         if len(sq_op.terms()) > 2:
    #             U.add(
    #                 compact_excitation_circuit(
    #                     tamp * sq_op.terms()[2][0],
    #                     sq_op.terms()[2][1],
    #                     sq_op.terms()[2][2],
    #                 )
    #             )
    #             U.add(
    #                 compact_excitation_circuit(
    #                     tamp * sq_op.terms()[3][0],
    #                     sq_op.terms()[3][1],
    #                     sq_op.terms()[3][2],
    #                 )
    #             )
    #         else:
    #             U.add(
    #                 compact_excitation_circuit(
    #                     tamp * sq_op.terms()[1][0],
    #                     sq_op.terms()[1][1],
    #                     sq_op.terms()[1][2],
    #                 )
    #             )
        
    #     return U