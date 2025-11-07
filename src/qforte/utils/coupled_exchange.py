"""
Utility functions for handeling exponentials of QubitOperators
"""

import qforte
import numpy as np

def generate_ceo_tags(self):
    pool_size = len(self._pool_obj.terms())
    self._ceo_qubits = []
    for i in range(pool_size):
        if self._pool_type == "CEO":
            temp_sqop = self._pool_obj.terms()[i][1]
            if len(temp_sqop.terms()) == 4: # CEO operators
                temp_qop = temp_sqop.jw_transform(True)

                temp_seq = []
                for idgate in range(temp_qop.terms()[0][1].size()):
                    temp_seq.append(temp_qop.terms()[0][1].gate(idgate).target())
                seq_set = set(temp_seq)
                q3 = -1
                q4 = -1

                for qt in temp_qop.terms():
                    qt_coeff = qt[0]
                    qt_circ = qt[1]
                    if qt_coeff.imag > 0:
                        num_x = 0
                        
                        for idgate in range(qt_circ.size()):
                            if qt_circ.gate(idgate).gate_id() == 'X':
                                num_x += 1
                        
                        if num_x == 3:
                            for idgate in range(qt_circ.size()):
                                if qt_circ.gate(idgate).gate_id() == 'Y':
                                    q4 = qt_circ.gate(idgate).target()
                                    seq_set.remove(q4)
                        elif num_x == 1:
                            for idgate in range(qt_circ.size()):
                                if qt_circ.gate(idgate).gate_id() == 'X':
                                    q3 = qt_circ.gate(idgate).target()
                                    seq_set.remove(q3)

                assert q3 >= 0
                assert q4 >= 0

                seq = []
                for q in temp_seq:
                    if q in seq_set:
                        seq.append(q)
                seq.append(q3)
                seq.append(q4)

                assert len(seq) == 4
                self._ceo_qubits.append(seq)
            else:
                self._ceo_qubits.append(None)
        else:
            self._ceo_qubits.append(None)

                                
                            