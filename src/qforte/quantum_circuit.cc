#include <algorithm>

#include "helpers.h"
#include "quantum_gate.h"
#include "quantum_circuit.h"

void QuantumCircuit::set_parameters(const std::vector<double>& params) {
    // need a loop over only gates in state preparation circuit that
    // have a parameter dependance (if gate_id == Rx, Ry, or Rz)
    // TODO: make a indexing funciton using a map (Nick)
    size_t param_idx = 0;
    for (auto& gate : gates_) {
        std::string gate_id = gate.gate_id();
        if (gate_id == "Rz") {
            size_t target_qubit = gate.target();
            gate = make_gate(gate_id, target_qubit, target_qubit, params[param_idx]);
            param_idx++;
        }
    }
}

void QuantumCircuit::add_circuit(const QuantumCircuit& circ) {
    for (const auto gate : circ.gates()) {
        gates_.push_back(gate);
    }
}

QuantumCircuit QuantumCircuit::adjoint() {
    QuantumCircuit qcirc_adjoint;
    for (auto& gate : gates_) {
        qcirc_adjoint.add_gate(gate.adjoint());
    }
    std::reverse(std::begin(qcirc_adjoint.gates_), std::end(qcirc_adjoint.gates_));
    return qcirc_adjoint;
}

std::complex<double> QuantumCircuit::canonical_order() {
    /// TODO: peraphs put ths somewhere else?
    using namespace std::complex_literals;
    std::map<
        std::pair<std::string,std::string> ,
        std::pair<std::complex<double>,std::string>
        > m = {
        {std::make_pair("X", "Y"), std::make_pair( 1.0i, "Z")},
        {std::make_pair("X", "Z"), std::make_pair(-1.0i, "Y")},
        {std::make_pair("Y", "X"), std::make_pair(-1.0i, "Z")},
        {std::make_pair("Y", "Z"), std::make_pair( 1.0i, "X")},
        {std::make_pair("Z", "X"), std::make_pair( 1.0i, "Y")},
        {std::make_pair("Z", "Y"), std::make_pair(-1.0i, "X")},
        {std::make_pair("X", "X"), std::make_pair( 1.0,  "I")},
        {std::make_pair("Y", "Y"), std::make_pair( 1.0,  "I")},
        {std::make_pair("Z", "Z"), std::make_pair( 1.0,  "I")},
        {std::make_pair("I", "X"), std::make_pair( 1.0,  "X")},
        {std::make_pair("I", "Y"), std::make_pair( 1.0,  "Y")},
        {std::make_pair("I", "Z"), std::make_pair( 1.0,  "Z")}
    };

    std::sort(gates_.begin(), gates_.end(),
        [&](const QuantumGate& a, const QuantumGate& b) {
            return (a.target() < b.target());
        }
    );

    int n_gates = gates_.size();
    QuantumCircuit simplified_circ;
    std::complex<double> coeff = 1.0;
    std::string s;
    // counts contractions of gates for the same qubit
    int counter = 0;

    for (int i=0; i<n_gates-1; i++) {
        if(gates_[i].target() == gates_[i+1].target()) {
            if (counter == 0) {
                s = gates_[i].gate_id();
            }
            coeff *= m[std::make_pair(
                s,
                gates_[i+1].gate_id()
            )].first;

            s = m[std::make_pair(
                s,
                gates_[i+1].gate_id()
            )].second;
            counter++;

        } else {
            // ith gate is only gate to act on target qubit
            if (counter == 0) {
                simplified_circ.add_gate(
                    make_gate(gates_[i].gate_id(), gates_[i].target(), gates_[i].target())
                );
            } else {
                simplified_circ.add_gate(
                    make_gate(s, gates_[i].target(), gates_[i].target())
                );
                counter = 0;
            }
        }
    }

    std::cout << "i: " << "f" << " s: " << s << " c: " << counter << std::endl;

    // imples last elemet is different qbit that the rest
    if (counter == 0){
        simplified_circ.add_gate(
            make_gate(gates_[n_gates-1].gate_id(),
            gates_[n_gates-1].target(),
            gates_[n_gates-1].target())
        );
    // } else /*if (s != "I")*/{
    } else if (s != "I"){
        simplified_circ.add_gate(
            make_gate(s,
            gates_[n_gates-1].target(),
            gates_[n_gates-1].target())
        );
    }

    //copy simplified terms_
    // maybe copy, mc simp_circ may not be continuous in mem.
    gates_ = std::move(simplified_circ.gates());
    return coeff;
}

std::string QuantumCircuit::str() const {
    std::vector<std::string> s;
    for (const auto& gate : gates_) {
        s.push_back(gate.str());
    }
    return "[" + join(s, " ") + "]";
}

// std::vector<double> QuantumCircuit::get_parameters() {
//     // need a loop over only gates in state preparation circuit that
//     // have a parameter dependance (if gate_id == Rx, Ry, or Rz)
//     // TODO: make a indexing funciton using a map (Nick)
//     size_t param_idx = 0;
//     std::vector<double> params
//     for (auto& gate : gates_) {
//         std::string gate_id = gate.gate_id();
//         if (gate_id == "Rz") {
//
//             double param = gate.gate()[][];
//             gate = make_gate(gate_id, target_qubit, target_qubit, params[param_idx]);
//             param_idx++;
//         }
//     }
// }
