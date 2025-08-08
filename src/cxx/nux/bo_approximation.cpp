/*
 * Copyright 2024 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "nux_modules.hpp"

namespace nux {

using simde::type::chemical_system;
using simde::type::hamiltonian;
using simde::type::many_electrons;
using pt = simde::Convert<hamiltonian, chemical_system>;

const auto desc = R"(
Born-Oppenheimer Hamiltonian Driver
-----------------------------------

This module defines the algorithm used throughout NWChemEx to convert a
ChemicalSystem object into a Hamiltonian object using the Born-Oppenheimer
approximation.

Specifically this driver will convert:
- ManyElectrons to a T_e and V_ee terms.
- Nuclei to a V_nn
- ManyElectrons + Nuclei to a V_en term.

The terms will be added following the "usual" ordering:

- Electronic Hamiltonian
   - Electronic kinetic energy
   - Electron-nucleus attraction
   - Electron-electron repulsion
- Nuclear-nuclear repulsion

Terms that are zero will not appear in the Hamiltonian.
)";

MODULE_CTOR(BOApproximation) {
    satisfies_property_type<pt>();
    description(desc);
}

MODULE_RUN(BOApproximation) {
    const auto& [sys] = pt::unwrap_inputs(inputs);

    // Step 0: Decompose system into pieces
    const auto electrons = sys.molecule().electrons();
    const auto nuclei    = sys.molecule().nuclei().as_nuclei();
    bool has_electrons   = electrons.size();
    bool has_nuclei      = nuclei.size();

    using simde::type::T_e_type;
    using simde::type::V_ee_type;
    using simde::type::V_en_type;
    using simde::type::V_nn_type;

    // Step 1: Create non-zero terms and add them to the Hamiltonian
    // N.b. the logic is a bit convoluted to conform to "usual" ordering (see
    // module description)

    hamiltonian H;
    if(has_electrons) {
        H.emplace_back(1.0, std::make_unique<T_e_type>(electrons));
        if(has_nuclei)
            H.emplace_back(1.0, std::make_unique<V_en_type>(electrons, nuclei));
        if(electrons.size() > 1)
            H.emplace_back(1.0,
                           std::make_unique<V_ee_type>(electrons, electrons));
    }
    if(nuclei.size() > 1) {
        H.emplace_back(1.0, std::make_unique<V_nn_type>(nuclei, nuclei));
    }

    auto rv = results();
    return pt::wrap_results(rv, H);
}

} // namespace nux
