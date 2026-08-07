
/*
 * Copyright 2025 NWChemEx-Project
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
#include <fstream>
#include <sstream>

namespace nux {

MODULE_CTOR(XYZFileToMolecule) {
    satisfies_property_type<simde::MoleculeFromString>();
    add_submodule<simde::MoleculeFromString>("XYZ to molecule");
}

MODULE_RUN(XYZFileToMolecule) {
    const auto& [xyz_filename] =
      simde::MoleculeFromString::unwrap_inputs(inputs);

    auto& xyz_parser = submods.at("XYZ to molecule");

    std::ifstream file(xyz_filename);
    std::stringstream buffer;

    if(!file.is_open()) {
        throw std::runtime_error("File not found: " + xyz_filename);
    }
    buffer << file.rdbuf();
    file.close();

    chemist::Molecule mol =
      xyz_parser.run_as<simde::MoleculeFromString>(buffer.str());

    auto rv = results();
    return simde::MoleculeFromString::wrap_results(rv, mol);
}
} // namespace nux
