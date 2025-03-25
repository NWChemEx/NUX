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
#include <string>
#include <vector>

namespace nux {

  MODULE_CTOR(XYZToMolecule) {
    satisfies_property_type<simde::MoleculeFromString>();
    add_submodule<simde::ZFromSymbol>("Z from symbol");
    add_submodule<simde::AtomFromZ>("Atom from z");
  }

  MODULE_RUN(XYZToMolecule) {
    const auto& [filename] = simde::MoleculeFromString::unwrap_inputs(inputs);

    auto& z_from_sym = submods.at("Z from symbol");
    auto& atom_from_z = submods.at("Atom from z");

    chemist::Molecule mol;

    std::ifstream xyz_file(filename);

    if (!xyz_file) {
      throw std::runtime_error("File not found: " + filename);
    }
    
    std::string line;

    std::getline(xyz_file, line);
    std::getline(xyz_file, line);

    while (std::getline(xyz_file, line)) {
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string token;
      while (std::getline(iss, token, ' ')) {
        if (token.size()) { tokens.emplace_back(token); }
      }
      auto Z = z_from_sym.run_as<simde::ZFromSymbol>(tokens.at(0));
      auto x = std::stod(tokens.at(1));
      auto y = std::stod(tokens.at(2));
      auto z = std::stod(tokens.at(3));

      auto atm = atom_from_z.run_as<simde::AtomFromZ>(Z);
      atm.x()  = x;
      atm.y()  = y;
      atm.z()  = z;
      mol.push_back(atm);
    }

    xyz_file.close();
    
    auto rv = results();
    return simde::MoleculeFromString::wrap_results(rv, mol);
  }

} // namespace nux
