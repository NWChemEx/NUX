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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "chemcache/chemcache.hpp"
#include "chemical_system/Z_from_symbol.hpp"
#include "chemical_system/atom.hpp"
#include "chemical_system/molecule_from_string.hpp"
#include "module/macros.hpp"
#include "module_manager/module_manager_class.hpp"
#include "molecule/molecule_class.hpp"
#include "nux_modules.hpp"

namespace nux {

MODULE_CTOR(XYZToMolecule) {
  satisfies_property_type<simde::MoleculeFromString>();
  add_submodule<simde::ZFromSymbol>("Z");
  add_submodule<simde::AtomFromZ>("Atom");
}

MODULE_RUN(XYZToMolecule){
  std::cout << "GETTING FILENAME\n";
  const auto& [filename] = simde::MoleculeFromString::unwrap_inputs(inputs);

  std::cout << "OPENING FILE\n";
  std::ifstream xyz_file(filename);

  chemist::Molecule mol;
  
  std::string line;
  
  auto& z_from_symbol = submods.at("Z");
  auto& z_to_atom = submods.at("Atom");


  std::cout << "ASSESSING IF FILE IS OPEN\n";
  if (!xyz_file.is_open()){
    std::cerr << "Error opening file: " << filename;
  } else {
    std::cout << "GOING TO LINE 3\n";
    xyz_file.seekg(2);


    std::cout << "SETTING ATOMS\n";
    while (std::getline(xyz_file, line)) {
      std::vector<std::string> tokens;
      std::string token;
      std::istringstream iss(line);

      while (std::getline(iss, token, ' ')) {
        if (token.size()) { tokens.emplace_back(token); }
      }

      auto Z = z_from_symbol.run_as<simde::ZFromSymbol>(tokens.at(0));
      auto x = std::stod(tokens.at(1));
      auto y = std::stod(tokens.at(2));
      auto z = std::stod(tokens.at(3));

      auto atom = z_to_atom.run_as<simde::AtomFromZ>(Z);
      atom.x() = x;
      atom.y() = y;
      atom.z() = z;

      mol.push_back(atom);
    }
  std::cout << "PROCESS COMPLETED\n";
  }
  auto rv = results();
  
  return simde::MoleculeFromString::wrap_results(rv, mol);
}

} // namespace nux
