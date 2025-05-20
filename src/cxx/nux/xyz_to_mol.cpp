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

namespace nux {

MODULE_CTOR(XYZToMolecule) {
    satisfies_property_type<simde::MoleculeFromString>();
    add_submodule<simde::ZFromSymbol>("Z from symbol");
    add_submodule<simde::AtomFromZ>("Atom from z");
}

MODULE_RUN(XYZToMolecule) {
    const auto& [xyz_data] = simde::MoleculeFromString::unwrap_inputs(inputs);
    auto& z_from_sym       = submods.at("Z from symbol");
    auto& atom_from_z      = submods.at("Atom from z");

    std::string line;
    chemist::Molecule mol;

    std::stringstream buffer;

    buffer << xyz_data;
    if(buffer.str().empty()) {
        throw std::runtime_error("File or XYZ data not valid, empty string");
    }

    int i = 0;

    while(std::getline(buffer, line)) {
        i++;
        int num_atoms;
        if(i == 1) {
            try {
                num_atoms = std::stoi(line);
            } catch(const std::invalid_argument& e) {
                throw std::runtime_error("Formatting of XYZ data incorrect");
            }
            continue;
        }
        if(i == 2) { continue; }
        std::istringstream iss(line);
        std::string atom_string;
        double x, y, z;
        if(!(iss >> atom_string >> x >> y >> z)) {
            throw std::runtime_error("This is not a good coordinate");
        }

        auto Z    = z_from_sym.run_as<simde::ZFromSymbol>(atom_string);
        auto atom = atom_from_z.run_as<simde::AtomFromZ>(Z);
        atom.x()  = x;
        atom.y()  = y;
        atom.z()  = z;

        mol.push_back(atom);

        if(!(buffer.peek() == EOF)) {
            continue;
        } else {
            if(mol.size() == num_atoms) {
                continue;
            } else {
                throw std::runtime_error("Incorrect format for XYZ file: "
                                         "Number of atoms and number of atom "
                                         "coordinates do not match");
            }
        }
    }
    auto rv = results();
    return simde::MoleculeFromString::wrap_results(rv, mol);
}

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
    try {
        if(file.is_open()) {
            buffer << file.rdbuf();
            file.close();
        } else {
            throw std::runtime_error("");
        }
    } catch(const std::runtime_error& e) {
        std::cerr << "File not found: " << xyz_filename << std::endl;
    }

    chemist::Molecule mol =
      xyz_parser.run_as<simde::MoleculeFromString>(buffer.str());

    auto rv = results();
    return simde::MoleculeFromString::wrap_results(rv, mol);
}
} // namespace nux
