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

#include "../test_nux.hpp"
#include <nux/nux.hpp>

TEST_CASE("XYZToMolecule") {
    pluginplay::ModuleManager mm;
    nux::load_modules(mm);

    auto atom_mod = pluginplay::make_lambda<simde::AtomFromZ>([=](auto&& Z) {
        double h_mass = 1837.4260218693814;
        return simde::type::atom{"H", 1, h_mass, 0.0, 0.0, 0.0};
    });
    

    auto z_mod = pluginplay::make_lambda<simde::ZFromSymbol>([=](auto&& symbol) {
        return simde::type::atomic_number{1};
    });

    auto& xyz_mod = mm.at("XYZ To Molecule");
    xyz_mod.change_submod("Z from symbol", z_mod);
    xyz_mod.change_submod("Atom from z", atom_mod);

  SECTION("No XYZ file found") {
    std::string file = "file.xyz";
    REQUIRE_THROWS_AS(xyz_mod.run_as<simde::MoleculeFromString>(file), std::runtime_error);
  }

  SECTION("Full XYZ To Molecule run") {
    chemist::Molecule test_mol;
    chemist::Atom atom("H", 1, 1837.42602186938, 0, 0, 0);
    chemist::Atom atom2("H", 1, 1837.42602186938, 0, 0, 1);

    test_mol.push_back(atom);
    test_mol.push_back(atom2);

    std::ofstream xyz_file;
    xyz_file.open("h2.xyz");

    xyz_file << "2\n";
    xyz_file << "This is a comment!\n";
    xyz_file << "H 0 0 0\n";
    xyz_file << "H 0 0 1\n";
    xyz_file.close();
    
    std::string filename = "h2.xyz";

    auto mol = xyz_mod.run_as<simde::MoleculeFromString>(filename);

    remove("h2.xyz");
    std::ifstream del_file("h2.xyz");

    REQUIRE(mol == test_mol);
    REQUIRE(del_file.good() == false);
    
  }
}
