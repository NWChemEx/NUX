
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

TEST_CASE("XYZFileToMolecule") {
    pluginplay::ModuleManager mm;
    nux::load_modules(mm);

    auto make_atoms = [=](auto&& Z) {
        REQUIRE(Z == 1);
        double h_mass = 1837.4260218693814;
        return simde::type::atom{"H", 1, h_mass, 0.0, 0.0, 0.0};
    };

    auto atom_mod = pluginplay::make_lambda<simde::AtomFromZ>(make_atoms);

    auto z_mod = pluginplay::make_lambda<simde::ZFromSymbol>(
      [=](auto&& symbol) { return simde::type::atomic_number{1}; });

    auto& xyz_file_mod = mm.at("XYZ File To Molecule");
    auto& xyz_mod      = mm.at("XYZ to Molecule");
    xyz_mod.change_submod("Z from symbol", z_mod);
    xyz_mod.change_submod("Atom from z", atom_mod);

    SECTION("XYZ File Non-existent") {
        std::string file = "file.xyz";
        REQUIRE_THROWS_AS(xyz_file_mod.run_as<simde::MoleculeFromString>(file),
                          std::runtime_error);
    }

    SECTION("Full XYZ To Molecule run: File") {
        auto atom0{make_atoms(1)};
        auto atom1{make_atoms(1)};
        atom1.z() = 1;

        simde::type::molecule test_mol{atom0, atom1};

        std::ofstream xyz_file;
        xyz_file.open("h2.xyz");

        xyz_file << "2\n";
        xyz_file << "This is a comment!\n";
        xyz_file << "H 0 0 0\n";
        xyz_file << "H 0 0 1\n";
        xyz_file.close();

        std::string filename = "h2.xyz";

        auto& xyz_file_mod = mm.at("XYZ File to Molecule");
        auto mol = xyz_file_mod.run_as<simde::MoleculeFromString>(filename);

        remove("h2.xyz");
        std::ifstream del_file("h2.xyz");

        REQUIRE(mol == test_mol);
        REQUIRE(del_file.good() == false);
    }
}
