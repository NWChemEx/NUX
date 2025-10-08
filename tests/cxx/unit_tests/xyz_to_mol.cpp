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

    auto make_atoms = [=](auto&& Z) {
        REQUIRE(Z == 1);
        double h_mass = 1837.4260218693814;
        return simde::type::atom{"H", 1, h_mass, 0.0, 0.0, 0.0};
    };

    auto atom_mod = pluginplay::make_lambda<simde::AtomFromZ>(make_atoms);

    auto z_mod = pluginplay::make_lambda<simde::ZFromSymbol>(
      [=](auto&& symbol) { return simde::type::atomic_number{1}; });

    auto& xyz_mod = mm.at("XYZ To Molecule");
    xyz_mod.change_submod("Z from symbol", z_mod);
    xyz_mod.change_submod("Atom from z", atom_mod);

    SECTION("XYZ Atom Count Incorrect") {
        std::string data = "2\nComment\nH 0 0 0";
        REQUIRE_THROWS_AS(xyz_mod.run_as<simde::MoleculeFromString>(data),
                          std::runtime_error);
    }

    SECTION("XYZ Atom Coordinate Data Incorrect") {
        std::string data = "2\nComment\nH 0 0 K\nH 0 0 1";
        REQUIRE_THROWS_AS(xyz_mod.run_as<simde::MoleculeFromString>(data),
                          std::runtime_error);
    }

    SECTION("Full XYZ To Molecule run: Data") {
        auto atom0{make_atoms(1)};
        auto atom1{make_atoms(1)};
        atom1.z() = 1;

        simde::type::molecule test_mol{atom0, atom1};

        std::stringstream xyz_data;
        xyz_data << "2\n";
        xyz_data << "This is a comment!\n";
        xyz_data << "H 0 0 0\n";
        xyz_data << "H 0 0 1\n";

        auto mol = xyz_mod.run_as<simde::MoleculeFromString>(xyz_data.str());

        REQUIRE(mol == test_mol);
    }

    SECTION("XYZ to Molecule Unit Scaling: Angstroms to Bohr") {
        auto atom0{make_atoms(1)};
        auto atom1{make_atoms(1)};
        atom1.z() = 1.8897259886;

        simde::type::molecule test_mol{atom0, atom1};

        std::stringstream xyz_data;
        xyz_data << "2\n";
        xyz_data << "This is a comment!\n";
        xyz_data << "H 0 0 0\n";
        xyz_data << "H 0 0 1\n";

        xyz_mod.change_input("Unit scaling factor", 1.8897259886);

        auto mol = xyz_mod.run_as<simde::MoleculeFromString>(xyz_data.str());

        REQUIRE(mol == test_mol);
    }
}
