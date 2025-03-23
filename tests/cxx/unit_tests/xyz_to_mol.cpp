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
#include "catch2/catch_test_macros.hpp"
#include "module_manager/module_manager_class.hpp"
#include <chemcache/chemcache_mm.hpp>
#include <fstream>
#include <nux/nux.hpp>

TEST_CASE("XYZToMolecule") {

  SECTION("Is this working?") {
    pluginplay::ModuleManager mm_test;
    nux::load_modules(mm_test);

    REQUIRE(mm_test.size() > 0);
  }

  SECTION("Another One") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    nux::load_modules(mm);

    mm.change_submod("XYZ To Molecule", "Z", "Z from Symbol");
    mm.change_submod("XYZ To Molecule", "A", "Atom");

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

    auto mol = mm.at("XYZ To Molecule").run_as<simde::MoleculeFromString>(filename);
    
    std::cout << "TEST RAN!\n";
    REQUIRE(mol == test_mol);
  }
}
