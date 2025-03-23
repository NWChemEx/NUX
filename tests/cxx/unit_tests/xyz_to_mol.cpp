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
#include "chemcache_mm.hpp"
#include "chemical_system/molecule_from_string.hpp"
#include <nux/nux.hpp>

using namespace nux;

TEST_CASE("XYZToMolecule") {
  SECTION("TEST") {
    pluginplay::ModuleManager mm;
    chemcache::load_modules(mm);
    nux::load_modules(mm);
    std::cout << "TEST IS RUNNING\n";
    auto mol = mm.run_as<simde::MoleculeFromString>("XYZ To Molecule", "h2.xyz");
    std::cout << mol.nuclei();
    REQUIRE(1 == 1);
    }
  }
