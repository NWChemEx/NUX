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

namespace {

// Stands in for ChemCache, which NUX gets atomic data from
class FakeChemCache : public pluginplay::Plugin {
public:
    FakeChemCache() {
        set_name("chemcache");
        auto z_from_sym = pluginplay::make_lambda<simde::ZFromSymbol>(
          [](auto&&) { return simde::type::atomic_number{1}; });
        auto atom_from_z =
          pluginplay::make_lambda<simde::AtomFromZ>([](auto&&) {
              return simde::type::atom{"H", 1ul, 1836.15, 0.0, 0.0, 0.0};
          });
        modules().add_module("Z from Symbol", z_from_sym);
        modules().add_module("Atom", atom_from_z);
    }
};

} // namespace

TEST_CASE("NUX plugin") {
    nux::NUX plugin;
    const auto& modules = plugin.get_modules();

    SECTION("Identity") {
        REQUIRE(plugin.name() == "nux");
        using provider_list = nux::NUX::provider_list_type;
        REQUIRE(plugin.preferred_providers() == provider_list{"chemcache"});
    }

    SECTION("Modules") {
        REQUIRE(modules.size() > 1);
        REQUIRE(modules.count("Born-Oppenheimer Approximation"));
        REQUIRE(modules.count("XYZ To Molecule"));
        REQUIRE(modules.count("XYZ File to Molecule"));
    }

    SECTION("Defaults set within the plugin") {
        REQUIRE(modules.has_submod("XYZ File to Molecule", "XYZ to molecule"));
    }

    pluginplay::ModuleManager mm;
    mm.load_plugin<nux::NUX>();
    const auto& xyz_submods = mm.at("XYZ To Molecule").submods();

    SECTION("Can be loaded on its own") {
        REQUIRE(mm.has_plugin("nux"));
        REQUIRE(mm.size() == modules.size());
        mm.wire_defaults();
        REQUIRE_FALSE(xyz_submods.at("Z from symbol").has_module());
        REQUIRE_FALSE(xyz_submods.at("Atom from z").has_module());
    }

    SECTION("Gets atomic data from ChemCache") {
        mm.load_plugin<FakeChemCache>();
        mm.wire_defaults();
        REQUIRE(xyz_submods.at("Z from symbol").value() ==
                mm.at("Z from Symbol"));
        REQUIRE(xyz_submods.at("Atom from z").value() == mm.at("Atom"));
    }
}
