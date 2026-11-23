/*
 * Copyright 2026 NWChemEx-Project
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
#include <simde/derivative/derivative.hpp>

using namespace nux;

namespace {

using pt        = simde::EnergyNuclearGradientStdVectorD;
using energy_pt = simde::TotalEnergy;
using simde::type::chemical_system;
using point_set = chemist::PointSet<double>;

// E = sum of the cube of every nuclear coordinate, so dE/dx = 3x^2. The second
// order stencil is off by exactly h^2; higher orders are exact.
auto cubic_energy() {
    return pluginplay::make_lambda<energy_pt>([](const chemical_system& sys) {
        double e = 0.0;
        for(const auto& nucleus : sys.molecule().nuclei())
            for(std::size_t q = 0; q < 3; ++q)
                e += nucleus.coord(q) * nucleus.coord(q) * nucleus.coord(q);
        return simde::type::tensor(e);
    });
}

// Coordinates are distinct so the test is sensitive to the gradient's layout
const std::vector<double> coords{0.1, -0.2, 0.3, 0.4, 0.5, 1.6};

auto make_system() {
    simde::type::nuclei nuclei{
      test_nux::h_nucleus(coords[0], coords[1], coords[2]),
      test_nux::h_nucleus(coords[3], coords[4], coords[5])};
    return chemical_system(simde::type::molecule(0, 1, nuclei));
}

auto make_points() {
    point_set points;
    points.push_back(chemist::Point<double>(coords[0], coords[1], coords[2]));
    points.push_back(chemist::Point<double>(coords[3], coords[4], coords[5]));
    return points;
}

} // namespace

TEST_CASE("CentralDifferenceGradient") {
    pluginplay::ModuleManager mm;
    mm.load_plugin<NUX>();
    auto base = mm.at("Central Difference Gradient").unlocked_copy();
    base.change_submod("Energy", cubic_energy());

    const auto sys    = make_system();
    const auto points = make_points();

    auto run = [&](std::size_t order, double h) {
        auto mod = base.unlocked_copy();
        mod.change_input("Accuracy Order", order);
        mod.change_input("Step Size", h);
        return mod.run_as<pt>(sys, points);
    };

    SECTION("Defaults") {
        auto grad      = base.unlocked_copy().run_as<pt>(sys, points);
        const double h = 0.005;
        REQUIRE(grad.size() == coords.size());
        for(std::size_t i = 0; i < coords.size(); ++i)
            REQUIRE(grad[i] ==
                    Catch::Approx(3.0 * coords[i] * coords[i] + h * h));
    }

    SECTION("Order 2 error is h^2") {
        for(double h : {0.01, 0.1}) {
            auto grad = run(2, h);
            REQUIRE(grad.size() == coords.size());
            for(std::size_t i = 0; i < coords.size(); ++i)
                REQUIRE(grad[i] ==
                        Catch::Approx(3.0 * coords[i] * coords[i] + h * h));
        }
    }

    SECTION("Orders 4, 6, and 8 are exact for a cubic") {
        for(std::size_t order : {4, 6, 8}) {
            auto grad = run(order, 0.1);
            REQUIRE(grad.size() == coords.size());
            for(std::size_t i = 0; i < coords.size(); ++i)
                REQUIRE(
                  grad[i] ==
                  Catch::Approx(3.0 * coords[i] * coords[i]).margin(1.0E-10));
        }
    }

    SECTION("Points must match the chemical system") {
        auto mod = base.unlocked_copy();

        point_set too_few;
        too_few.push_back(points.at(0).as_point());
        REQUIRE_THROWS_AS(mod.run_as<pt>(sys, too_few), std::runtime_error);

        auto moved = points;
        moved.at(1).z() += 1.0E-12;
        REQUIRE_THROWS_AS(mod.run_as<pt>(sys, moved), std::runtime_error);
    }

    SECTION("Invalid inputs") {
        auto mod = base.unlocked_copy();
        REQUIRE_THROWS(mod.change_input("Accuracy Order", std::size_t{3}));
        REQUIRE_THROWS(mod.change_input("Accuracy Order", std::size_t{10}));
        REQUIRE_THROWS(mod.change_input("Step Size", 0.0));
        REQUIRE_THROWS(mod.change_input("Step Size", -0.005));
    }
}
