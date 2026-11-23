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
#include <cmath>
#include <nux/nux.hpp>
#include <simde/derivative/derivative.hpp>

using namespace nux;

namespace {

using pt        = simde::EnergyNuclearGradientTensor;
using energy_pt = simde::TotalEnergy;
using simde::type::chemical_system;
using simde::type::tensor;
using point_set = chemist::PointSet<double>;

// E = sum of the cube of every nuclear coordinate, so dE/dx = 3x^2. The second
// order stencil is off by exactly h^2; higher orders are exact.
double cubic(const chemical_system& sys) {
    double e = 0.0;
    for(const auto& nucleus : sys.molecule().nuclei())
        for(std::size_t q = 0; q < 3; ++q)
            e += nucleus.coord(q) * nucleus.coord(q) * nucleus.coord(q);
    return e;
}

auto cubic_energy() {
    return pluginplay::make_lambda<energy_pt>(
      [](const chemical_system& sys) { return tensor(cubic(sys)); });
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

// The 2 by 3 gradient's elements as doubles, in row-major order
std::vector<double> elements(const tensor& grad) {
    using tensorwrapper::buffer::make_contiguous;
    const auto& buffer = make_contiguous(grad.buffer());
    const auto shape   = buffer.shape();
    REQUIRE(shape.rank() == 2);
    REQUIRE(shape.extent(0) == 2);
    REQUIRE(shape.extent(1) == 3);
    std::vector<double> rv;
    for(std::size_t atom = 0; atom < 2; ++atom)
        for(std::size_t q = 0; q < 3; ++q)
            rv.push_back(
              wtf::fp::float_cast<double>(buffer.get_elem({atom, q})));
    return rv;
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
        auto grad      = elements(base.unlocked_copy().run_as<pt>(sys, points));
        const double h = 0.005;
        for(std::size_t i = 0; i < coords.size(); ++i)
            REQUIRE(grad[i] ==
                    Catch::Approx(3.0 * coords[i] * coords[i] + h * h));
    }

    SECTION("Order 2 error is h^2") {
        for(double h : {0.01, 0.1}) {
            auto grad = elements(run(2, h));
            for(std::size_t i = 0; i < coords.size(); ++i)
                REQUIRE(grad[i] ==
                        Catch::Approx(3.0 * coords[i] * coords[i] + h * h));
        }
    }

    SECTION("Orders 4, 6, and 8 are exact for a cubic") {
        for(std::size_t order : {4, 6, 8}) {
            auto grad = elements(run(order, 0.1));
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

    SECTION("System without nuclei") {
        auto mod = base.unlocked_copy();
        REQUIRE_THROWS_AS(mod.run_as<pt>(chemical_system{}, point_set{}),
                          std::runtime_error);
    }

    SECTION("Invalid inputs") {
        auto mod = base.unlocked_copy();
        REQUIRE_THROWS(mod.change_input("Accuracy Order", std::size_t{3}));
        REQUIRE_THROWS(mod.change_input("Accuracy Order", std::size_t{10}));
        REQUIRE_THROWS(mod.change_input("Step Size", 0.0));
        REQUIRE_THROWS(mod.change_input("Step Size", -0.005));
    }

#ifdef ENABLE_SIGMA
    SECTION("Uncertain energies give an uncertain gradient") {
        using tensorwrapper::types::udouble;
        // Each energy is an independent random variable with this std. dev.
        const double sigma = 1.0E-6;
        base.change_submod(
          "Energy",
          pluginplay::make_lambda<energy_pt>([=](const chemical_system& s) {
              tensorwrapper::shape::Smooth scalar{};
              std::vector<udouble> e{udouble(cubic(s), sigma)};
              tensorwrapper::buffer::Contiguous buffer(std::move(e), scalar);
              return tensor(scalar, std::move(buffer));
          }));

        const double h = 0.01;
        for(std::size_t order : {2, 4}) {
            auto grad = run(order, h);
            using tensorwrapper::buffer::get_raw_data;
            auto values = get_raw_data<udouble>(grad.buffer());
            REQUIRE(values.size() == coords.size());

            // Mean is the double result; std. dev. is sqrt(2 sum c_k^2) s / h
            const double c2    = order == 2 ? 0.25 : 4.0 / 9.0 + 1.0 / 144.0;
            const double sd    = std::sqrt(2.0 * c2) * sigma / h;
            const double error = order == 2 ? h * h : 0.0;
            for(std::size_t i = 0; i < coords.size(); ++i) {
                REQUIRE(values[i].mean() ==
                        Catch::Approx(3.0 * coords[i] * coords[i] + error)
                          .margin(1.0E-10));
                REQUIRE(values[i].sd() == Catch::Approx(sd));
            }
        }
    }
#endif
}
