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

#include "nux_modules.hpp"
#include <map>
#include <simde/derivative/derivative.hpp>
#include <stdexcept>
#include <vector>

namespace nux {
namespace {

using pt        = simde::EnergyNuclearGradientStdVectorD;
using energy_pt = simde::TotalEnergy;

/// Central-difference coefficients c_k, k = 1, ..., order / 2, keyed by order
const std::map<std::size_t, std::vector<double>> coefficients{
  {2, {1.0 / 2.0}},
  {4, {2.0 / 3.0, -1.0 / 12.0}},
  {6, {3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0}},
  {8, {4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0}}};

/// Unwraps the scalar energy stored in the rank 0 tensor @p e
double to_double(const simde::type::tensor& e) {
    using tensorwrapper::buffer::make_contiguous;
    return wtf::fp::float_cast<double>(
      make_contiguous(e.buffer()).get_elem({}));
}

const auto desc = R"(
Central Difference Nuclear Gradient
-----------------------------------

This module computes the gradient of the total energy, with respect to the
nuclear coordinates, by central finite differences. The energy is computed by
the "Energy" submodule, which is called on copies of the input chemical system
where a single nuclear coordinate has been displaced.

For a coordinate :math:`x_i` the gradient component is:

.. math::

   \frac{\partial E}{\partial x_i} \approx
   \frac{1}{h}\sum_{k=1}^{p/2} c_k
   \left[E(x_i + kh) - E(x_i - kh)\right]

where :math:`h` is the "Step Size" input (in bohr) and :math:`p` is the
"Accuracy Order" input (the error is :math:`\mathcal{O}(h^p)`). Supported
orders and their coefficients are:

- 2: 1/2
- 4: 2/3, -1/12
- 6: 3/4, -3/20, 1/60
- 8: 4/5, -1/5, 4/105, -1/280

Computing the gradient of an :math:`N` atom system requires :math:`3Np`
energy evaluations.

The gradient is returned in atom-major order, i.e., x, y, z of atom 0 followed
by x, y, z of atom 1, etc. Derivatives are taken at the points in the PointSet
input, which must exactly match the nuclear coordinates of the chemical system.
)";

} // namespace

MODULE_CTOR(CentralDifferenceGradient) {
    satisfies_property_type<pt>();
    description(desc);

    add_submodule<energy_pt>("Energy").set_description(
      "Computes the energy of each displaced chemical system");

    add_input<std::size_t>("Accuracy Order")
      .set_default(std::size_t{2})
      .set_description("Order of the error in the step size. Determines the "
                       "number of points in the stencil (order + 1).")
      .add_check(
        std::function<bool(const std::size_t&)>(
          [](const std::size_t& order) { return coefficients.count(order); }),
        "Accuracy order is 2, 4, 6, or 8");

    add_input<double>("Step Size")
      .set_default(0.005)
      .set_description("Displacement of each coordinate in bohr")
      .add_check(pluginplay::bounds_checking::GreaterThan<double>(0.0));
}

MODULE_RUN(CentralDifferenceGradient) {
    const auto& [sys, points] = pt::unwrap_inputs(inputs);
    const auto order = inputs.at("Accuracy Order").value<std::size_t>();
    const auto h     = inputs.at("Step Size").value<double>();

    // Step 0: Derivatives are only meaningful at the system's geometry
    const auto nuclei   = sys.molecule().nuclei();
    const auto n_nuclei = nuclei.size();
    const std::string error_msg =
      "Derivative must be computed at the chemical system's geometry";
    if(points.size() != n_nuclei) throw std::runtime_error(error_msg);
    for(std::size_t atom = 0; atom < n_nuclei; ++atom)
        for(std::size_t q = 0; q < 3; ++q)
            if(points.at(atom).coord(q) != nuclei.at(atom).coord(q))
                throw std::runtime_error(error_msg);

    // Step 1: Energy of sys with nuclear coordinate q of atom displaced by dx
    auto& energy_mod      = submods.at("Energy");
    auto displaced_energy = [&](std::size_t atom, std::size_t q, double dx) {
        auto displaced = sys;
        displaced.molecule().at(atom).coord(q) += dx;
        return to_double(energy_mod.run_as<energy_pt>(displaced));
    };

    // Step 2: Central differences for each coordinate (atom-major order)
    const auto& c = coefficients.at(order);
    std::vector<double> grad(3 * n_nuclei, 0.0);
    for(std::size_t atom = 0; atom < n_nuclei; ++atom) {
        for(std::size_t q = 0; q < 3; ++q) {
            double g = 0.0;
            for(std::size_t k = 1; k <= c.size(); ++k) {
                const auto dx = static_cast<double>(k) * h;
                g += c[k - 1] * (displaced_energy(atom, q, dx) -
                                 displaced_energy(atom, q, -dx));
            }
            grad[3 * atom + q] = g / h;
        }
    }

    auto rv = results();
    return pt::wrap_results(rv, grad);
}

} // namespace nux
