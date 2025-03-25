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
#pragma once
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <simde/simde.hpp>

namespace test_nux {

inline auto h_nucleus(double x, double y, double z) {
    return simde::type::nucleus("H", 1ul, 1836.15, x, y, z);
}

template<typename ResultType>
ResultType make_h2() {
    using simde::type::chemical_system;
    using simde::type::molecule;
    using simde::type::nuclei;
    if constexpr(std::is_same_v<ResultType, nuclei>) {
        auto h0 = h_nucleus(0.0, 0.0, 0.0);
        auto h1 = h_nucleus(0.0, 0.0, 1.3984);
        return nuclei{h0, h1};
    } else if constexpr(std::is_same_v<ResultType, molecule>) {
        return molecule(0, 1, make_h2<nuclei>());
    } else if constexpr(std::is_same_v<ResultType, chemical_system>) {
        return chemical_system(make_h2<molecule>());
    } else {
        // We know this assert fails if we're in the else statement
        // Getting here means you provided a bad type.
        static_assert(std::is_same_v<ResultType, nuclei>);
    }
}
} // namespace test_nux