#include "../test_nux.hpp"
#include <nux/nux.hpp>

using namespace nux;

TEST_CASE("BOApproximation") {
  pluginplay::ModuleManager mm;
  load_modules(mm);

  using simde::type::chemical_system;
  using simde::type::hamiltonian;
  using simde::type::T_e_type;
  using simde::type::V_ee_type;
  using simde::type::V_en_type;
  using simde::type::V_nn_type;

  using pt = simde::Convert<hamiltonian, chemical_system>;

  auto &mod = mm.at("Born-Oppenheimer Approximation");

  SECTION("Empty chemical system") {
    chemical_system sys;
    auto H = mod.run_as<pt>(sys);
    REQUIRE(H == hamiltonian{});
  }

  SECTION("H2") {
    auto h2 = test_nux::make_h2<chemical_system>();
    auto electrons = h2.molecule().electrons();
    auto nuclei = h2.molecule().nuclei().as_nuclei();

    auto H = mod.run_as<pt>(h2);
    hamiltonian H_corr;
    H_corr.emplace_back(1.0, std::make_unique<T_e_type>(electrons));
    H_corr.emplace_back(1.0, std::make_unique<V_en_type>(electrons, nuclei));
    H_corr.emplace_back(1.0, std::make_unique<V_ee_type>(electrons, electrons));
    H_corr.emplace_back(1.0, std::make_unique<V_nn_type>(nuclei, nuclei));
    REQUIRE(H == H_corr);
  }

  SECTION("H atom") {
    simde::type::nuclei nuclei{test_nux::h_nucleus(0.0, 0.0, 0.0)};
    simde::type::chemical_system h(simde::type::molecule(0, 2, nuclei));

    auto electrons = h.molecule().electrons();

    auto H = mod.run_as<pt>(h);
    hamiltonian H_corr;
    H_corr.emplace_back(1.0, std::make_unique<T_e_type>(electrons));
    H_corr.emplace_back(1.0, std::make_unique<V_en_type>(electrons, nuclei));
    REQUIRE(H == H_corr);
  }

  SECTION("H anion") {
    simde::type::nuclei nuclei{test_nux::h_nucleus(0.0, 0.0, 0.0)};
    simde::type::chemical_system h(simde::type::molecule(-1, 1, nuclei));

    auto electrons = h.molecule().electrons();

    auto H = mod.run_as<pt>(h);
    hamiltonian H_corr;
    H_corr.emplace_back(1.0, std::make_unique<T_e_type>(electrons));
    H_corr.emplace_back(1.0, std::make_unique<V_en_type>(electrons, nuclei));
    H_corr.emplace_back(1.0, std::make_unique<V_ee_type>(electrons, electrons));
    REQUIRE(H == H_corr);
  }

  SECTION("H2 cation") {
    auto nuclei = test_nux::make_h2<simde::type::nuclei>();
    simde::type::chemical_system h(simde::type::molecule(1, 2, nuclei));

    auto electrons = h.molecule().electrons();

    auto H = mod.run_as<pt>(h);
    hamiltonian H_corr;
    H_corr.emplace_back(1.0, std::make_unique<T_e_type>(electrons));
    H_corr.emplace_back(1.0, std::make_unique<V_en_type>(electrons, nuclei));
    H_corr.emplace_back(1.0, std::make_unique<V_nn_type>(nuclei, nuclei));
    REQUIRE(H == H_corr);
  }
}