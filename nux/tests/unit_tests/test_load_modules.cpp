#include <catch2/catch.hpp>
#include <nux/nux.hpp>

TEST_CASE("load_modules") {
    pluginplay::ModuleManager mm;
    nux::load_modules(mm);
}
