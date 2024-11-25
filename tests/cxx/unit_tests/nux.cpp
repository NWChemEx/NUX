#include "../test_nux.hpp"
#include <nux/nux.hpp>
TEST_CASE("load_plugin") {
  pluginplay::ModuleManager mm;
  nux::load_modules(mm);
}