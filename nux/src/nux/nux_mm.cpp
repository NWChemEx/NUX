#include "nux/nux_mm.hpp"

namespace nux {

inline void set_defaults(pluginplay::ModuleManager& mm) {
    // Default submodules between collections can be set here
}

DECLARE_PLUGIN(nux) {
    // Add subcollection load calls here

    // Assign default submodules
    set_defaults(mm);
}

} // namespace nux
