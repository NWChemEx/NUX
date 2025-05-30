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

#include "nux_modules.hpp"
#include <nux/nux.hpp>

namespace nux {

void load_modules(pluginplay::ModuleManager& mm) {
    mm.add_module<BOApproximation>("Born-Oppenheimer Approximation");
    mm.add_module<XYZToMolecule>("XYZ To Molecule");
    mm.add_module<XYZFileToMolecule>("XYZ File to Molecule");

    set_defaults(mm);
}

} // namespace nux
