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
#include <pluginplay/pluginplay.hpp>

namespace nux {
/** @brief The NUX plugin.
 *
 *  Load it into a ModuleManager with `mm.load_plugin<nux::NUX>()`
 *  (`mm.load_plugin(nux.NUX())` from Python). The plugin's name is "nux". If
 *  the ChemCache plugin is also loaded, ModuleManager::wire_defaults uses its
 *  modules for the atomic data needed by "XYZ To Molecule".
 */
DECLARE_PLUGIN(NUX);

} // namespace nux
