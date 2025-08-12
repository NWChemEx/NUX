# Copyright 2024 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import chemist


def h_nucleus(x, y, z):
    return chemist.Nucleus("H", 1, 1836.15, x, y, z)


def h2_nuclei():
    h0 = h_nucleus(0.0, 0.0, 0.0)
    h1 = h_nucleus(0.0, 0.0, 1.3984)

    return chemist.Nuclei(h0, h1)
