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

from pluginplay import ModuleManager
import nux
import simde
import chemist
import unittest


class TestBOApproximation(unittest.TestCase):

    def test_empty_chemical_system(self):
        sys = chemist.ChemicalSystem()
        # Hamiltonian isn't exposed to Python yet...
        #H = self.mm.run_as(self.pt, 'Born-Oppenheimer Approximation', sys)
        #H_corr = chemist.qm_operator.Hamiltonian()
        #self.assertEqual(H, H_corr)

    def setUp(self):
        self.mm = ModuleManager()
        nux.load_modules(self.mm)
        self.pt = simde.HamiltonianFromChemicalSystem()
