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
import unittest


class TestLoadModules(unittest.TestCase):

    def test_load_modules(self):
        self.assertGreater(self.mm.size(), 0)

    def setUp(self):
        self.mm = ModuleManager()
        nux.load_modules(self.mm)
