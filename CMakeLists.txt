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

cmake_minimum_required(VERSION 3.14)

# Downloads common CMake modules used throughout NWChemEx

#Sets the version to whatever git thinks it is
#include(get_version_from_git)
#get_version_from_git(scf_version "${CMAKE_CURRENT_LIST_DIR}")
project(nux VERSION "1.0.0" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
include(cmake/get_nwx_cmake.cmake)

include(get_cmaize)
include(nwx_cxx_api_docs)

set(NUX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(NUX_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(NUX_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

nwx_cxx_api_docs("${NUX_SOURCE_DIR}" "${NUX_INCLUDE_DIR}")

### Options ###
cmaize_option_list(
    BUILD_TESTING OFF "Should we build the tests?"
    BUILD_PYBIND11_PYBINDINGS ON "Build Python bindings with pybind11?"
    INTEGRATION_TESTING OFF "Should we build the integration tests?"
)

set(DEPENDENCIES simde)

foreach(dependency_i ${DEPENDENCIES})
    include(get_${dependency_i})
endforeach()

cmaize_add_library(
    ${PROJECT_NAME}
    SOURCE_DIR "${NUX_SOURCE_DIR}/cxx/nux"
    INCLUDE_DIRS "${NUX_INCLUDE_DIR}/nux"
    DEPENDS "${DEPENDENCIES}"
)

include(nwx_pybind11)
nwx_add_pybind11_module(
    ${PROJECT_NAME}
    SOURCE_DIR "${NUX_SOURCE_DIR}/python"
    DEPENDS "${PROJECT_NAME}"
)

if("${BUILD_TESTING}")
    include(CTest)
    set(PYTHON_TEST_DIR "${NUX_TESTS_DIR}/python")
    set(CXX_TEST_DIR "${NUX_TESTS_DIR}/cxx")

    include(get_catch2)
    cmaize_add_tests(
        unit_test_${PROJECT_NAME}
        SOURCE_DIR "${CXX_TEST_DIR}/unit_tests"
        INCLUDE_DIRS "${NUX_SOURCE_DIR}/cxx/nux"
        DEPENDS Catch2 ${PROJECT_NAME}
    )
    nwx_pybind11_tests(
        py_unit_test_${PROJECT_NAME}
        "${PYTHON_TEST_DIR}/unit_tests/run_unit_tests.py"
        SUBMODULES nux simde chemist pluginplay parallelzone
    )

    if("${INTEGRATION_TESTING}")
        include(get_nwchemex)
        nwx_pybind11_tests(
            py_integration_test_${PROJECT_NAME}
            "${PYTHON_TEST_DIR}/integration_tests/run_integration_tests.py"
            SUBMODULES nwchemex nux chemcache simde chemist pluginplay
                       parallelzone scf friendzone
        )
    endif()
endif()

cmaize_add_package(${PROJECT_NAME} NAMESPACE nwx::)
