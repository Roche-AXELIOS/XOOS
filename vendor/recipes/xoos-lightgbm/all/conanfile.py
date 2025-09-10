import os

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import (
    copy,
    get,
    apply_conandata_patches,
    export_conandata_patches,
)
from conan.tools.env import VirtualBuildEnv
from conan.tools.scm import Version


class LightGBMConan(ConanFile):
    name = "xoos-lightgbm"
    description = (
        "A fast, distributed, high performance gradient boosting "
        "(GBT, GBDT, GBRT, GBM or MART) framework based on decision tree algorithms, "
        "used for ranking, classification and many other machine learning tasks."
    )
    license = "MIT"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://github.com/microsoft/LightGBM"
    topics = ("machine-learning", "boosting")

    settings = "os", "arch", "compiler", "build_type"

    @property
    def _build_testing(self):
        return not self.conf.get("tools.build:skip_test", default=False)

    def export_sources(self):
        export_conandata_patches(self)
        copy(self, "CMakeLists.txt", self.recipe_folder, self.export_sources_folder)

    def layout(self):
        cmake_layout(self, src_folder="src")

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires(
            "fast_double_parser/0.8.0", transitive_headers=True, transitive_libs=True
        )
        self.requires("fmt/11.0.2", transitive_headers=True, transitive_libs=True)
        self.requires(
            "llvm-openmp/18.1.8", transitive_headers=True, transitive_libs=True
        )
        self.test_requires("gtest/1.15.0")
        if Version(self.version) >= "4.3.0":
            self.tool_requires("cmake/[>=3.18 <4]")

    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 11)

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["BUILD_STATIC_LIB"] = True
        tc.cache_variables["USE_DEBUG"] = self.settings.build_type in [
            "Debug",
            "RelWithDebInfo",
        ]
        tc.cache_variables["USE_OPENMP"] = True
        tc.cache_variables["BUILD_CLI"] = False
        #         TODO
        #         tc.cache_variables["BUILD_CPP_TEST"] = self._build_testing

        tc.variables["_MAJOR_VERSION"] = Version(self.version).major
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()
        venv = VirtualBuildEnv(self)
        venv.generate(scope="build")

    def build(self):
        apply_conandata_patches(self)
        cmake = CMake(self)
        cmake.configure(build_script_folder=self.source_path.parent)
        cmake.build()

    #         if self._build_testing:

    def package(self):
        copy(
            self,
            "LICENSE",
            dst=os.path.join(self.package_folder, "licenses"),
            src=self.source_folder,
        )
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "xoos-lightgbm")
        self.cpp_info.set_property("cmake_target_name", "xoos-lightgbm::xoos-lightgbm")

        self.cpp_info.libs = ["_lightgbm"]
        self.cpp_info.system_libs.append("pthread")
