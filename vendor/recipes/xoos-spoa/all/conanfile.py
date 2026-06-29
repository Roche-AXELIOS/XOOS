from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout, CMakeDeps, CMakeToolchain
from conan.tools.files import get, copy, rm, rmdir, download, unzip
import os


class spoaRecipe(ConanFile):
    name = "xoos-spoa"

    license = "The MIT/Expat License"
    homepage = "https://github.com/rvaser/spoa"

    settings = "os", "compiler", "build_type", "arch"

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        cmake_deps = CMakeDeps(self)
        cmake_deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["spoa_build_tests"] = not self.conf.get(
            "tools.build:skip_test", default=False
        )
        tc.variables["spoa_build_exe"] = False
        tc.variables["spoa_optimize_for_native"] = False
        tc.variables["spoa_optimize_for_portability"] = True
        tc.variables["spoa_use_simde"] = True
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build(target=["all"])

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )

        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["spoa"]
