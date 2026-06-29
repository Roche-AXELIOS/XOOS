from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout, CMakeDeps, CMakeToolchain
from conan.tools.files import copy, get, rmdir
import os


class kseqppRecipe(ConanFile):
    name = "xoos-kseqpp"

    license = "MIT"
    homepage = "https://github.com/cartoonist/kseqpp"

    settings = "os", "compiler", "build_type", "arch"

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        self.requires("zlib-ng/2.2.2", options={"zlib_compat": True})
        self.requires("bzip2/1.0.8", options={"build_executable": False})

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        cmake_deps = CMakeDeps(self)
        cmake_deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["BUILD_TESTING"] = not self.conf.get(
            "tools.build:skip_test", default=False
        )
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()

        targets_to_build = ["all"]
        if not self.conf.get("tools.build:skip_test", default=False):
            targets_to_build.append("test")

        cmake.build(target=targets_to_build)

    def package(self):
        copy(
            self,
            "*LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )

        cmake = CMake(self)
        cmake.install()

        rmdir(self, os.path.join(self.package_folder, "lib"))

    def package_info(self):
        self.cpp_info.system_libs = ["pthread"]
