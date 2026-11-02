from conan import ConanFile
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain
from conan.tools.files import get, copy, rm, rmdir
import os


class LibBigWigRecipe(ConanFile):
    name = "xoos-libbigwig"

    license = "MIT"
    homepage = "https://github.com/dpryan79/libBigWig"
    description = (
        "A C library for reading/parsing local and remote bigWig and bigBed files"
    )

    settings = "os", "compiler", "build_type", "arch"

    options = {"with_curl": [True, False]}
    default_options = {"with_curl": False}

    def requirements(self):
        self.requires("zlib/[>=1.2.11 <2]")
        if self.options.with_curl:
            self.requires("libcurl/[>=7.0]")

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        cmake_deps = CMakeDeps(self)
        cmake_deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["WITH_CURL"] = self.options.with_curl
        tc.variables["BUILD_SHARED_LIBS"] = False
        tc.variables["ENABLE_TESTING"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

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
        self.cpp_info.libs = ["BigWig"]
        self.cpp_info.includedirs = ["include"]
        if not self.options.with_curl:
            self.cpp_info.defines.append("NOCURL")
        if self.settings.os in ["Linux", "FreeBSD"]:
            self.cpp_info.system_libs.append("m")
