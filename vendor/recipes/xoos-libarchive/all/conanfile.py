from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout, CMakeDeps, CMakeToolchain
from conan.tools.files import copy, get, rm, rmdir, download, unzip
import os


class libarchiveRecipe(ConanFile):
    name = "xoos-libarchive"

    license = "BSD-2-Clause"
    homepage = "https://libarchive.org"

    settings = "os", "compiler", "build_type", "arch"

    def layout(self):
        cmake_layout(self, src_folder="src")

    def requirements(self):
        self.requires("zlib-ng/2.2.2", options={"zlib_compat": True})

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        cmake_deps = CMakeDeps(self)
        cmake_deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["ENABLE_NETTLE"] = False
        tc.variables["ENABLE_OPENSSL"] = False
        tc.variables["ENABLE_LIBB2"] = False
        tc.variables["ENABLE_LZ4"] = False
        tc.variables["ENABLE_LZO"] = False
        tc.variables["ENABLE_LZMA"] = False
        tc.variables["ENABLE_ZSTD"] = False
        tc.variables["ENABLE_ZLIB"] = True
        tc.variables["ENABLE_BZip2"] = False
        tc.variables["ENABLE_LIBXML2"] = False
        tc.variables["ENABLE_ICONV"] = False
        tc.variables["ENABLE_EXPAT"] = False
        tc.variables["ENABLE_PCREPOSIX"] = False
        tc.variables["ENABLE_PCRE2POSIX"] = False
        tc.variables["ENABLE_LibGCC"] = False
        tc.variables["ENABLE_CNG"] = False
        tc.variables["ENABLE_ACL"] = False
        tc.variables["ENABLE_TAR"] = False
        tc.variables["ENABLE_CPIO"] = False
        tc.variables["ENABLE_CAT"] = False
        tc.variables["ENABLE_TEST"] = False
        tc.variables["ENABLE_UNZIP"] = False
        tc.variables["ENABLE_WERROR"] = False
        tc.variables["ENABLE_MBEDTLS"] = False
        tc.variables["ENABLE_XATTR"] = False
        tc.variables["BUILD_SHARED_LIBS"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(
            self,
            "COPYING",
            self.source_folder,
            os.path.join(self.package_folder, "licenses"),
        )

        cmake = CMake(self)
        cmake.install()
        rmdir(self, os.path.join(self.package_folder, "lib", "pkgconfig"))
        rmdir(self, os.path.join(self.package_folder, "share"))

    def package_info(self):
        self.cpp_info.libs = ["archive"]
