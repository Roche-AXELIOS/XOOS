from conan import ConanFile
from conan.tools.gnu import Autotools, AutotoolsToolchain, AutotoolsDeps
from conan.tools.files import get, copy, rm, rmdir, download, unzip
import os


class htslibRecipe(ConanFile):
    name = "xoos-htslib"

    license = "The MIT/Expat License"
    homepage = "https://github.com/samtools/htslib"

    settings = "os", "compiler", "build_type", "arch"

    options = {"with_zlibng": [True, False]}

    default_options = {"with_zlibng": False}

    def requirements(self):
        self.requires("libdeflate/1.22")
        self.requires(
            "bzip2/1.0.8", options={"shared": False, "build_executable": False}
        )
        if self.options.with_zlibng:
            self.requires(
                "zlib-ng/2.2.2", options={"zlib_compat": True, "shared": False}
            )
        else:
            self.requires("zlib/[>=1.2.11 <2]")

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def generate(self):
        tc = AutotoolsToolchain(self)
        tc.configure_args.extend(
            [
                "--disable-lzma",
                "--disable-libcurl",
                "--with-libdeflate",
            ]
        )
        tc.generate()

        ad = AutotoolsDeps(self)
        ad.generate()

    def build(self):
        autotools = Autotools(self)
        autotools.autoreconf()
        autotools.configure()
        autotools.make()
        if not self.conf.get("tools.build:skip_test", default=False):
            autotools.make(target="test-shlib-exports")
            autotools.make(target="test")

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )

        autotools = Autotools(self)
        autotools.install()
        rmdir(self, os.path.join(self.package_folder, "share"))
        rmdir(self, os.path.join(self.package_folder, "bin"))
        rmdir(self, os.path.join(self.package_folder, "lib/pkgconfig"))
        rm(self, "*.so", os.path.join(self.package_folder, "lib"))
        rm(self, "*.so.*", os.path.join(self.package_folder, "lib"))
        rm(self, "*.dylib", os.path.join(self.package_folder, "lib"))
        rm(self, "*.dylib.*", os.path.join(self.package_folder, "lib"))

    def package_info(self):
        self.cpp_info.libs = ["hts"]
