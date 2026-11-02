import os

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.files import download, copy
from conan.tools.layout import basic_layout


class VincentlaucsbCsvParserRecipe(ConanFile):
    name = "xoos-vincentlaucsb-csv-parser"
    description = "Vince's CSV Parser with simple and intuitive syntax"
    license = "MIT"
    homepage = "https://github.com/vincentlaucsb/csv-parser"
    topics = ("csv", "rfc 4180", "parser", "generator", "header-only")
    package_type = "header-library"
    settings = "os", "arch", "compiler", "build_type"
    no_copy_source = True

    @property
    def _min_cppstd(self):
        # 5.x uses std::string_view, std::from_chars, and optional C++20 features.
        return 17

    def layout(self):
        basic_layout(self, src_folder="src")

    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def source(self):
        # Upstream ships a pre-generated single-header amalgamation as a release
        # asset. The in-repo single_include/csv.hpp is a compatibility shim, so we
        # download the canonical header (and license) directly.
        sources = self.conan_data["sources"][self.version]
        download(self, filename="csv.hpp", **sources["header"])
        download(self, filename="LICENSE", **sources["license"])

    def package(self):
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )
        copy(
            self,
            "csv.hpp",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "include"),
        )

    def package_info(self):
        # Use the default xoos-prefixed CMake package/target names
        # (xoos-vincentlaucsb-csv-parser::xoos-vincentlaucsb-csv-parser) so this
        # vendored package cannot be shadowed by a system-installed
        # vincentlaucsb-csv-parser, matching the other xoos-* vendor recipes.
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        if self.settings.os in ("Linux", "FreeBSD"):
            # Multi-threaded speculative parsing (enabled by default in 5.x) uses
            # std::thread, which requires linking pthread on Linux.
            self.cpp_info.system_libs.append("pthread")

    def package_id(self):
        self.info.clear()
