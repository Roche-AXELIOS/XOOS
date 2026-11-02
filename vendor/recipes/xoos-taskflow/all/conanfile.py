from conan import ConanFile
from conan.tools.files import copy, get
import os


class TaskflowRecipe(ConanFile):
    name = "xoos-taskflow"
    description = (
        "A fast C++ header-only library to help you quickly write parallel "
        "programs with complex task dependencies."
    )
    license = "MIT"
    homepage = "https://github.com/taskflow/taskflow"
    package_type = "header-library"
    settings = "os", "compiler", "build_type", "arch"
    no_copy_source = True

    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def package(self):
        copy(
            self,
            "LICENSE",
            self.source_folder,
            os.path.join(self.package_folder, "licenses"),
        )
        copy(
            self,
            "*",
            os.path.join(self.source_folder, "taskflow"),
            os.path.join(self.package_folder, "include", "taskflow"),
        )

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        if self.settings.os in ["Linux", "FreeBSD"]:
            self.cpp_info.system_libs.append("pthread")

    def package_id(self):
        self.info.clear()
