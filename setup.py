import pybind11
from setuptools import Extension, setup

ext_modules = [
    Extension(
        "partition_function",  # The name of the Python extension
        ["cpp/wrapper.cpp", "cpp/partition_function.cpp"],  # Include both C++ files
        include_dirs=[
            pybind11.get_include(),  # pybind11 include directory
            pybind11.get_include(
                True
            ),  # pybind11 include directory for the "pybind11" submodule
        ],
        language="c++",  # Specify that it's C++ code
        extra_compile_args=["-std=c++11"],  # Use C++11 standard
    ),
]

# Set up the package
setup(
    name="partition_function",  # Package name
    ext_modules=ext_modules,  # Extension modules to be compiled
    zip_safe=False,  # Prevent packaging as a zip file
)
