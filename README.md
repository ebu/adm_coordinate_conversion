# ADM Coordinate Conversion Library

A C++ library for converting between polar and Cartesian coordinates used in ADM metadata, as specified by ITU-R BS.2127-0.

## Usage

CMake support has been provided. The library can be built as a static (default) or shared library using the CMake `BUILD_SHARED_LIBS` variable. 

Tests are also included. These can be built using the `BUILD_TESTS` variable, which is off by default. Catch2 is used as the testing framework and is included as a submodule of this repository.

All of the exposed functions and types are declared in `include/adm_coord_conv/adm_coord_conv.hpp` and are fairly self explanatory.

### Examples

To convert only position coordinates, the simpliest solution is;

```
// Polar to Cartesian

double azimuth = 0.5;
double elevation = 0.5;
double distance = 1.0;

CartesianPosition cart = convert(PolarPosition{azimuth, elevation, distance});

std::cout << "X:" << cart.x << "\n";
std::cout << "Y:" << cart.y << "\n";
std::cout << "Z:" << cart.z << "\n";
```

```
// Cartesian to Polar

double x = 0.5;
double y = -0.5;
double z = 0.0;

PolarPosition polar = convert(CartesianPosition{x, y, z});

std::cout << "Azimuth:" << polar.azimuth << "\n";
std::cout << "Elevation:" << polar.elevation << "\n";
std::cout << "Distance:" << polar.distance << "\n";
```

To convert position coordinates along with extent;

```
// Polar to Cartesian with Extent

double azimuth = 0.5;
double elevation = 0.5;
double distance = 1.0;
double extentWidth = 0.5;
double extentHeight = 0.5;
double extentDepth = 0.5;

CartesianSource cart = convert(
    PolarPosition{azimuth, elevation, distance},
    PolarExtent{extentWidth, extentHeight, extentDepth}
);

std::cout << "Position X:" << cart.position.x << "\n";
std::cout << "Position Y:" << cart.position.y << "\n";
std::cout << "Position Z:" << cart.position.z << "\n";
std::cout << "Extent X:" << cart.extent->x << "\n";
std::cout << "Extent Y:" << cart.extent->y << "\n";
std::cout << "Extent Z:" << cart.extent->z << "\n";
```

```
// Cartesian to Polar with Extent

double positionX = 1.0;
double positionY = -1.0;
double positionZ = 0.0;
double extentX = 0.5;
double extentY = 0.5;
double extentZ = 0.5;

PolarSource polar = convert(
    CartesianPosition{positionX, positionY, positionZ},
    CartesianExtent{extentX, extentY, extentZ}
);

std::cout << "Position Azimuth:" << polar.position.x << "\n";
std::cout << "Position Elevation:" << polar.position.y << "\n";
std::cout << "Position Distance:" << polar.position.z << "\n";
std::cout << "Extent Width:" << polar.extent->width << "\n";
std::cout << "Extent Height:" << polar.extent->height << "\n";
std::cout << "Extent Depth:" << polar.extent->depth << "\n";
```

The `convert` method can also be called directly with a `PolarSource` or `CartesianSource` struct, both of which support optional extent data.