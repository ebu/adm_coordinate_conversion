// main.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <array>
#include <random>

#include <adm_coord_conv/adm_coord_conv.hpp>
#include <catch2/catch_all.hpp> // Include the Catch2 header

using Catch::Approx;
using namespace adm::coords;

namespace {
  std::vector<std::vector<double>> getRandomCoordinates(double min, double max, int count) {
    std::vector<std::vector<double>> coordinates;
    coordinates.reserve(count);
    std::random_device r;
    auto seed = r();
    INFO("Seed: " + std::to_string(seed));
    std::default_random_engine engine{ seed };
    std::uniform_real_distribution<double> dist{ min, max };
    for (auto i = 0; i != count; ++i) {
      coordinates.push_back({ dist(engine), dist(engine), dist(engine) });
    }
    return coordinates;
  }
}

TEST_CASE("Test conversion corners") {

  using Azimuth = double;
  using Elevation = double;
  using Distance = double;
  using X = double;
  using Y = double;
  using Z = double;

  std::array<std::tuple<Elevation, Z>, 3> el_z{
          std::make_tuple(-30, -1),
          std::make_tuple(0, 0),
          std::make_tuple(30,1)
  };

  std::array<std::tuple<Azimuth, X, Y>, 5> az_x_y{
          std::make_tuple(0, 0, 1),
          std::make_tuple(-30, 1, 1),
          std::make_tuple(30, -1, 1),
          std::make_tuple(-110, 1, -1),
          std::make_tuple(110, -1, -1)
  };

  std::array<Distance, 3> distances{
          0.5, 1.0, 2.0
  };

  {
    Elevation el;
    Z z;
    for (auto const& el_z_pair : el_z) {
      std::tie(el, z) = el_z_pair;
      Azimuth az;
      X x;
      Y y;
      for (auto const& az_x_y_tuple : az_x_y) {
        std::tie(az, x, y) = az_x_y_tuple;
        for (auto const& d : distances) {
          if (el == Approx(0) || az != Approx(0)) {

            auto cartRet = convert(PolarPosition{ az, el, d });

            double expectedCartX = x * d;
            double expectedCartY = y * d;
            double expectedCartZ = z * d;

            CHECK(cartRet.x == Approx(expectedCartX).margin(1e-10));
            CHECK(cartRet.y == Approx(expectedCartY).margin(1e-10));
            CHECK(cartRet.z == Approx(expectedCartZ).margin(1e-10));

            auto polarRet = convert(CartesianPosition{ cartRet.x, cartRet.y, cartRet.z });

            CHECK(polarRet.azimuth == Approx(az).margin(1e-10));
            CHECK(polarRet.elevation == Approx(el).margin(1e-10));
            CHECK(polarRet.distance == Approx(d).margin(1e-10));
          }
        }
      }
    }
  }
}

TEST_CASE("Test conversion poles") {
  std::array<double, 2> minusPlus{ -1.0, 1.0 };
  std::array<double, 3> distances{ 0.5, 1.0, 2.0 };

  double az = 0.f;
  for (auto sign : minusPlus) {
    double el = sign * 90.0f;
    for (auto d : distances) {

      auto cartRet = convert(PolarSource{ {az, el, d} });

      double expectedCartX = 0.;
      double expectedCartY = 0.;
      double expectedCartZ = sign * d;

      CHECK(cartRet.position.x == Approx(expectedCartX).margin(1e-10));
      CHECK(cartRet.position.y == Approx(expectedCartY).margin(1e-10));
      CHECK(cartRet.position.z == Approx(expectedCartZ).margin(1e-10));

      auto polarRet = convert(CartesianSource{ {cartRet.position.x, cartRet.position.y, cartRet.position.z} });

      CHECK(polarRet.position.azimuth == Approx(az).margin(1e-10));
      CHECK(polarRet.position.elevation == Approx(el).margin(1e-10));
      CHECK(polarRet.position.distance == Approx(d).margin(1e-10));
    }
  }
}

TEST_CASE("Test conversion centre") {
  std::array<double, 3> azimuths{ -90.0, 0.0, 90.0 };
  std::array<double, 3> elevations{ -90.0, 0.0, 90.0 };

  double d = 0.0;
  for (auto az : azimuths) {
    for (auto el : elevations) {
      auto cartRet = convert(PolarSource{ {az, el, d} });

      CHECK(cartRet.position.x == Approx(0.).margin(1e-10));
      CHECK(cartRet.position.y == Approx(0.).margin(1e-10));
      CHECK(cartRet.position.z == Approx(0.).margin(1e-10));
    }
  }

  auto polarRet = convert(CartesianSource{ {0., 0., 0.} });
  CHECK(polarRet.position.distance == Approx(0.0));
}

TEST_CASE("Test conversion reversible") {
  auto cartPositions = getRandomCoordinates(-2, 1, 1000);
  for (auto const& pos : cartPositions) {
    INFO("(x, y, z): " << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")");
    PolarSource polar;
    REQUIRE_NOTHROW(polar = convert(CartesianSource{ {pos[0], pos[1], pos[2]} }));
    INFO("(az, el, d): " << "(" << polar.position.azimuth << ", " << polar.position.elevation << ", " << polar.position.distance << ")");
    CartesianSource cart;
    REQUIRE_NOTHROW(cart = convert(polar));
    CHECK(cart.position.x == Approx(pos[0]).margin(1e-10));
    CHECK(cart.position.y == Approx(pos[1]).margin(1e-10));
    CHECK(cart.position.z == Approx(pos[2]).margin(1e-10));
  }
}

TEST_CASE("Test cartesian to polar extent conversion consistent with python EAR") {
  std::ifstream inputFile("test_data-pos_extent_cart2polar.txt");
  REQUIRE(inputFile.good());

  std::string line;
  while (std::getline(inputFile, line)) {
    std::vector<double> row;
    std::stringstream ss(line);
    std::string value;

    // Split the line by commas and convert to doubles
    while (std::getline(ss, value, ',')) {
      try {
        row.push_back(std::stod(value)); // Convert to double and add to row
      }
      catch (const std::invalid_argument& e) {
        std::cerr << "Invalid number found: " << value << std::endl;
      }
    }

    auto polar = convert(CartesianPosition{row[0], row[1], row[2]}, CartesianExtent{row[3], row[4], row[5]});

    CHECK(polar.position.azimuth == Approx(row[6]).margin(10e-5));
    CHECK(polar.position.elevation == Approx(row[7]).margin(10e-5));
    CHECK(polar.position.distance == Approx(row[8]).margin(10e-5));
    CHECK(polar.extent.has_value());
    if (polar.extent.has_value()) {
      CHECK(polar.extent->width == Approx(row[9]).margin(10e-5));
      CHECK(polar.extent->height == Approx(row[10]).margin(10e-5));
      CHECK(polar.extent->depth == Approx(row[11]).margin(10e-5));
    }
  }

  inputFile.close();
}