#include <adm_coord_conv/adm_coord_conv.hpp>
#include <array>
#include <cmath>
#include <stdexcept>
#include <tuple>

namespace {
  using namespace adm::coords;

  template<typename Param, typename T>
  float getValueOrZero(T const& element) {
    if (element.template has<Param>()) {
      return element.template get<Param>().get();
    }
    return 0;
  }

  struct AzCartMapping {
    double az;
    double x;
    double y;
    double z;
  };

  CartesianPosition operator*(CartesianPosition const& lhs, double rhs) {
    return { lhs.x * rhs, lhs.y * rhs, lhs.z * rhs };
  }

  CartesianPosition operator*(double lhs, CartesianPosition const& rhs) {
    return rhs * lhs;
  }

  struct ElevationSpec {
    double top{ 30.0 };
    double top_tilde{ 45.0 };
  };

  std::array<AzCartMapping, 5> mapping{ {{0, 0, 1, 0},
                                                {-30, 1, 1, 0},
                                                {-110, 1, -1, 0},
                                                {110, -1, -1, 0},
                                                {30, -1, 1, 0}} };

  constexpr double pi = 3.14159265358979323846;

  constexpr double radians(double degrees) {
    return degrees * pi / 180.0;
  }

  constexpr double degrees(double radians) {
    return radians * 180.0 / pi;
  }

  bool insideAngleRange(double x, double start, double end, double tol = 0.0) {
    while ((end - 360.0) > start) {
      end -= 360.0;
    }
    while (end < start) {
      end += 360.0;
    }
    auto start_tol = start - tol;
    while ((x - 360.0) >= start_tol) {
      x -= 360.0;
    }
    while (x < start_tol) {
      x += 360.0;
    }
    return x <= end + tol;
  }

  std::pair<AzCartMapping, AzCartMapping> findSector(double az) {
    for (auto i = 0u; i != mapping.size(); ++i) {
      auto j = (i + 1) % mapping.size();
      if (insideAngleRange(az, mapping[j].az, mapping[i].az)) {
        return { mapping[i], mapping[j] };
      }
    }
    throw std::runtime_error("azimuth not found in any sector");
  }

  double toAzimuth(double x, double y, double z) {
    return -(degrees(std::atan2(x, y)));
  }

  std::pair<AzCartMapping, AzCartMapping> findCartSector(double az) {
    for (auto i = 0u; i != mapping.size(); ++i) {
      auto j = (i + 1) % mapping.size();
      if (insideAngleRange(az,
        toAzimuth(mapping[j].x, mapping[j].y, mapping[j].z),
        toAzimuth(mapping[i].x, mapping[i].y, mapping[i].z))) {
        return { mapping[i], mapping[j] };
      }
    }
    throw std::runtime_error("azimuth not found in any sector");
  }

  double relativeAngle(double x, double y) {
    while ((y - 360.0) >= x) {
      y -= 360.0;
    }
    while (y < x) {
      y += 360.0;
    }
    return y;
  }

  double mapAzToLinear(double left_az, double right_az, double azimuth) {
    auto mid_az = (left_az + right_az) / 2.0;
    auto az_range = right_az - mid_az;
    auto rel_az = azimuth - mid_az;
    auto gain_r = 0.5 + 0.5 * std::tan(radians(rel_az)) / tan(radians(az_range));
    return std::atan2(gain_r, 1 - gain_r) * (2 / pi);
  }

  double mapLinearToAz(double left_az, double right_az, double x) {
    auto mid_az = (left_az + right_az) / 2.0;
    auto az_range = right_az - mid_az;
    auto gain_l_ = std::cos(x * pi / 2.0);
    auto gain_r_ = std::sin(x * pi / 2.0);
    auto gain_r = gain_r_ / (gain_l_ + gain_r_);
    auto rel_az = degrees(std::atan(2.0 * (gain_r - 0.5) * std::tan(radians(az_range))));
    return mid_az + rel_az;
  }

  std::pair<double, double> calculateXY(double az, double r_xy) {
    AzCartMapping leftSector{};
    AzCartMapping rightSector{};
    std::tie(leftSector, rightSector) = findSector(az);
    auto relAz = relativeAngle(rightSector.az, az);
    auto relLeftAz = relativeAngle(rightSector.az, leftSector.az);
    auto p = mapAzToLinear(relLeftAz, rightSector.az, relAz);
    auto x = r_xy * (leftSector.x + (rightSector.x - leftSector.x) * p);
    auto y = r_xy * (leftSector.y + (rightSector.y - leftSector.y) * p);
    return { x, y };
  }

  std::pair<double, double> calcZ_r_xy(PolarPosition const& polar, ElevationSpec const& elSpec) {
    double z{}, r_xy{};
    if (std::abs(polar.elevation) > elSpec.top) {
      auto el_tilde = elSpec.top_tilde + (90.0 - elSpec.top_tilde) * (std::abs(polar.elevation) - elSpec.top) / (90.0 - elSpec.top);
      z = std::copysign(polar.distance, polar.elevation);
      r_xy = polar.distance * std::tan(radians(90.0 - el_tilde));
    }
    else {
      auto el_tilde = elSpec.top_tilde * (polar.elevation / elSpec.top);
      z = tan(radians(el_tilde)) * polar.distance;
      r_xy = polar.distance;
    }
    return { z, r_xy };
  }

  CartesianPosition pointPolarToCart(PolarPosition const& polar) {
    CartesianPosition cart{};
    double r_xy{};
    std::tie(cart.z, r_xy) = calcZ_r_xy(polar, ElevationSpec());
    std::tie(cart.x, cart.y) = calculateXY(polar.azimuth, r_xy);
    return cart;
  }

  bool xAndYNearZero(CartesianPosition const& cart, double epsilon) {
    return (std::abs(cart.x) < epsilon && std::abs(cart.y) < epsilon);
  }

  bool zNearZero(CartesianPosition const& cart, double epsilon) {
    return std::abs(cart.z) < epsilon;
  }

  std::optional<PolarPosition> nearZeroConversion(CartesianPosition const& cart) {
    double const epsilon{ 1e-10 };
    if (xAndYNearZero(cart, epsilon)) {
      if (zNearZero(cart, epsilon)) {
        return { {0, 0, 0} };
      }
      else {
        return { {0.0, std::copysign(90, cart.z), std::abs(cart.z)} };
      }
    }
    return {};
  }

  std::array<double, 4> invert(std::array<double, 4> const& matrix2by2) {
    auto det = (matrix2by2[0] * matrix2by2[3] - matrix2by2[1] * matrix2by2[2]);
    auto mul = 1.0 / det;
    auto inv = std::array<double, 4> {matrix2by2[3] * mul, -matrix2by2[1] * mul,
      -matrix2by2[2] * mul, matrix2by2[0] * mul};
    return inv;
  }

  std::array<double, 2> dotProduct(std::array<double, 2> const& rowVec, std::array<double, 4> const& matrix2by2) {
    return { rowVec[0] * matrix2by2[0] + rowVec[1] * matrix2by2[2],
            rowVec[0] * matrix2by2[1] + rowVec[1] * matrix2by2[3] };
  }

  // didn't want to add a dependency to do one linear algebra calc
  std::pair<double, double> calculate_g_lr(CartesianPosition const& cart, std::pair<AzCartMapping, AzCartMapping> const& sectors) {
    auto const& left = sectors.first;
    auto const& right = sectors.second;
    auto inverse = invert({ left.x, left.y, right.x, right.y });
    auto dot = dotProduct({ cart.x, cart.y }, inverse);
    return { dot[0], dot[1] };
  }

  PolarPosition convert(CartesianPosition const& cart,
    ElevationSpec elSpec) {
    auto sectors = findCartSector(toAzimuth(cart.x, cart.y, 0));
    auto const& left = sectors.first;
    auto const& right = sectors.second;
    auto rel_left_az = relativeAngle(right.az, left.az);
    auto g_lr = calculate_g_lr(cart, sectors);
    auto r_xy = g_lr.first + g_lr.second;
    PolarPosition polar{};
    auto relAz = mapLinearToAz(rel_left_az, right.az, g_lr.second / r_xy);
    polar.azimuth = relativeAngle(-180, relAz);
    auto el_tilde = degrees(std::atan(cart.z / r_xy));

    if (std::abs(el_tilde) > elSpec.top_tilde) {
      auto abs_el = elSpec.top + ((90.0 - elSpec.top) * (std::abs(el_tilde) - elSpec.top_tilde))
        / (90.0 - elSpec.top_tilde);
      polar.elevation = std::copysign(abs_el, el_tilde);
      polar.distance = std::abs(cart.z);
    }
    else {
      polar.elevation = elSpec.top * el_tilde / elSpec.top_tilde;
      polar.distance = r_xy;
    }
    return polar;
  }

  PolarPosition pointCartToPolar(CartesianPosition const& cart) {
    auto converted = nearZeroConversion(cart);
    if (converted) {
      return *converted;
    }
    else {
      return convert(cart, ElevationSpec());
    }
  }


  // Converts a polar extent value into a cartesian extent assuming position in directly in front of listener
  // and a radius of 1. See BS.2127, section 10.2
  CartesianExtent whd2xyz(PolarExtent const& polarExtent) {
    double x_size_width{}, y_size_width{}, z_size_height{}, y_size_height{}, y_size_depth{};
    if (polarExtent.width < 180.0) {
      x_size_width = std::sin(radians(polarExtent.width / 2.0));
    }
    else {
      x_size_width = 1.0;
    }
    y_size_width = (1.0 - std::cos(radians(polarExtent.width / 2.0))) / 2.0;

    if (polarExtent.height < 180.0) {
      z_size_height = std::sin(radians(polarExtent.height / 2.0));
    }
    else {
      z_size_height = 1.0;
    }
    y_size_height = (1.0 - std::cos(radians(polarExtent.height / 2.0))) / 2.0;
    y_size_depth = polarExtent.depth;
    return { x_size_width, std::max(y_size_width, std::max(y_size_height, y_size_depth)), z_size_height };
  }

  PolarExtent xyz2whd(CartesianExtent const& cartExtent) {
    PolarExtent polarExtent{};
    auto width_from_sx = 2.0 * degrees(std::asin(cartExtent.x));
    auto width_from_sy = 2.0 * degrees(std::acos(1.0 - (2 * cartExtent.y)));
    polarExtent.width = width_from_sx + cartExtent.x * std::max(width_from_sy - width_from_sx, 0.0);

    auto height_from_sz = 2.0 * degrees(std::asin(cartExtent.z));
    auto height_from_sy = 2.0 * degrees(std::acos(1.0 - (2.0 * cartExtent.y)));
    polarExtent.height = height_from_sz + cartExtent.z * std::max(height_from_sy - height_from_sz, 0.0);

    auto equiv_y = whd2xyz({ polarExtent.width, polarExtent.height, 0.0 }).y;
    polarExtent.depth = std::max(0.0, cartExtent.y - equiv_y);
    return polarExtent;
  }

  double euclidianNorm(double a, double b, double c) {
    return std::sqrt(a * a + b * b + c * c);
  }

  CartesianPosition toCart(PolarPosition const& polar) {
    CartesianPosition cart{};
    cart.x = std::sin(-pi * polar.azimuth / 180.0) * std::cos(pi * polar.elevation / 180.0) * polar.distance;
    cart.y = std::cos(-pi * polar.azimuth / 180.0) * std::cos(pi * polar.elevation / 180.0) * polar.distance;
    cart.z = std::sin(pi * polar.elevation / 180.0) * polar.distance;
    return cart;
  }

  // returns a rotation matrix that maps a forward vector to a given azimuth an elevation
  // See BS.2127, section 10.2
  std::array<CartesianPosition, 3> localCoordinateSystem(double az, double el) {
    return { toCart({az - 90.0, 0, 1}),
            toCart({az, el, 1}),
            toCart({az, el + 90.0, 1}) };
  }

  std::array<CartesianPosition, 3> transpose(std::array<CartesianPosition, 3> const& mat) {
    return{ CartesianPosition{mat[0].x, mat[1].x, mat[2].x},
           CartesianPosition{mat[0].y, mat[1].y, mat[2].y},
           CartesianPosition{mat[0].z, mat[1].z, mat[2].z} };
  }

  std::array<CartesianPosition, 3> calcMPolarToCart(CartesianExtent const& forwardExtent, PolarPosition const& position) {
    std::array<CartesianPosition, 3> M{};
    auto rotation = localCoordinateSystem(position.azimuth, position.elevation);
    M[0] = { forwardExtent.x * rotation[0] };
    M[1] = { forwardExtent.y * rotation[1] };
    M[2] = { forwardExtent.z * rotation[2] };
    return M;
  }

  std::array<CartesianPosition, 3> calcMCartToPolar(CartesianExtent const& forwardExtent, PolarPosition const& position) {
    std::array<CartesianPosition, 3> M{};
    auto rotation = transpose(localCoordinateSystem(position.azimuth, position.elevation));
    M[0] = { forwardExtent.x * rotation[0] };
    M[1] = { forwardExtent.y * rotation[1] };
    M[2] = { forwardExtent.z * rotation[2] };
    return M;
  }

  CartesianExtent polarToCartExtent(PolarPosition const& polarPosition, PolarExtent const& polarExtent) {
    auto cartesianForwardExtent = whd2xyz(polarExtent);
    auto M = calcMPolarToCart(cartesianForwardExtent, polarPosition);
    return { euclidianNorm(M[0].x, M[1].x, M[2].x),
            euclidianNorm(M[0].y, M[1].y, M[2].y),
            euclidianNorm(M[0].z, M[1].z, M[2].z) };
  }

  PolarExtent cartToPolarExtent(PolarPosition const& polarPosition, CartesianExtent const& cartExtent) {
    auto M = calcMCartToPolar(cartExtent, polarPosition);
    auto cartForwardExtent = CartesianExtent{ euclidianNorm(M[0].x, M[1].x, M[2].x),
                                             euclidianNorm(M[0].y, M[1].y, M[2].y),
                                             euclidianNorm(M[0].z, M[1].z, M[2].z) };
    return xyz2whd(cartForwardExtent);
  }
}

namespace adm {
  namespace coords {

    ADMCOORDCONV_EXPORT CartesianPosition convert(PolarPosition const& input)
    {
      return pointPolarToCart(input);
    }

    ADMCOORDCONV_EXPORT CartesianSource convert(PolarSource const& input)
    {
      if (input.extent.has_value()) {
        return convert(input.position, input.extent.value());
      }
      return CartesianSource{ convert(input.position) };
    }

    ADMCOORDCONV_EXPORT CartesianSource convert(PolarPosition const& position, PolarExtent const& extent) {
      return CartesianSource{ 
        pointPolarToCart(position), 
        polarToCartExtent(position, extent) 
      };
    }

    ADMCOORDCONV_EXPORT PolarPosition convert(CartesianPosition const& input)
    {
      return pointCartToPolar(input);
    }

    ADMCOORDCONV_EXPORT PolarSource convert(CartesianSource const& input)
    {
      if (input.extent.has_value()) {
        return convert(input.position, input.extent.value());
      }
      return PolarSource{ convert(input.position) };
    }

    ADMCOORDCONV_EXPORT PolarSource convert(CartesianPosition const& position, CartesianExtent const& extent) {
      PolarSource ret;
      ret.position = pointCartToPolar(position);
      ret.extent = cartToPolarExtent(ret.position, extent);
      return ret;
    }

  }
}