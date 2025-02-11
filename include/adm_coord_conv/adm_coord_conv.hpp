#pragma once

#include "adm_coord_conv_export.hpp"
#include <optional>

namespace adm {
  namespace coords {

    struct ADMCOORDCONV_EXPORT CartesianPosition {
      double x;
      double y;
      double z;
    };

    struct ADMCOORDCONV_EXPORT CartesianExtent {
      double x;
      double y;
      double z;
    };

    struct ADMCOORDCONV_EXPORT CartesianSource {
      CartesianPosition position;
      std::optional<CartesianExtent> extent;
    };

    struct ADMCOORDCONV_EXPORT PolarPosition {
      double azimuth;
      double elevation;
      double distance = 1.0;
    };

    struct ADMCOORDCONV_EXPORT PolarExtent {
      double width;
      double height;
      double depth;
    };

    struct ADMCOORDCONV_EXPORT PolarSource {
      PolarPosition position;
      std::optional<PolarExtent> extent;
    };

    ADMCOORDCONV_EXPORT CartesianPosition convert(PolarPosition const& position);
    ADMCOORDCONV_EXPORT CartesianSource convert(PolarSource const& source);
    ADMCOORDCONV_EXPORT CartesianSource convert(PolarPosition const& position, PolarExtent const& extent);

    ADMCOORDCONV_EXPORT PolarPosition convert(CartesianPosition const& position);
    ADMCOORDCONV_EXPORT PolarSource convert(CartesianSource const& source);
    ADMCOORDCONV_EXPORT PolarSource convert(CartesianPosition const& position, CartesianExtent const& extent);

  }
}
