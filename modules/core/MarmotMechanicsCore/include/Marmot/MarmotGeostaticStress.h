#pragma once

#include "Marmot/MarmotJournal.h"
#include <Eigen/Core>
#include <functional>

/**
 * @file MarmotGeostaticStress.h
 * @brief Utilities for evaluating geostatic (gravity-induced) initial stress states.
 *
 * Provides functions for computing the in-situ stress components from a
 * linearly-varying stress distribution defined along a vertical coordinate axis.
 */

namespace Marmot::GeostaticStress {

  /**
   * @brief Computes the three principal geostatic stress components at a given depth.
   *
   * The stress distribution is assumed to vary linearly with the vertical coordinate
   * @p coordinate_y as specified by @p geostaticStressDefintion.
   *
   * @param geostaticStressDefintion  Pointer to the array defining the linear distribution.
   * @param coordinate_y              Vertical coordinate of the evaluation point.
   * @return Tuple of three stress components (vertical, horizontal-1, horizontal-2).
   */
  std::tuple< double, double, double > getGeostaticStressFromLinearDistribution( const double* geostaticStressDefintion,
                                                                                 double        coordinate_y );

} // namespace Marmot::GeostaticStress
