/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Magdalena Schreter magdalena.schreter@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include <map>

/**
 * @class MarmotGeometryElement
 * @brief Statically-sized geometry base class for all isoparametric Marmot elements.
 *
 * @tparam nDim   Number of spatial dimensions (1, 2 or 3).
 * @tparam nNodes Number of element nodes.
 *
 * Provides shape functions @c N, their natural-coordinate derivatives @c dNdXi,
 * the B-operator @c B, and the element Jacobian.  The element shape is determined
 * automatically from the template parameters.
 *
 * Concrete MarmotElement classes inherit from this class to access the geometric
 * infrastructure without code duplication.
 */
template < int nDim, int nNodes >
class MarmotGeometryElement {

public:
  /*Typedefs*/
  /* static constexpr VoigtSize voigtSize = ( ( ( nDim * nDim ) + nDim ) / 2 ); */
  static constexpr Marmot::ContinuumMechanics::VoigtNotation::VoigtSize
    voigtSize = Marmot::ContinuumMechanics::VoigtNotation::voigtSizeFromDimension( nDim );

  typedef Eigen::Matrix< double, nDim, 1 >                  XiSized;
  typedef Eigen::Matrix< double, nDim * nNodes, 1 >         CoordinateVector;
  typedef Eigen::Matrix< double, nDim, nDim >               JacobianSized;
  typedef Eigen::Matrix< double, 1, nNodes >                NSized;
  typedef Eigen::Matrix< double, nDim, nNodes * nDim >      NBSized;
  typedef Eigen::Matrix< double, nDim, nNodes >             dNdXiSized;
  typedef Eigen::Matrix< double, voigtSize, nNodes * nDim > BSized;
  typedef Eigen::Matrix< double, 4, nNodes * nDim >         BSizedAxisymmetric;

  /*Properties*/
  Eigen::Map< const CoordinateVector >       coordinates;
  const Marmot::FiniteElement::ElementShapes shape;

  /*Methods*/
  MarmotGeometryElement()
    : coordinates( nullptr ), shape( Marmot::FiniteElement::getElementShapeByMetric( nDim, nNodes ) ){};

  std::string getElementShape() const
  {
    using namespace Marmot::FiniteElement;
    static std::map< ElementShapes, std::string > shapes = { { Bar2, "bar2" },
                                                             { Quad4, "quad4" },
                                                             { Quad8, "quad8" },
                                                             { Tetra4, "tetra4" },
                                                             { Tetra10, "tetra10" },
                                                             { Hexa8, "hexa8" },
                                                             { Hexa20, "hexa20" } };

    return shapes[this->shape];
  }

  void assignNodeCoordinates( const double* coords )
  {
    new ( &coordinates ) Eigen::Map< const CoordinateVector >( coords );
  }

  /*Please specialize these functions for each element individially
   *.cpp file.
   *Fully specialized templates are precompiled in marmotMechanics (rather than the unspecialized and
   *partially specialized templates)
   * */
  NSized             N( const XiSized& xi ) const;
  dNdXiSized         dNdXi( const XiSized& xi ) const;
  BSized             B( const dNdXiSized& dNdX ) const;
  BSizedAxisymmetric B_axisymmetric( const dNdXiSized& dNdX, const NSized& N, const XiSized& x_gauss ) const;
  BSized             B_bar( const dNdXiSized& dNdX, const dNdXiSized& dNdX0 ) const;
  BSized             BGreen( const dNdXiSized& dNdX, const JacobianSized& F ) const;

  /*These functions are equal for each element and independent of node number and  nDimension*/
  NBSized NB( const NSized& N ) const { return Marmot::FiniteElement::NB< nDim, nNodes >( N ); }

  JacobianSized Jacobian( const dNdXiSized& dNdXi ) const
  {
    return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdXi, coordinates );
  }

  dNdXiSized dNdX( const dNdXiSized& dNdXi, const JacobianSized& JacobianInverse ) const
  {
    return ( dNdXi.transpose() * JacobianInverse ).transpose();
  }

  JacobianSized F( const dNdXiSized& dNdX, const CoordinateVector& Q ) const
  {
    return Marmot::FiniteElement::Jacobian< nDim, nNodes >( dNdX, Q ) + JacobianSized::Identity();
  }
};
