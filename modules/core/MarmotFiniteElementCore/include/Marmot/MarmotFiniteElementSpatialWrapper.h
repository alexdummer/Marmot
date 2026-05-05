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
#include "Eigen/Sparse"
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include <functional>
#include <memory>

/**
 * @class MarmotElementSpatialWrapper
 * @brief Wrapper that embeds a lower-dimensional child element (e.g.\ a truss) into a
 *        higher-dimensional ambient space (2-D or 3-D).
 *
 * The projection transformation is constructed automatically from the supplied
 * nodal coordinates.  The child element is created via a user-provided factory
 * functor so that the wrapper remains independent of the concrete child type.
 */
class MarmotElementSpatialWrapper : public MarmotElement {

public:
  const int                                 nDim;
  const int                                 nDimChild;
  const int                                 nNodes;
  const int                                 nRhsChild;
  const Eigen::Map< const Eigen::VectorXi > rhsIndicesToBeProjected;
  const int                                 projectedSize, unprojectedSize;

  std::unique_ptr< MarmotElement > childElement;
  Eigen::MatrixXd                  T;
  Eigen::MatrixXd                  P;
  Eigen::MatrixXd                  projectedCoordinates;

  MarmotElementSpatialWrapper( int                              nDim,
                               int                              nChildDim,
                               int                              nNodes,
                               int                              sizeRhsChild,
                               const int                        rhsIndicesToBeWrapped_[],
                               int                              nRhsIndicesToBeWrapped,
                               std::unique_ptr< MarmotElement > childElement );

  int getNumberOfRequiredStateVars();

  std::vector< std::vector< std::string > > getNodeFields();

  std::vector< int > getDofIndicesPermutationPattern();

  int getNNodes();

  int getNSpatialDimensions();

  int getNDofPerElement();

  std::string getElementShape();

  void assignStateVars( double* stateVars, int nStateVars );

  void assignProperty( const ElementProperties& property );

  void assignProperty( const MarmotMaterialSection& property );

  void assignNodeCoordinates( const double* coordinates );

  void initializeYourself();

  void computeYourself( const double* QTotal,
                        const double* dQ,
                        double*       Pe,
                        double*       Ke,
                        const double* time,
                        double        dT,
                        double&       pNewdT );

  void setInitialConditions( StateTypes state, const double* values );

  void computeDistributedLoad( DistributedLoadTypes loadType,
                               double*              P,
                               double*              K,
                               int                  elementFace,
                               const double*        load,
                               const double*        QTotal,
                               const double*        time,
                               double               dT );

  void computeBodyForce( double*       P,
                         double*       K,
                         const double* load,
                         const double* QTotal,
                         const double* time,
                         double        dT );

  StateView getStateView( const std::string& stateName, int quadraturePoint );

  std::vector< double > getCoordinatesAtCenter();

  std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

  int getNumberOfQuadraturePoints();
};
