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
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElastic.h"
#include "Marmot/MarmotMaterialGeneralGradientEnhancedHypoElasticFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  template < int nDim, int nNodes, int nNonlocalVariables = 1, int nNonLocalNodes = nNodes >
  class GeneralGradientEnhancedDisplacementFiniteElement : public MarmotElement {

  public:
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim;
    static constexpr int nDofPerNodeK = nNonlocalVariables;

    static constexpr int sizeLoadVector = nNodes * nDofPerNodeU + nNonLocalNodes * nDofPerNodeK;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int sizeDoFU = nNodes * nDofPerNodeU;
    static constexpr int sizeDoFK = nNonLocalNodes * nDofPerNodeK;

    MarmotGeometryElement< nDim, nNodes >         localGeometryElement;
    MarmotGeometryElement< nDim, nNonLocalNodes > nonLocalGeometryElement;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using NSized                = typename ParentGeometryElement::NSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using BSized                = typename ParentGeometryElement::BSized;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using RhsSized              = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector          = Matrix< double, sizeDoFU, 1 >;
    using KSizedVector          = Matrix< double, sizeDoFK, 1 >;
    using CSized                = Matrix< double, ParentGeometryElement::voigtSize, ParentGeometryElement::voigtSize >;
    using Voigt                 = Matrix< double, ParentGeometryElement::voigtSize, 1 >;

    using ParentGeometryElementK = MarmotGeometryElement< nDim, nNonLocalNodes >;
    using JacobianSizedK         = typename ParentGeometryElementK::JacobianSized;
    using NSizedK                = typename ParentGeometryElementK::NSized;
    using dNdXiSizedK            = typename ParentGeometryElementK::dNdXiSized;
    using BSizedK                = typename ParentGeometryElementK::BSized;

    Map< const VectorXd > elementProperties;
    const int             elLabel;

    const SectionType sectionType;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double      detJ;
      double      J0xW;
      NSized      N;
      dNdXiSized  dNdX;
      BSized      B;
      NSizedK     N_K;
      dNdXiSizedK dNdX_K;
      BSizedK     B_K;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 6 },
          { .name = "strain", .length = 6 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        mVector6d                     stress;
        mVector6d                     strain;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() ){};
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables > > material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      };

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      };

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
        /* material->assignStateVars( managedStateVars->materialStateVars.data(), */
        /*                            managedStateVars->materialStateVars.size() ); */
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          J0xW( 0.0 ),
          dNdX( dNdXiSized::Zero() ),
          B( BSized::Zero() ),
          dNdX_K( dNdXiSizedK::Zero() ),
          B_K( BSizedK::Zero() ){};
    };

    std::vector< QuadraturePoint > qps;

    GeneralGradientEnhancedDisplacementFiniteElement( int                                         elementID,
                                                      FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                      SectionType                                 sectionType );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return localGeometryElement.getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const ElementProperties& marmotElementProperty );

    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 const double*                       time,
                                 double                              dT );

    void computeBodyForce( double* P,
                           double* K,

                           const double* load,
                           const double* QTotal,
                           const double* time,
                           double        dT );

    void computeYourself( const double* QTotal,
                          const double* dQ,
                          double*       Pe,
                          double*       Ke,
                          const double* time,
                          double        dT,
                          double&       pNewdT );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( stateName == "sdv" ) {
        std::cout << __PRETTY_FUNCTION__ << " on 'sdv' is discouraged and deprecated, please use precise state name";
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }
      else {
        return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
      }
    }

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    GeneralGradientEnhancedDisplacementFiniteElement( int                                         elementID,
                                                      FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                      SectionType                                 sectionType )
    : elementProperties( nullptr, 0 ), elLabel( elementID ), sectionType( sectionType )
  {
    for ( const auto& qpInfo :
          FiniteElement::Quadrature::getGaussPointInfo( localGeometryElement.shape, integrationType ) ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  int GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties )
      Map< const VectorXd >( elementPropertiesInfo.elementProperties, elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables > >(
        MarmotLibrary::MarmotMaterialGeneralGradientEnhancedHypoElasticFactory< nNonlocalVariables >::
          createMaterial( section.materialName, section.materialProperties, section.nMaterialProperties, elLabel ) );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< std::vector< std::string > > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        if ( i < nNonLocalNodes ) {
          if constexpr ( nNonlocalVariables == 6 )
            nodeFields[i].push_back( "strain symmetric" );
          else {
            nodeFields[i].push_back( "nonlocal damage" );
            for ( int j = 1; j < nNonlocalVariables; j++ )
              nodeFields[i].push_back( "nonlocal damage " + to_string( j + 1 ) );
          }
        }
      }
    return nodeFields;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< int > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * nDim + nNonlocalVariables * ( i < nNonLocalNodes ? i : nNonLocalNodes ) +
                                        j );
      for ( int j = 0; j < nNonlocalVariables; j++ )
        for ( int i = 0; i < nNonLocalNodes; i++ )
          permutationPattern.push_back( i * ( nDim + nNonlocalVariables ) + nDim + j );
    }
    return permutationPattern;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    assignNodeCoordinates( const double* coordinates )
  {
    localGeometryElement.assignNodeCoordinates( coordinates );
    nonLocalGeometryElement.assignNodeCoordinates( coordinates );
    // This assumes that the corner nodes (vertices) are listed before the mid-edge nodes!
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {
      const auto          dNdXi = localGeometryElement.dNdXi( qp.xi );
      const JacobianSized J     = localGeometryElement.Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      qp.detJ                   = J.determinant();
      qp.N                      = localGeometryElement.N( qp.xi );
      qp.dNdX                   = localGeometryElement.dNdX( dNdXi, JInv );
      qp.B                      = localGeometryElement.B( qp.dNdX );

      const auto           dNdXi_K = nonLocalGeometryElement.dNdXi( qp.xi );
      const JacobianSizedK J_K     = nonLocalGeometryElement.Jacobian( dNdXi_K );
      const JacobianSizedK JInv_K  = J_K.inverse();
      qp.N_K                       = nonLocalGeometryElement.N( qp.xi );
      qp.dNdX_K                    = nonLocalGeometryElement.dNdX( dNdXi_K, JInv_K );
      qp.B_K                       = nonLocalGeometryElement.B( qp.dNdX_K );

      if constexpr ( nDim == 3 ) {
        qp.J0xW = qp.weight * qp.detJ;
      }
      if constexpr ( nDim == 2 ) {
        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * qp.detJ * thickness;
      }
      if constexpr ( nDim == 1 ) {
        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * qp.detJ * crossSection;
      }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeYourself( const double* QTotal_,
                     const double* dQ_,
                     double*       Pe_,
                     double*       Ke_,
                     const double* time,
                     double        dT,
                     double&       pNewDT )
  {

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    // incremental nodal displacements Q and nodal internal parameters qK
    const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
    const Ref< const KSizedVector > dQK( dQ.tail( sizeDoFK ) );
    const Ref< const USizedVector > qU( QTotal.head( sizeDoFU ) );
    const Ref< const KSizedVector > qK( QTotal.tail( sizeDoFK ) );

    // substiffness matrices
    // [k_qq    k_qQK
    //  k_QKQ   k_QKQK ]
    Ref< Matrix< double, sizeDoFU, sizeDoFU > > kUU( Ke.topLeftCorner( sizeDoFU, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFU, sizeDoFK > > kUK( Ke.topRightCorner( sizeDoFU, sizeDoFK ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFU > > kKU( Ke.bottomLeftCorner( sizeDoFK, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFK > > kKK( Ke.bottomRightCorner( sizeDoFK, sizeDoFK ) );

    // Righthandside
    Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
    Ref< KSizedVector > fK( Pe.tail( sizeDoFK ) );

    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Voigt  S, dE, dSdK, dK_Local_dDE;
    CSized C;
    using response  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::response;
    using tangents  = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::tangents;
    using increment = typename MarmotMaterialGeneralGradientEnhancedHypoElastic< nNonlocalVariables >::increment;

    for ( size_t i = 0; i < this->qps.size(); i++ ) {

      QuadraturePoint& qp = this->qps[i];

      const BSized&      B      = qp.B;
      const NSizedK&     N_K    = qp.N_K;
      const dNdXiSizedK& dNdX_K = qp.dNdX_K;

      dE = B * dQU;

      // delta of _K at Gausspoint
      Vector< double, nNonlocalVariables > K;
      Vector< double, nNonlocalVariables > dK;

      for ( size_t n = 0; n < nNonlocalVariables; n++ ) {
        double idx = n * nNonLocalNodes;
        K( n )     = N_K * qK.segment( idx, nNonLocalNodes );
        dK( n )    = N_K * dQK.segment( idx, nNonLocalNodes );
      }

      response  res;
      tangents  tan;
      increment inc;
      try {
        if constexpr ( nDim == 2 ) {
          Vector6d dE6  = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
          res.stress    = qp.managedStateVars->stress;
          res.stateVars = qp.managedStateVars->materialStateVars.data();
          inc           = { dE6, K, dK, time[1], dT };
          Matrix3d C;
          Vector3d S;

          if ( sectionType == SectionType::PlaneStress ) {
            qp.material->computePlaneStress( res, tan, inc );
            S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            C = ContinuumMechanics::PlaneStress::getPlaneStressTangent( tan.dStressddStrain );
          }
          else if ( sectionType == SectionType::PlaneStrain ) {
            qp.material->computeStress( res, tan, inc );
            S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            C = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( tan.dStressddStrain );
          }

          fU -= B.transpose() * S * qp.J0xW;
          kUU += B.transpose() * C * B * qp.J0xW;

          for ( size_t n = 0; n < nNonlocalVariables; n++ ) {
            double idx = n * nNonLocalNodes;
            fK.segment( idx, nNonLocalNodes ) -= ( N_K.transpose() * K( n ) +
                                                   res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                     qK.segment( idx, nNonLocalNodes ) -
                                                   N_K.transpose() * res.KLocal( n ) ) *
                                                 qp.J0xW;
            const auto dSdK         = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddK.col( n ) );
            const auto dK_Local_dDE = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt(
              tan.dKLocalddStrain.row( n ) );

            kUK.block( 0, idx, sizeDoFU, nNonLocalNodes ) += B.transpose() * dSdK * N_K * qp.J0xW;
            kKU.block( idx, 0, nNonLocalNodes, sizeDoFU ) += N_K.transpose() * -dK_Local_dDE.transpose() * B * qp.J0xW;
            kKK.block( idx, idx, nNonLocalNodes, nNonLocalNodes ) += ( N_K.transpose() * N_K +
                                                                       res.c( n ) * dNdX_K.transpose() * dNdX_K +
                                                                       tan.dcddK( n ) * dNdX_K.transpose() * dNdX_K *
                                                                         qK.segment( idx, nNonLocalNodes ) * N_K -
                                                                       N_K.transpose() * tan.dKLocalddK( n, n ) *
                                                                         N_K ) *
                                                                     qp.J0xW;
          }
        }

        else if ( nDim == 3 ) {
          if ( sectionType == Solid ) {

            S             = qp.managedStateVars->stress;
            res.stress    = S;
            res.stateVars = qp.managedStateVars->materialStateVars.data();
            inc           = { dE, K, dK, time[1], dT };
            qp.material->computeStress( res, tan, inc );

            fU -= B.transpose() * res.stress * qp.J0xW;
            kUU += B.transpose() * tan.dStressddStrain * B * qp.J0xW;

            for ( size_t n = 0; n < nNonlocalVariables; n++ ) {
              double idx = n * nNonLocalNodes;
              fK.segment( idx, nNonLocalNodes ) -= ( N_K.transpose() * K( n ) +
                                                     res.c( n ) * dNdX_K.transpose() * dNdX_K *
                                                       qK.segment( idx, nNonLocalNodes ) -
                                                     N_K.transpose() * res.KLocal( n ) ) *
                                                   qp.J0xW;

              kUK.block( 0, idx, sizeDoFU, nNonLocalNodes ) += B.transpose() * tan.dStressddK.col( n ) * N_K * qp.J0xW;
              kKU.block( idx, 0, nNonLocalNodes, sizeDoFU ) += N_K.transpose() * -tan.dKLocalddStrain.row( n ) * B *
                                                               qp.J0xW;

              kKK.block( idx, idx, nNonLocalNodes, nNonLocalNodes ) += ( N_K.transpose() * N_K +
                                                                         res.c( n ) * dNdX_K.transpose() * dNdX_K +
                                                                         tan.dcddK( n ) * dNdX_K.transpose() * dNdX_K *
                                                                           qK.segment( idx, nNonLocalNodes ) * N_K -
                                                                         N_K.transpose() * tan.dKLocalddK( n, n ) *
                                                                           N_K ) *
                                                                       qp.J0xW;
            }
          }
        }
      }
      catch ( std::exception& e ) {
        pNewDT = 0.25;
        return;
      }
      qp.managedStateVars->stress = res.stress;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                            double*                             P,
                            double*                             K,
                            const int                           elementFace,
                            const double*                       load,
                            const double*                       QTotal,
                            const double*                       time,
                            double                              dT )
  {

    Map< RhsSized >     P_( P );
    Ref< USizedVector > fU( P_.head( sizeDoFU ) );

    switch ( loadType ) {

    case MarmotElement::Pressure: {
      const double p = load[0];

      FiniteElement::BoundaryElement boundaryEl( localGeometryElement.shape,
                                                 elementFace,
                                                 nDim,
                                                 localGeometryElement.coordinates );

      VectorXd Pk = -p * boundaryEl.computeSurfaceNormalVectorialLoadVector();

      if ( nDim == 2 )
        Pk *= elementProperties[0]; // thickness

      boundaryEl.assembleIntoParentVectorial( Pk, fU );

      break;
    }
    default: {
      throw std::invalid_argument( "Invalid Load Type specified" );
    }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    setInitialConditions( StateTypes state, const double* values )
  {
    if constexpr ( nDim > 1 ) {
      switch ( state ) {
      case MarmotElement::GeostaticStress: {
        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = localGeometryElement.NB( localGeometryElement.N( qp.xi ) ) *
                                 localGeometryElement.coordinates;

          const double sigY1 = values[0];
          const double sigY2 = values[2];
          const double y1    = values[1];
          const double y2    = values[3];

          using namespace Math;
          qp.managedStateVars->stress( 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 ); // sigma_y
          qp.managedStateVars->stress( 0 ) = values[4] * qp.managedStateVars->stress( 1 );                 // sigma_x
          qp.managedStateVars->stress( 2 ) = values[5] * qp.managedStateVars->stress( 1 );                 // sigma_z
        }
        break;
      }
      case MarmotElement::MarmotMaterialStateVars: {
        throw std::invalid_argument( "Please use initializeStateVars directly on material" );
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  void GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    computeBodyForce( double*       P_,
                      double*       K,
                      const double* load,

                      const double* QTotal,
                      const double* time,
                      double        dT )
  {
    Map< RhsSized >                              P( P_ );
    Ref< USizedVector >                          Pe( P.head( sizeDoFU ) );
    const Map< const Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps )
      Pe += localGeometryElement.NB( localGeometryElement.N( qp.xi ) ).transpose() * f * qp.J0xW;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< double > GeneralGradientEnhancedDisplacementFiniteElement< nDim,
                                                                          nNodes,
                                                                          nNonlocalVariables,
                                                                          nNonLocalNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap = localGeometryElement.NB( localGeometryElement.N( centerXi ) ) * localGeometryElement.coordinates;
    return coords;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  std::vector< std::vector< double > > GeneralGradientEnhancedDisplacementFiniteElement<
    nDim,
    nNodes,
    nNonlocalVariables,
    nNonLocalNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = localGeometryElement.NB( localGeometryElement.N( qp.xi ) ) * localGeometryElement.coordinates;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes, int nNonlocalVariables, int nNonLocalNodes >
  int GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, nNonlocalVariables, nNonLocalNodes >::
    getNumberOfQuadraturePoints()
  {
    return qps.size();
  }
} // namespace Marmot::Elements
