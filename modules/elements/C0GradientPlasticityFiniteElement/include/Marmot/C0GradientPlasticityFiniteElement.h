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
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialGradientPlasticityHypoElastic.h"
#include "Marmot/MarmotMaterialGradientPlasticityHypoElasticFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;
using namespace FastorStandardTensors;

namespace Marmot::Elements {

  template < int nDim, int nNodesU, int nNodesL = nNodesU, int nNodesG = nNodesL >
  class C0GradientPlasticityFiniteElement : public MarmotElement {

  public:
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim;
    static constexpr int nDofPerNodeL = 1;
    static constexpr int nDofPerNodeG = nDim;

    static constexpr int nNodes = nNodesU >= nNodesL ? ( nNodesU >= nNodesG ? nNodesU : nNodesG )
                                                     : ( nNodesL >= nNodesG ? nNodesL : nNodesG );

    static constexpr int sizeDoFU       = nNodesU * nDofPerNodeU;
    static constexpr int sizeDoFL       = nNodesL * nDofPerNodeL;
    static constexpr int sizeDoFG       = nNodesG * nDofPerNodeG;
    static constexpr int sizeLoadVector = sizeDoFU + sizeDoFL + sizeDoFG;

    using ParentGeometryElementU = MarmotGeometryElement< nDim, nNodesU >;
    using ParentGeometryElementL = MarmotGeometryElement< nDim, nNodesL >;
    using ParentGeometryElementG = MarmotGeometryElement< nDim, nNodesG >;
    using ParentGeometryElement  = MarmotGeometryElement< nDim, nNodes >;
    using JacobianSizedU         = typename ParentGeometryElementU::JacobianSized;
    using JacobianSizedL         = typename ParentGeometryElementL::JacobianSized;
    using JacobianSizedG         = typename ParentGeometryElementG::JacobianSized;
    using NSizedU                = typename ParentGeometryElementU::NSized;
    using NSizedL                = typename ParentGeometryElementL::NSized;
    using NSizedG                = typename ParentGeometryElementG::NSized;
    using dNdXiSizedU            = typename ParentGeometryElementU::dNdXiSized;
    using dNdXiSizedL            = typename ParentGeometryElementL::dNdXiSized;
    using dNdXiSizedG            = typename ParentGeometryElementG::dNdXiSized;
    using BSizedU                = typename ParentGeometryElementU::BSized;
    using XiSized                = typename ParentGeometryElementU::XiSized;
    using RhsSized               = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix          = Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector           = Matrix< double, sizeDoFU, 1 >;
    using LSizedVector           = Matrix< double, sizeDoFL, 1 >;
    using GSizedVector           = Matrix< double, sizeDoFG, 1 >;
    using CSized = Matrix< double, ParentGeometryElementU::voigtSize, ParentGeometryElementU::voigtSize >;
    using Voigt  = Matrix< double, ParentGeometryElementU::voigtSize, 1 >;

    ParentGeometryElementU localGeometryElementU;
    ParentGeometryElementL localGeometryElementL;
    ParentGeometryElementG localGeometryElementG;
    ParentGeometryElement  localGeometryElement;

    Map< const VectorXd > elementProperties;
    const int             elLabel;
    const SectionType     sectionType;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double      detJ;
      double      J0xW;
      NSizedU     N_U;
      dNdXiSizedU dNdX_U;
      BSizedU     B_U;
      NSizedL     N_L;
      dNdXiSizedL dNdX_L;
      NSizedG     N_G;
      dNdXiSizedG dNdX_G;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 9 },
          { .name = "strain", .length = 9 },
          { .name = "total strain energy", .length = 1 },
          { .name = "elastic strain energy", .length = 1 },
          { .name = "dissipation", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        TensorMap33d                  stress;
        TensorMap33d                  strain;
        double&                       totalStrainEnergy;
        double&                       elasticStrainEnergy;
        double&                       dissipation;
        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; };

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            stress( &find( "stress" ) ),
            strain( &find( "strain" ) ),
            totalStrainEnergy( find( "total strain energy" ) ),
            elasticStrainEnergy( find( "elastic strain energy" ) ),
            dissipation( find( "dissipation" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      std::unique_ptr< MarmotMaterialGradientPlasticityHypoElastic< 1 > > material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      }

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          J0xW( 0.0 ),
          N_U( NSizedU::Zero() ),
          dNdX_U( dNdXiSizedU::Zero() ),
          B_U( BSizedU::Zero() ),
          N_L( NSizedL::Zero() ),
          dNdX_L( dNdXiSizedL::Zero() ),
          N_G( NSizedG::Zero() ),
          dNdX_G( dNdXiSizedG::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    C0GradientPlasticityFiniteElement( int                                         elementID,
                                       FiniteElement::Quadrature::IntegrationTypes integrationType,
                                       SectionType                                 sectionType )
      : elementProperties( nullptr, 0 ), elLabel( elementID ), sectionType( sectionType )
    {
      for ( const auto& qpInfo :
            FiniteElement::Quadrature::getGaussPointInfo( localGeometryElement.shape, integrationType ) ) {
        qps.emplace_back( qpInfo.xi, qpInfo.weight );
      }
    }

    int getNumberOfRequiredStateVars() { return qps[0].getNumberOfRequiredStateVars() * qps.size(); }

    std::vector< std::vector< std::string > > getNodeFields()
    {
      static std::vector< std::vector< std::string > > nodeFields;
      if ( nodeFields.empty() ) {
        nodeFields.resize( nNodes );
        for ( int i = 0; i < nNodes; i++ ) {
          if ( i < nNodesU )
            nodeFields[i].push_back( "displacement" );
          if ( i < nNodesL )
            nodeFields[i].push_back( "plastic multiplier" );
          if ( i < nNodesG )
            nodeFields[i].push_back( "plastic multiplier gradient" );
        }
      }
      return nodeFields;
    }

    std::vector< int > getDofIndicesPermutationPattern()
    {
      static std::vector< int > permutationPattern;
      if ( permutationPattern.empty() ) {
        int totalNodes = std::max( { nNodesU, nNodesL, nNodesG } );

        // 1. Calculate how many elements are in each global block
        // so we can figure out where the 'l' and 'g' global layouts start.
        int globalLStart = nNodesU * nDim;
        int globalGStart = globalLStart + nNodesL; // (1 DoF per node for L)

        // Tracks the next available global index for each DoF type
        int globalUIdx = 0;
        int globalLIdx = globalLStart;
        int globalGIdx = globalGStart;

        // 2. Pre-calculate the local memory offset for every node.
        // Each node only takes up space for the DoFs it actually contains.
        std::vector< int > localNodeOffsets( totalNodes, 0 );
        int                currentLocalOffset = 0;
        for ( int i = 0; i < totalNodes; ++i ) {
          localNodeOffsets[i] = currentLocalOffset;

          if ( i < nNodesU )
            currentLocalOffset += nDim;
          if ( i < nNodesL )
            currentLocalOffset += 1;
          if ( i < nNodesG )
            currentLocalOffset += nDim;
        }

        // 3. Size the permutation pattern array to hold all global DoFs
        int totalGlobalDofs = ( nNodesU * nDim ) + nNodesL + ( nNodesG * nDim );
        permutationPattern.resize( totalGlobalDofs );

        // 4. Map each node's local components to their true global positions
        for ( int i = 0; i < totalNodes; ++i ) {
          int nodeLocalStart  = localNodeOffsets[i];
          int currentLocalDof = 0;

          // Map U DoFs if this node has them
          if ( i < nNodesU ) {
            for ( int j = 0; j < nDim; ++j ) {
              permutationPattern[globalUIdx++] = nodeLocalStart + currentLocalDof++;
            }
          }

          // Map L DoF if this node has it
          if ( i < nNodesL ) {
            permutationPattern[globalLIdx++] = nodeLocalStart + currentLocalDof++;
          }

          // Map G DoFs if this node has them
          if ( i < nNodesG ) {
            for ( int j = 0; j < nDim; ++j ) {
              permutationPattern[globalGIdx++] = nodeLocalStart + currentLocalDof++;
            }
          }
        }
      }
      return permutationPattern;
    }

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return localGeometryElement.getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars )
    {
      const int nQpStateVars = nStateVars / qps.size();

      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto&   qp          = qps[i];
        double* qpStateVars = stateVars + ( i * nQpStateVars );
        qp.assignStateVars( qpStateVars, nQpStateVars );
      }
    }

    void assignProperty( const ElementProperties& elementPropertiesInfo )
    {
      new ( &elementProperties )
        Map< const VectorXd >( elementPropertiesInfo.elementProperties, elementPropertiesInfo.nElementProperties );
    }

    void assignProperty( const MarmotMaterialSection& section )
    {
      for ( auto& qp : qps ) {
        qp.material = std::unique_ptr< MarmotMaterialGradientPlasticityHypoElastic< 1 > >(
          MarmotLibrary::MarmotMaterialGradientPlasticityHypoElasticFactory< 1 >::
            createMaterial( section.materialName, section.materialProperties, section.nMaterialProperties, elLabel ) );
      }
    }

    void assignNodeCoordinates( const double* coordinates )
    {
      localGeometryElement.assignNodeCoordinates( coordinates );
      localGeometryElementU.assignNodeCoordinates( coordinates );
      localGeometryElementL.assignNodeCoordinates( coordinates );
      localGeometryElementG.assignNodeCoordinates( coordinates );
    }

    void initializeYourself()
    {
      for ( QuadraturePoint& qp : qps ) {
        const auto           dNdXi_U = localGeometryElementU.dNdXi( qp.xi );
        const JacobianSizedU J_U     = localGeometryElementU.Jacobian( dNdXi_U );
        const JacobianSizedU JInv_U  = J_U.inverse();
        const auto           dNdXi_L = localGeometryElementL.dNdXi( qp.xi );
        const JacobianSizedL J_L     = localGeometryElementL.Jacobian( dNdXi_L );
        const JacobianSizedL JInv_L  = J_L.inverse();
        const auto           dNdXi_G = localGeometryElementG.dNdXi( qp.xi );
        const JacobianSizedG J_G     = localGeometryElementG.Jacobian( dNdXi_G );
        const JacobianSizedG JInv_G  = J_G.inverse();
        qp.detJ                      = J_U.determinant();
        qp.N_U                       = localGeometryElementU.N( qp.xi );
        qp.dNdX_U                    = localGeometryElementU.dNdX( dNdXi_U, JInv_U );
        qp.B_U                       = localGeometryElementU.B( qp.dNdX_U );
        qp.N_L                       = localGeometryElementL.N( qp.xi );
        qp.dNdX_L                    = localGeometryElementL.dNdX( dNdXi_L, JInv_L );
        qp.N_G                       = localGeometryElementG.N( qp.xi );
        qp.dNdX_G                    = localGeometryElementG.dNdX( dNdXi_G, JInv_G );

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

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
    {
      using namespace Marmot::FastorIndices;

      using Bijk    = Fastor::Index< B_, i_, j_, k_ >;
      using to_iAkB = Fastor::OIndex< i_, A_, k_, B_ >;
      using to_iAB  = Fastor::OIndex< i_, A_, B_ >;
      using to_AiB  = Fastor::OIndex< A_, i_, B_ >;
      using to_iAjB = Fastor::OIndex< i_, A_, j_, B_ >;

      Map< const RhsSized > QTotal( QTotal_ );
      Map< const RhsSized > dQ( dQ_ );
      Map< KeSizedMatrix >  Ke( Ke_ );
      Map< RhsSized >       Pe( Pe_ );

      const auto qU_np = Fastor::TensorMap< const double, nNodesU, nDim >( dQ.data() );
      const auto qL_np = Fastor::TensorMap< const double, nNodesL >( dQ.data() + sizeDoFU );
      const auto qG_np = Fastor::TensorMap< const double, nNodesG, nDim >( dQ.data() + sizeDoFU + sizeDoFL );

      const auto QG_np = Fastor::TensorMap< const double, nNodesG, nDim >( QTotal.data() + sizeDoFU + sizeDoFL );
      const auto QL_np = Fastor::TensorMap< const double, nNodesL >( QTotal.data() + sizeDoFU );

      Ref< Matrix< double, sizeDoFU, sizeDoFU > > KUU( Ke.topLeftCorner( sizeDoFU, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFU, sizeDoFL > > KUL( Ke.block( 0, sizeDoFU, sizeDoFU, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFU, sizeDoFG > > KUG( Ke.topRightCorner( sizeDoFU, sizeDoFG ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFU > > KLU( Ke.block( sizeDoFU, 0, sizeDoFL, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFL > > KLL( Ke.block( sizeDoFU, sizeDoFU, sizeDoFL, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFG > > KLG( Ke.block( sizeDoFU, sizeDoFU + sizeDoFL, sizeDoFL, sizeDoFG ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFU > > KGU( Ke.block( sizeDoFU + sizeDoFL, 0, sizeDoFG, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFL > > KGL( Ke.block( sizeDoFU + sizeDoFL, sizeDoFU, sizeDoFG, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFG > > KGG( Ke.bottomRightCorner( sizeDoFG, sizeDoFG ) );

      Fastor::TensorMap< double, nNodesU, nDim > r_U( Pe.data() );
      Fastor::TensorMap< double, nNodesL >       r_L( Pe.data() + sizeDoFU );
      Fastor::TensorMap< double, nNodesG, nDim > r_G( Pe.data() + sizeDoFU + sizeDoFL );

      Fastor::Tensor< double, nDim, nNodesU, nDim, nNodesU > k_UU( 0.0 );
      Fastor::Tensor< double, nDim, nNodesU, nNodesL >       k_UL( 0.0 );
      Fastor::Tensor< double, nDim, nNodesU, nDim, nNodesG > k_UG( 0.0 );
      Fastor::Tensor< double, nNodesL, nDim, nNodesU >       k_LU( 0.0 );
      Fastor::Tensor< double, nNodesL, nNodesL >             k_LL( 0.0 );
      Fastor::Tensor< double, nNodesL, nDim, nNodesG >       k_LG( 0.0 );
      Fastor::Tensor< double, nDim, nNodesG, nDim, nNodesU > k_GU( 0.0 );
      Fastor::Tensor< double, nDim, nNodesG, nNodesL >       k_GL( 0.0 );
      Fastor::Tensor< double, nDim, nNodesG, nDim, nNodesG > k_GG( 0.0 );

      using response  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::response;
      using tangents  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::tangents;
      using increment = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::increment;

      const double penalty = 1e11;

      for ( auto& qp : qps ) {
        const auto dNdX_U = Fastor::Tensor< double, nDim, nNodesU >( qp.dNdX_U.data(), Fastor::ColumnMajor );
        const auto N_L    = Fastor::Tensor< double, nNodesL >( qp.N_L.data() );
        const auto dNdX_L = Fastor::Tensor< double, nDim, nNodesL >( qp.dNdX_L.data(), Fastor::ColumnMajor );
        const auto N_G    = Fastor::Tensor< double, nNodesG >( qp.N_G.data() );
        const auto dNdX_G = Fastor::Tensor< double, nDim, nNodesG >( qp.dNdX_G.data(), Fastor::ColumnMajor );

        const auto gradU = Fastor::evaluate( Fastor::einsum< Ai, jA >( qU_np, dNdX_U ) );

        Tensor33d dStrain = Tensor33d( 0.0 );
        for ( int i = 0; i < nDim; i++ ) {
          for ( int j = 0; j < nDim; j++ ) {
            dStrain( i, j ) = 0.5 * ( gradU( i, j ) + gradU( j, i ) );
          }
        }

        Fastor::Tensor< double, nDim, nDim, nNodesU, nDim > dStrain_dU;
        for ( int i = 0; i < nDim; i++ ) {
          for ( int j = 0; j < nDim; j++ ) {
            for ( int n = 0; n < nNodesU; n++ ) {
              for ( int d = 0; d < nDim; d++ ) {
                dStrain_dU( i, j, n, d ) = 0.5 * ( dNdX_U( d, n ) * Spatial3D::I( i, j ) +
                                                   dNdX_U( i, n ) * Spatial3D::I( j, d ) );
              }
            }
          }
        }

        const double lambdaIncrement        = Fastor::evaluate( Fastor::einsum< A, A >( N_L, qL_np ) ).toscalar();
        const double laplaceLambdaIncrement = Fastor::evaluate( Fastor::einsum< iA, Ai >( dNdX_G, qG_np ) ).toscalar();

        auto N_times_qG    = Fastor::evaluate( Fastor::einsum< A, Ai >( N_G, qG_np ) );
        auto dNdX_times_qL = Fastor::evaluate( Fastor::einsum< iA, A >( dNdX_L, qL_np ) );

        Fastor::Tensor< double, nDim > cConstraint    = N_times_qG - dNdX_times_qL;
        Fastor::Tensor< double, nDim > augmentedForce = penalty * cConstraint;

        response  res;
        tangents  tan;
        increment inc;

        res.stress               = qp.managedStateVars->stress;
        res.f                    = 0.0;
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();

        inc.dLambda( 0 )        = lambdaIncrement;
        inc.laplaceDLambda( 0 ) = laplaceLambdaIncrement;
        inc.time                = time;
        inc.dT                  = dT;
        inc.dStrain             = dStrain;

        if constexpr ( nDim == 2 ) {
          if ( sectionType == SectionType::PlaneStress ) {
            qp.material->computePlaneStress( res, tan, inc );
          }
          else if ( sectionType == SectionType::PlaneStrain ) {
            qp.material->computeStress( res, tan, inc );
          }
          else {
            throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
          }
        }
        else if constexpr ( nDim == 3 ) {
          if ( sectionType != SectionType::Solid )
            throw std::invalid_argument( "Invalid section type for 3D element, expected Solid" );
          qp.material->computeStress( res, tan, inc );
        }

        Fastor::Tensor< double, nDim, nDim > stress_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            stress_nDim( i, j ) = res.stress( i, j );

        Fastor::Tensor< double, nDim, nDim, nDim, nDim > C_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            for ( int k = 0; k < nDim; ++k )
              for ( int l = 0; l < nDim; ++l )
                C_nDim( i, j, k, l ) = tan.dStressddStrain( i, j, k, l );

        Fastor::Tensor< double, nDim, nDim > dSdLambda_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            dSdLambda_nDim( i, j ) = tan.dStressddLambda( i, j, 0 );

        Fastor::Tensor< double, nDim, nDim > dSdLaplacian_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            dSdLaplacian_nDim( i, j ) = tan.dStressddLaplacian( i, j, 0 );

        Fastor::Tensor< double, nDim, nDim > dFddE_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            dFddE_nDim( i, j ) = tan.dFddStrain( 0, i, j );

        const double f          = res.f( 0 );
        const double dFddLambda = tan.dFddLambda( 0, 0 );
        const double dFddLap    = tan.dFddLaplacian( 0, 0 );

        r_U += Fastor::evaluate( Fastor::einsum< iA, ij >( dNdX_U, stress_nDim ) ) * qp.J0xW;
        r_L += ( N_L * f - Fastor::evaluate( Fastor::einsum< iA, i >( dNdX_L, augmentedForce ) ) ) * qp.J0xW;
        r_G += Fastor::evaluate( Fastor::einsum< A, i >( N_G, augmentedForce ) ) * qp.J0xW;

        auto dStress_dU = Fastor::evaluate( Fastor::einsum< ijkl, klmn >( C_nDim, dStrain_dU ) );
        k_UU += Fastor::evaluate( Fastor::einsum< iA, ijkB, to_jAkB >( dNdX_U, dStress_dU ) ) * qp.J0xW;

        k_UL += Fastor::evaluate( Fastor::einsum< jA, ij, B, to_iAB >( dNdX_U, dSdLambda_nDim, N_L ) ) * qp.J0xW;

        auto divOp_G_trans = Fastor::evaluate( Fastor::einsum< jB, ij >( dNdX_G, dSdLaplacian_nDim ) );
        k_UG += Fastor::evaluate( Fastor::einsum< jA, iB, to_iAkB >( dNdX_U, divOp_G_trans ) ) * qp.J0xW;

        auto dF_dU = Fastor::evaluate( Fastor::einsum< ij, jB >( dFddE_nDim, dNdX_U ) );
        k_LU += Fastor::evaluate( Fastor::einsum< A, iB, to_AiB >( N_L, dF_dU ) ) * qp.J0xW;

        k_LL += ( Fastor::evaluate( Fastor::einsum< A, B >( N_L, N_L ) ) * dFddLambda +
                  penalty * Fastor::evaluate( Fastor::einsum< iA, iB >( dNdX_L, dNdX_L ) ) ) *
                qp.J0xW;

        k_LG += ( Fastor::evaluate( Fastor::einsum< A, iB >( N_L, dNdX_G ) ) * dFddLap -
                  penalty * Fastor::evaluate( Fastor::einsum< iA, B, to_AiB >( dNdX_L, N_G ) ) ) *
                qp.J0xW;

        k_GU += Fastor::Tensor< double, nDim, nNodesG, nDim, nNodesU >( 0.0 );

        k_GL += ( -penalty * Fastor::evaluate( Fastor::einsum< A, iB, to_iAB >( N_G, dNdX_L ) ) ) * qp.J0xW;

        k_GG += ( penalty *
                  Fastor::evaluate( Fastor::einsum< A, B, ij, to_iAjB >( N_G,
                                                                         N_G,
                                                                         Fastor::Tensor< double, nDim, nDim >(
                                                                           Marmot::FastorStandardTensors::Spatial3D::
                                                                             I( Fastor::fseq< 0, nDim >(),
                                                                                Fastor::fseq< 0, nDim >() ) ) ) ) ) *
                qp.J0xW;

        qp.managedStateVars->stress              = res.stress;
        qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
        qp.managedStateVars->strain += dStrain;
      }

      KUU += Map< Matrix< double, sizeDoFU, sizeDoFU > >( Fastor::torowmajor( k_UU ).data() );
      KUL += Map< Matrix< double, sizeDoFU, sizeDoFL > >( Fastor::torowmajor( k_UL ).data() );
      KUG += Map< Matrix< double, sizeDoFU, sizeDoFG > >( Fastor::torowmajor( k_UG ).data() );
      KLU += Map< Matrix< double, sizeDoFL, sizeDoFU > >( Fastor::torowmajor( k_LU ).data() );
      KLL += Map< Matrix< double, sizeDoFL, sizeDoFL > >( Fastor::torowmajor( k_LL ).data() );
      KLG += Map< Matrix< double, sizeDoFL, sizeDoFG > >( Fastor::torowmajor( k_LG ).data() );
      KGU += Map< Matrix< double, sizeDoFG, sizeDoFU > >( Fastor::torowmajor( k_GU ).data() );
      KGL += Map< Matrix< double, sizeDoFG, sizeDoFL > >( Fastor::torowmajor( k_GL ).data() );
      KGG += Map< Matrix< double, sizeDoFG, sizeDoFG > >( Fastor::torowmajor( k_GG ).data() );
    }

    void computeKernelsExplicit( const double* QTotal_, const double* dQ_, double* Pe_, double time, double dT )
    {
      using namespace Marmot::FastorIndices;

      using Bijk    = Fastor::Index< B_, i_, j_, k_ >;
      using to_iAkB = Fastor::OIndex< i_, A_, k_, B_ >;
      using to_iAB  = Fastor::OIndex< i_, A_, B_ >;
      using to_AiB  = Fastor::OIndex< A_, i_, B_ >;
      using to_iAjB = Fastor::OIndex< i_, A_, j_, B_ >;

      Map< const RhsSized > QTotal( QTotal_ );
      Map< const RhsSized > dQ( dQ_ );
      Map< RhsSized >       Pe( Pe_ );

      const auto qU_np = Fastor::TensorMap< const double, nNodesU, nDim >( dQ.data() );
      const auto qL_np = Fastor::TensorMap< const double, nNodesL >( dQ.data() + sizeDoFU );
      const auto qG_np = Fastor::TensorMap< const double, nNodesG, nDim >( dQ.data() + sizeDoFU + sizeDoFL );

      const auto QG_np = Fastor::TensorMap< const double, nNodesG, nDim >( QTotal.data() + sizeDoFU + sizeDoFL );
      const auto QL_np = Fastor::TensorMap< const double, nNodesL >( QTotal.data() + sizeDoFU );

      Fastor::TensorMap< double, nNodesU, nDim > r_U( Pe.data() );
      Fastor::TensorMap< double, nNodesL >       r_L( Pe.data() + sizeDoFU );
      Fastor::TensorMap< double, nNodesG, nDim > r_G( Pe.data() + sizeDoFU + sizeDoFL );

      using response  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::response;
      using increment = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::increment;

      const double penalty = 1e11;

      for ( auto& qp : qps ) {
        const auto dNdX_U = Fastor::Tensor< double, nDim, nNodesU >( qp.dNdX_U.data(), Fastor::ColumnMajor );
        const auto N_L    = Fastor::Tensor< double, nNodesL >( qp.N_L.data() );
        const auto dNdX_L = Fastor::Tensor< double, nDim, nNodesL >( qp.dNdX_L.data(), Fastor::ColumnMajor );
        const auto N_G    = Fastor::Tensor< double, nNodesG >( qp.N_G.data() );
        const auto dNdX_G = Fastor::Tensor< double, nDim, nNodesG >( qp.dNdX_G.data(), Fastor::ColumnMajor );

        const auto gradU = Fastor::evaluate( Fastor::einsum< Ai, jA >( qU_np, dNdX_U ) );

        Tensor33d dStrain = Tensor33d( 0.0 );
        for ( int i = 0; i < nDim; i++ ) {
          for ( int j = 0; j < nDim; j++ ) {
            dStrain( i, j ) = 0.5 * ( gradU( i, j ) + gradU( j, i ) );
          }
        }

        const double lambdaIncrement        = Fastor::evaluate( Fastor::einsum< A, A >( N_L, qL_np ) ).toscalar();
        const double laplaceLambdaIncrement = Fastor::evaluate( Fastor::einsum< iA, Ai >( dNdX_G, qG_np ) ).toscalar();

        auto N_times_qG    = Fastor::evaluate( Fastor::einsum< A, Ai >( N_G, QG_np ) );
        auto dNdX_times_qL = Fastor::evaluate( Fastor::einsum< iA, A >( dNdX_L, QL_np ) );

        Fastor::Tensor< double, nDim > cConstraint    = N_times_qG - dNdX_times_qL;
        Fastor::Tensor< double, nDim > augmentedForce = penalty * cConstraint;

        response  res;
        increment inc;

        res.stress               = qp.managedStateVars->stress;
        res.f                    = 0.0;
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();

        inc.dLambda( 0 )        = lambdaIncrement;
        inc.laplaceDLambda( 0 ) = laplaceLambdaIncrement;
        inc.time                = time;
        inc.dT                  = dT;
        inc.dStrain             = dStrain;

        if constexpr ( nDim == 2 ) {
          if ( sectionType == SectionType::PlaneStress ) {
            qp.material->computePlaneStressExplicit( res, inc );
          }
          else if ( sectionType == SectionType::PlaneStrain ) {
            qp.material->computeStressExplicit( res, inc );
          }
        }
        else if constexpr ( nDim == 3 ) {
          qp.material->computeStressExplicit( res, inc );
        }

        Fastor::Tensor< double, nDim, nDim > stress_nDim;
        for ( int i = 0; i < nDim; ++i )
          for ( int j = 0; j < nDim; ++j )
            stress_nDim( i, j ) = res.stress( i, j );

        const double f = res.f( 0 );

        r_U -= Fastor::evaluate( Fastor::einsum< jA, ij >( dNdX_U, stress_nDim ) ) * qp.J0xW;
        r_L -= ( N_L * f - Fastor::evaluate( Fastor::einsum< iA, i >( dNdX_L, augmentedForce ) ) ) * qp.J0xW;
        r_G -= Fastor::evaluate( Fastor::einsum< A, i >( N_G, augmentedForce ) ) * qp.J0xW;

        qp.managedStateVars->stress              = res.stress;
        qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
        qp.managedStateVars->strain += dStrain;
      }
    }
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 int                                 elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
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
          Pk *= elementProperties[0];

        boundaryEl.assembleIntoParentVectorial( Pk, fU );

        break;
      }
      case MarmotElement::SurfaceTraction: {
        FiniteElement::BoundaryElement boundaryEl( localGeometryElement.shape,
                                                   elementFace,
                                                   nDim,
                                                   localGeometryElement.coordinates );

        const XiSized tractionVector( load );

        auto Pk = boundaryEl.computeVectorialLoadVector( tractionVector );
        if ( nDim == 2 )
          Pk *= elementProperties[0];
        boundaryEl.assembleIntoParentVectorial( Pk, fU );

        break;
      }
      default: throw std::invalid_argument( "Invalid load type specified" );
      }
    }

    void computeBodyForce( double* P_, double* K, const double* load, const double* QTotal, double time, double dT )
    {
      Map< RhsSized >                              P( P_ );
      Ref< USizedVector >                          Pe( P.head( sizeDoFU ) );
      const Map< const Matrix< double, nDim, 1 > > f( load );

      for ( const auto& qp : qps )
        Pe += localGeometryElementU.NB( qp.N_U ).transpose() * f * qp.J0xW;
    }

    void computeConsistentInertia( double* M )
    {
      Map< KeSizedMatrix > Me( M );
      Me.setZero();

      for ( const auto& qp : qps ) {
        const auto                  N_u  = localGeometryElementU.NB( qp.N_U );
        const auto                  N_g  = localGeometryElementG.NB( qp.N_G );
        const double                rho  = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
        const std::vector< double > etaV = qp.material->getNonlocalViscosity(
          qp.managedStateVars->materialStateVars.data() );
        const double eta = etaV.empty() ? 0.0 : etaV[0];

        Me.topLeftCorner( sizeDoFU, sizeDoFU ) += N_u.transpose() * N_u * qp.J0xW * rho;
        Me.block( sizeDoFU, sizeDoFU, sizeDoFL, sizeDoFL ) += qp.N_L.transpose() * qp.N_L * qp.J0xW * eta;
        Me.bottomRightCorner( sizeDoFG, sizeDoFG ) += N_g.transpose() * N_g * qp.J0xW * eta;
      }
    }

    void computeLumpedInertia( double* M )
    {
      Map< RhsSized > LMM( M );
      LMM.setZero();

      constexpr int nNodesLinear  = ( 1 << nDim );
      auto          linGeometryEl = MarmotGeometryElement< nDim, nNodesLinear >();

      auto computeWeightedN = [&]( const auto& NField, const XiSized& xi ) {
        VectorXd   N_weighted = 0.5 * NField;
        const auto N_lin      = linGeometryEl.N( xi );
        N_weighted.head( nNodesLinear ) += 0.5 * N_lin;
        return N_weighted;
      };

      for ( const auto& qp : qps ) {
        const VectorXd N_weightedU = computeWeightedN( qp.N_U, qp.xi );
        const VectorXd N_weightedL = computeWeightedN( qp.N_L, qp.xi );
        const VectorXd N_weightedG = computeWeightedN( qp.N_G, qp.xi );

        const double                rho  = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
        const std::vector< double > etaV = qp.material->getNonlocalViscosity(
          qp.managedStateVars->materialStateVars.data() );
        const double eta = etaV.empty() ? 0.0 : etaV[0];

        const VectorXd mU = N_weightedU * qp.J0xW * rho;
        const VectorXd mL = N_weightedL * qp.J0xW * eta;
        const VectorXd mG = N_weightedG * qp.J0xW * eta;

        for ( int i = 0; i < nNodesU; i++ ) {
          for ( int d = 0; d < nDim; d++ )
            LMM( i * nDim + d ) += mU( i );
        }

        for ( int i = 0; i < nNodesL; i++ ) {
          LMM( sizeDoFU + i ) += mL( i );
        }

        for ( int i = 0; i < nNodesG; i++ ) {
          for ( int d = 0; d < nDim; d++ )
            LMM( sizeDoFU + sizeDoFL + i * nDim + d ) += mG( i );
        }
      }
    }

    void computeCriticalTimeStepForExplicitDynamics( double& criticalTimeStep, const double* QTotal )
    {
      using response = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::response;

      criticalTimeStep = std::numeric_limits< double >::max();
      for ( const auto& qp : qps ) {
        double characteristicElementLength = 0.0;
        if constexpr ( nDim == 3 )
          characteristicElementLength = std::cbrt( 8 * qp.detJ );
        if constexpr ( nDim == 2 )
          characteristicElementLength = std::sqrt( 4 * qp.detJ );
        if constexpr ( nDim == 1 )
          characteristicElementLength = 2 * qp.detJ;

        response waveSpeedResponse;
        waveSpeedResponse.stress               = Tensor33d( qp.managedStateVars->stress.data() );
        waveSpeedResponse.f                    = 0;
        waveSpeedResponse.stateVars            = qp.managedStateVars->materialStateVars.data();
        waveSpeedResponse.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        waveSpeedResponse.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;

        const double c = qp.material->getMaximumWaveSpeed( waveSpeedResponse );
        if ( c <= 0.0 )
          throw std::runtime_error( "Material returned non-positive wave speed, cannot compute critical time step" );
        const double dt = characteristicElementLength / c;
        if ( dt < criticalTimeStep )
          criticalTimeStep = dt;
      }
    }

    void computeInternalEnergy( double& internalEnergy )
    {
      internalEnergy = 0.0;
      for ( const auto& qp : qps ) {
        internalEnergy += qp.managedStateVars->totalStrainEnergy;
      }
    }

    void setInitialConditions( StateTypes state, const double* values )
    {
      if constexpr ( nDim > 1 ) {
        switch ( state ) {
        case MarmotElement::GeostaticStress: {
          for ( QuadraturePoint& qp : qps ) {

            XiSized coordAtGauss = localGeometryElementU.NB( localGeometryElementU.N( qp.xi ) ) *
                                   localGeometryElementU.coordinates;

            const double sigY1 = values[0];
            const double sigY2 = values[2];
            const double y1    = values[1];
            const double y2    = values[3];

            using namespace Math;
            qp.managedStateVars->stress( 1, 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 );
            qp.managedStateVars->stress( 0, 0 ) = values[4] * qp.managedStateVars->stress( 1, 1 );
            qp.managedStateVars->stress( 2, 2 ) = values[5] * qp.managedStateVars->stress( 1, 1 );
          }
          break;
        }
        case MarmotElement::MarmotMaterialStateVars: {
          throw std::invalid_argument( "Please use initializeStateVars directly on material" );
        }
        case MarmotElement::MarmotMaterialInitialization: {
          for ( QuadraturePoint& qp : qps ) {
            qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                             qp.managedStateVars->materialStateVars.size() );
          }
          break;
        }
        default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
        }
      }
    }

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( stateName == "sdv" ) {
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

    std::vector< double > getCoordinatesAtCenter()
    {
      std::vector< double > coords( nDim );

      Eigen::Map< XiSized > coordsMap( &coords[0] );
      const auto            centerXi = XiSized::Zero();
      coordsMap = localGeometryElementU.NB( localGeometryElementU.N( centerXi ) ) * localGeometryElementU.coordinates;
      return coords;
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints()
    {
      std::vector< std::vector< double > > listedCoords;

      std::vector< double > coords( nDim );
      Eigen::Map< XiSized > coordsMap( &coords[0] );

      for ( const auto& qp : qps ) {
        coordsMap = localGeometryElementU.NB( localGeometryElementU.N( qp.xi ) ) * localGeometryElementU.coordinates;
        listedCoords.push_back( coords );
      }

      return listedCoords;
    }

    int getNumberOfQuadraturePoints() { return qps.size(); }
  };

} // namespace Marmot::Elements