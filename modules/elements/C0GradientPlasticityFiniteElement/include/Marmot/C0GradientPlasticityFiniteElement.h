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

namespace Marmot::Elements {

  template < int nDim, int nNodes >
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

    static constexpr int sizeDoFU       = nNodes * nDofPerNodeU;
    static constexpr int sizeDoFL       = nNodes * nDofPerNodeL;
    static constexpr int sizeDoFG       = nNodes * nDofPerNodeG;
    static constexpr int sizeLoadVector = sizeDoFU + sizeDoFL + sizeDoFG;

    using ParentGeometryElement = MarmotGeometryElement< nDim, nNodes >;
    using JacobianSized         = typename ParentGeometryElement::JacobianSized;
    using NSized                = typename ParentGeometryElement::NSized;
    using dNdXiSized            = typename ParentGeometryElement::dNdXiSized;
    using BSized                = typename ParentGeometryElement::BSized;
    using XiSized               = typename ParentGeometryElement::XiSized;
    using RhsSized              = Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix         = Matrix< double, sizeLoadVector, sizeLoadVector >;
    using USizedVector          = Matrix< double, sizeDoFU, 1 >;
    using LSizedVector          = Matrix< double, sizeDoFL, 1 >;
    using GSizedVector          = Matrix< double, sizeDoFG, 1 >;
    using CSized                = Matrix< double, ParentGeometryElement::voigtSize, ParentGeometryElement::voigtSize >;
    using Voigt                 = Matrix< double, ParentGeometryElement::voigtSize, 1 >;

    ParentGeometryElement localGeometryElement;

    Map< const VectorXd > elementProperties;
    const int             elLabel;
    const SectionType     sectionType;

    struct QuadraturePoint {

      const XiSized xi;
      const double  weight;

      double     detJ;
      double     J0xW;
      NSized     N;
      dNdXiSized dNdX;
      BSized     B;

      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "stress", .length = 6 },
          { .name = "strain", .length = 6 },
          { .name = "total strain energy", .length = 1 },
          { .name = "elastic strain energy", .length = 1 },
          { .name = "dissipation", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        mVector6d                     stress;
        mVector6d                     strain;
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
          N( NSized::Zero() ),
          dNdX( dNdXiSized::Zero() ),
          B( BSized::Zero() )
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
          nodeFields[i].push_back( "displacement" );
          nodeFields[i].push_back( "plastic multiplier" );
          nodeFields[i].push_back( "plastic multiplier gradient" );
        }
      }

      return nodeFields;
    }

    std::vector< int > getDofIndicesPermutationPattern()
    {
      static std::vector< int > permutationPattern;

      if ( permutationPattern.empty() ) {
        for ( int i = 0; i < nNodes; i++ ) {
          for ( int j = 0; j < nDim; j++ ) {
            permutationPattern.push_back( i * nDim + j );
          }

          permutationPattern.push_back( sizeDoFU + i );

          for ( int j = 0; j < nDim; j++ ) {
            permutationPattern.push_back( sizeDoFU + sizeDoFL + i * nDim + j );
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
    }

    void initializeYourself()
    {
      for ( QuadraturePoint& qp : qps ) {
        const auto          dNdXi = localGeometryElement.dNdXi( qp.xi );
        const JacobianSized J     = localGeometryElement.Jacobian( dNdXi );
        const JacobianSized JInv  = J.inverse();
        qp.detJ                   = J.determinant();
        qp.N                      = localGeometryElement.N( qp.xi );
        qp.dNdX                   = localGeometryElement.dNdX( dNdXi, JInv );
        qp.B                      = localGeometryElement.B( qp.dNdX );

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

    void computeYourself( const double* QTotal_,
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

      const Ref< const USizedVector > qU( QTotal.head( sizeDoFU ) );
      const Ref< const LSizedVector > qL( QTotal.segment( sizeDoFU, sizeDoFL ) );
      const Ref< const GSizedVector > qG( QTotal.tail( sizeDoFG ) );

      const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
      const Ref< const LSizedVector > dQL( dQ.segment( sizeDoFU, sizeDoFL ) );
      const Ref< const GSizedVector > dQG( dQ.tail( sizeDoFG ) );

      Ref< Matrix< double, sizeDoFU, sizeDoFU > > KUU( Ke.topLeftCorner( sizeDoFU, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFU, sizeDoFL > > KUL( Ke.block( 0, sizeDoFU, sizeDoFU, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFU, sizeDoFG > > KUG( Ke.topRightCorner( sizeDoFU, sizeDoFG ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFU > > KLU( Ke.block( sizeDoFU, 0, sizeDoFL, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFL > > KLL( Ke.block( sizeDoFU, sizeDoFU, sizeDoFL, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFL, sizeDoFG > > KLG( Ke.block( sizeDoFU, sizeDoFU + sizeDoFL, sizeDoFL, sizeDoFG ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFU > > KGU( Ke.block( sizeDoFU + sizeDoFL, 0, sizeDoFG, sizeDoFU ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFL > > KGL( Ke.block( sizeDoFU + sizeDoFL, sizeDoFU, sizeDoFG, sizeDoFL ) );
      Ref< Matrix< double, sizeDoFG, sizeDoFG > > KGG( Ke.bottomRightCorner( sizeDoFG, sizeDoFG ) );

      Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
      Ref< LSizedVector > fL( Pe.segment( sizeDoFU, sizeDoFL ) );
      Ref< GSizedVector > fG( Pe.tail( sizeDoFG ) );

      using response  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::response;
      using tangents  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::tangents;
      using increment = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::increment;

      const double penalty = 1e14;

      for ( auto& qp : qps ) {
        const BSized&    B    = qp.B;
        const NSized&    N    = qp.N;
        const dNdXiSized dNdX = qp.dNdX;
        const auto       NVec = localGeometryElement.NB( N );

        Matrix< double, 1, sizeDoFG > divOperator = Matrix< double, 1, sizeDoFG >::Zero();
        for ( int a = 0; a < nNodes; a++ ) {
          for ( int i = 0; i < nDim; i++ ) {
            divOperator( 0, a * nDim + i ) = dNdX( i, a );
          }
        }

        Voigt dE = B * dQU;

        const double lambdaIncrement        = N * dQL;
        const double laplaceLambdaIncrement = divOperator * dQG;

        const Matrix< double, nDim, 1 > gradLambdaInc = NVec * dQG;
        const Matrix< double, nDim, 1 > cConstraint   = NVec * dQG - dNdX * dQL;

        response  res;
        tangents  tan;
        increment inc;

        res.stress = qp.managedStateVars->stress;
        res.f.setZero();
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();

        inc.dLambda( 0 )        = lambdaIncrement;
        inc.laplaceDLambda( 0 ) = laplaceLambdaIncrement;
        inc.time                = time[1];
        inc.dT                  = dT;

        Voigt  stress       = Voigt::Zero();
        CSized C            = CSized::Zero();
        Voigt  dSdLambda    = Voigt::Zero();
        Voigt  dSdLaplacian = Voigt::Zero();
        Matrix< double, 1, ParentGeometryElement::voigtSize >
          dFddE = Matrix< double, 1, ParentGeometryElement::voigtSize >::Zero();

        try {
          if constexpr ( nDim == 2 ) {
            Vector6d dE6 = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
            inc.dStrain  = dE6;

            if ( sectionType == SectionType::PlaneStress ) {
              qp.material->computePlaneStress( res, tan, inc );
              stress       = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
              C            = ContinuumMechanics::PlaneStress::getPlaneStressTangent( tan.dStressddStrain );
              dSdLambda    = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddLambda.col( 0 ) );
              dSdLaplacian = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddLaplacian.col( 0 ) );
              dFddE        = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dFddStrain.row( 0 ).transpose() )
                        .transpose();
            }
            else if ( sectionType == SectionType::PlaneStrain ) {
              qp.material->computeStress( res, tan, inc );
              stress       = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
              C            = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( tan.dStressddStrain );
              dSdLambda    = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddLambda.col( 0 ) );
              dSdLaplacian = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressddLaplacian.col( 0 ) );
              dFddE        = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dFddStrain.row( 0 ).transpose() )
                        .transpose();
            }
            else {
              throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
            }
          }
          else if ( nDim == 3 ) {
            if ( sectionType != SectionType::Solid )
              throw std::invalid_argument( "Invalid section type for 3D element, expected Solid" );

            inc.dStrain = dE;
            qp.material->computeStress( res, tan, inc );
            stress       = res.stress;
            C            = tan.dStressddStrain;
            dSdLambda    = tan.dStressddLambda.col( 0 );
            dSdLaplacian = tan.dStressddLaplacian.col( 0 );
            dFddE        = tan.dFddStrain.row( 0 );
          }
        }
        catch ( const StressUpdateFailed& ) {
          pNewDT = 0.25;
          return;
        }

        // cout all tangents for debugging
        using namespace std;
        // cout << "dSdE: " << endl << C << endl;
        // cout << "dSdLambda: " << endl << dSdLambda.transpose() << endl;
        // cout << "dSdLaplacian: " << endl << dSdLaplacian.transpose() << endl;
        // cout << "dFddE: " << endl << dFddE.transpose() << endl;

        const double fYield     = res.f( 0 );
        const double dFddLambda = tan.dFddLambda( 0, 0 );
        const double dFddLap    = tan.dFddLaplacian( 0, 0 );

        // fischer-burmeister formulation for yield function and plastic multiplier constraints
        const double mu                    = 1e-14;
        const double scale                 = 1e3;
        const double lambdaIncrementScaled = scale * lambdaIncrement;

        const double fb = std::sqrt( fYield * fYield + lambdaIncrementScaled * lambdaIncrementScaled + 2 * mu ) +
                          fYield - lambdaIncrementScaled;
        // Corrected Fischer-Burmeister derivatives
        const double common_denom = std::sqrt( fYield * fYield + lambdaIncrementScaled * lambdaIncrementScaled +
                                               2 * mu );

        const double dfb_df       = ( fYield / common_denom ) + 1.0;
        const double dfb_ddLambda = ( lambdaIncrementScaled / common_denom ) - 1.0;

        fU -= B.transpose() * stress * qp.J0xW;
        fL -= ( N.transpose() * fb - penalty * dNdX.transpose() * cConstraint ) * qp.J0xW;
        fG -= ( penalty * NVec.transpose() * cConstraint ) * qp.J0xW;

        KUU += B.transpose() * C * B * qp.J0xW;
        KUL += B.transpose() * dSdLambda * N * qp.J0xW;
        KUG += B.transpose() * dSdLaplacian * divOperator * qp.J0xW;

        KLU += N.transpose() * dfb_df * dFddE * B * qp.J0xW;
        KLL += ( N.transpose() * ( dfb_df * dFddLambda + dfb_ddLambda * scale ) * N +
                 penalty * dNdX.transpose() * dNdX ) *
               qp.J0xW;
        KLG += ( N.transpose() * dfb_df * dFddLap * divOperator - penalty * dNdX.transpose() * NVec ) * qp.J0xW;

        KGU += Matrix< double, sizeDoFG, sizeDoFU >::Zero();
        KGL += ( -penalty * NVec.transpose() * dNdX ) * qp.J0xW;
        KGG += ( penalty * NVec.transpose() * NVec ) * qp.J0xW;

        // // / 2. Corrected Residuals
        // fU -= B.transpose() * stress * qp.J0xW;
        // fL += ( N.transpose() * fb ) * qp.J0xW; // <-- REMOVED penalty force term here!
        // fG -= ( NVec.transpose() * cConstraint ) * qp.J0xW; // <-- Pure L2 projection

        // // 3. Corrected Tangents (Remove penalty terms from KLL and KLG)
        // KUU += B.transpose() * C * B * qp.J0xW;
        // KUL += B.transpose() * dSdLambda * N * qp.J0xW;
        // KUG += B.transpose() * dSdLaplacian * divOperator * qp.J0xW;

        // KLU -= N.transpose() * dfb_df * dFddE * B * qp.J0xW;
        // KLL -= ( N.transpose() * ( dfb_df * dFddLambda + dfb_ddLambda * scale ) * N ) * qp.J0xW; // <-- REMOVED
        // penalty KLG -= ( N.transpose() * dfb_df * dFddLap * divOperator ) * qp.J0xW;           // <-- REMOVED penalty

        // KGU += Matrix< double, sizeDoFG, sizeDoFU >::Zero();
        // KGL += ( -NVec.transpose() * dNdX ) * qp.J0xW; // <-- Pure coupling
        // KGG += ( NVec.transpose() * NVec ) * qp.J0xW;  // <-- Standard vector mass matrix

        qp.managedStateVars->stress              = res.stress;
        qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
        qp.managedStateVars
          ->strain += ContinuumMechanics::VoigtNotation::make3DVoigt< ParentGeometryElement::voigtSize >( dE );
      }
    }

    void computeYourselfExplicit( const double* QTotal_,
                                  const double* dQ_,
                                  double*       Pe_,
                                  const double* time,
                                  double        dT,
                                  double&       pNewDT )
    {
      Map< const RhsSized > QTotal( QTotal_ );
      Map< const RhsSized > dQ( dQ_ );
      Map< RhsSized >       Pe( Pe_ );

      const Ref< const USizedVector > qU( QTotal.head( sizeDoFU ) );
      const Ref< const LSizedVector > qL( QTotal.segment( sizeDoFU, sizeDoFL ) );
      const Ref< const GSizedVector > qG( QTotal.tail( sizeDoFG ) );

      const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
      const Ref< const LSizedVector > dQL( dQ.segment( sizeDoFU, sizeDoFL ) );
      const Ref< const GSizedVector > dQG( dQ.tail( sizeDoFG ) );

      Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
      Ref< LSizedVector > fL( Pe.segment( sizeDoFU, sizeDoFL ) );
      Ref< GSizedVector > fG( Pe.tail( sizeDoFG ) );

      using response  = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::response;
      using increment = typename MarmotMaterialGradientPlasticityHypoElastic< 1 >::increment;

      const double penalty = elementProperties.size() > 1 ? elementProperties[1] : 1.0;

      for ( auto& qp : qps ) {
        const BSized&    B    = qp.B;
        const NSized&    N    = qp.N;
        const dNdXiSized dNdX = qp.dNdX;
        const auto       NVec = localGeometryElement.NB( N );

        Matrix< double, 1, sizeDoFG > divOperator = Matrix< double, 1, sizeDoFG >::Zero();
        for ( int a = 0; a < nNodes; a++ ) {
          for ( int i = 0; i < nDim; i++ ) {
            divOperator( 0, a * nDim + i ) = dNdX( i, a );
          }
        }

        Voigt dE = B * dQU;

        const double lambdaIncrement        = ( N * dQL )( 0 );
        const double laplaceLambdaIncrement = divOperator * dQG;

        const Matrix< double, nDim, 1 > gradLambdaInc = NVec * dQG;
        const Matrix< double, nDim, 1 > cConstraint   = gradLambdaInc - dNdX * dQL;

        response  res;
        increment inc;

        res.stress = qp.managedStateVars->stress;
        res.f.setZero();
        res.elasticEnergyDensity = qp.managedStateVars->elasticStrainEnergy / qp.J0xW;
        res.dissipation          = qp.managedStateVars->dissipation / qp.J0xW;
        res.stateVars            = qp.managedStateVars->materialStateVars.data();

        inc.dLambda( 0 )        = lambdaIncrement;
        inc.laplaceDLambda( 0 ) = laplaceLambdaIncrement;
        inc.time                = time[1];
        inc.dT                  = dT;

        Voigt stress = Voigt::Zero();

        try {
          if constexpr ( nDim == 2 ) {
            Vector6d dE6 = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
            inc.dStrain  = dE6;

            if ( sectionType == SectionType::PlaneStress ) {
              qp.material->computePlaneStressExplicit( res, inc );
              stress = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            }
            else if ( sectionType == SectionType::PlaneStrain ) {
              qp.material->computeStressExplicit( res, inc );
              stress = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            }
            else {
              throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
            }
          }
          else if ( nDim == 3 ) {
            if ( sectionType != SectionType::Solid )
              throw std::invalid_argument( "Invalid section type for 3D element, expected Solid" );

            inc.dStrain = dE;
            qp.material->computeStressExplicit( res, inc );
            stress = res.stress;
          }
        }
        catch ( const StressUpdateFailed& ) {
          pNewDT = 0.25;
          return;
        }

        const double fYield = res.f( 0 );

        fU -= B.transpose() * stress * qp.J0xW;
        fL -= ( N.transpose() * fYield - penalty * dNdX.transpose() * cConstraint ) * qp.J0xW;
        fG -= ( penalty * NVec.transpose() * cConstraint ) * qp.J0xW;

        qp.managedStateVars->stress              = res.stress;
        qp.managedStateVars->elasticStrainEnergy = res.elasticEnergyDensity * qp.J0xW;
        qp.managedStateVars->dissipation         = res.dissipation * qp.J0xW;
        qp.managedStateVars->totalStrainEnergy   = ( res.elasticEnergyDensity + res.dissipation ) * qp.J0xW;
        qp.managedStateVars
          ->strain += ContinuumMechanics::VoigtNotation::make3DVoigt< ParentGeometryElement::voigtSize >( dE );
      }
    }

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 int                                 elementFace,
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

    void computeBodyForce( double*       P_,
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

    void computeConsistentInertia( double* M )
    {
      Map< KeSizedMatrix > Me( M );
      Me.setZero();

      for ( const auto& qp : qps ) {
        const auto                  N    = qp.N;
        const auto                  N_u  = localGeometryElement.NB( N );
        const auto                  N_g  = localGeometryElement.NB( N );
        const double                rho  = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
        const std::vector< double > etaV = qp.material->getNonlocalViscosity(
          qp.managedStateVars->materialStateVars.data() );
        const double eta = etaV.empty() ? 0.0 : etaV[0];

        Me.topLeftCorner( sizeDoFU, sizeDoFU ) += N_u.transpose() * N_u * qp.J0xW * rho;
        Me.block( sizeDoFU, sizeDoFU, sizeDoFL, sizeDoFL ) += N.transpose() * N * qp.J0xW * eta;
        Me.bottomRightCorner( sizeDoFG, sizeDoFG ) += N_g.transpose() * N_g * qp.J0xW * eta;
      }
    }

    void computeLumpedInertia( double* M )
    {
      Map< RhsSized > LMM( M );
      LMM.setZero();

      constexpr int nNodesLinear  = ( 1 << nDim );
      auto          linGeometryEl = MarmotGeometryElement< nDim, nNodesLinear >();

      for ( const auto& qp : qps ) {
        const auto N     = qp.N;
        const auto N_lin = linGeometryEl.N( qp.xi );

        VectorXd N_weighted = 0.5 * N;
        N_weighted.head( nNodesLinear ) += 0.5 * N_lin;

        const double                rho  = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
        const std::vector< double > etaV = qp.material->getNonlocalViscosity(
          qp.managedStateVars->materialStateVars.data() );
        const double eta = etaV.empty() ? 0.0 : etaV[0];

        const VectorXd mU = N_weighted * qp.J0xW * rho;
        const VectorXd mL = N_weighted * qp.J0xW * eta;
        const VectorXd mG = N_weighted * qp.J0xW * eta;

        for ( int i = 0; i < nNodes; i++ ) {
          for ( int d = 0; d < nDim; d++ )
            LMM( i * nDim + d ) += mU( i );

          LMM( sizeDoFU + i ) += mL( i );

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
        waveSpeedResponse.stress = qp.managedStateVars->stress;
        waveSpeedResponse.f.setZero();
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

            XiSized coordAtGauss = localGeometryElement.NB( localGeometryElement.N( qp.xi ) ) *
                                   localGeometryElement.coordinates;

            const double sigY1 = values[0];
            const double sigY2 = values[2];
            const double y1    = values[1];
            const double y2    = values[3];

            using namespace Math;
            qp.managedStateVars->stress( 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 );
            qp.managedStateVars->stress( 0 ) = values[4] * qp.managedStateVars->stress( 1 );
            qp.managedStateVars->stress( 2 ) = values[5] * qp.managedStateVars->stress( 1 );
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
      coordsMap = localGeometryElement.NB( localGeometryElement.N( centerXi ) ) * localGeometryElement.coordinates;
      return coords;
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints()
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

    int getNumberOfQuadraturePoints() { return qps.size(); }
  };

} // namespace Marmot::Elements