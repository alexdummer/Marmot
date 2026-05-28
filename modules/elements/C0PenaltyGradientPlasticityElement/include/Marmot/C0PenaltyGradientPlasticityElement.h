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
 * Modified for C0-Continuous Penalty-Enhanced Gradient Plasticity
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMaterialC0PenaltyGradientPlasticity.h"
#include "Marmot/MarmotMaterialC0PenaltyGradientPlasticityFactory.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::C0PenaltyGradientPlasticityElement
   * @tparam nDim Number of spatial dimensions (1, 2, or 3).
   * @tparam nNodes Number of element nodes.
   * @brief C0-continuous penalty-enhanced gradient plasticity finite element.
   *
   * @details This element implements the C0-continuous penalty-enhanced gradient
   * plasticity formulation. It uses C0 shape functions for both the displacement
   * field and the nonlocal cumulative plastic strain field \f$ \bar{\kappa} \f$.
   *
   * The element has the following DOFs per node:
   *  - nDim displacement DOFs
   *  - 1 nonlocal cumulative plastic strain DOF \f$ \bar{\kappa} \f$
   *
   * The weak form consists of:
   * 1. Balance of linear momentum (standard):
   * \f[
   *   \int_\Omega \mathbf{B}^\mathsf{T} \boldsymbol{\sigma} \, \mathrm{d}\Omega = \mathbf{f}_{\mathrm{ext}}
   * \f]
   *
   * 2. Penalty-enhanced gradient equation:
   * \f[
   *   \int_\Omega \mathbf{N}_\kappa^\mathsf{T} \beta (\bar{\kappa} - \kappa) \, \mathrm{d}\Omega
   *   + \int_\Omega \nabla \mathbf{N}_\kappa^\mathsf{T} \, l^2 \, \nabla \bar{\kappa} \, \mathrm{d}\Omega = 0
   * \f]
   *
   * where \f$ \beta \f$ is the penalty parameter and \f$ l \f$ is the internal length scale.
   * Both \f$ \beta \f$ and \f$ l \f$ are provided as element properties.
   */
  template < int nDim, int nNodes >
  class C0PenaltyGradientPlasticityElement : public MarmotElement, public MarmotGeometryElement< nDim, nNodes > {

  public:
    enum SectionType {
      PlaneStress,
      PlaneStrain,
      Solid,
    };

    static constexpr int nDofPerNodeU = nDim;
    static constexpr int nDofPerNodeK = 1; ///< Scalar nonlocal cumulative plastic strain

    static constexpr int sizeLoadVector = nNodes * ( nDofPerNodeU + nDofPerNodeK );
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int sizeDoFU = nNodes * nDofPerNodeU;
    static constexpr int sizeDoFK = nNodes * nDofPerNodeK;

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

    /** Element-level properties. For 2D: [thickness, penalty, lengthScale].
     *  For 3D: [penalty, lengthScale]. */
    Map< const VectorXd > elementProperties;
    /** Element label (ID). */
    const int elLabel;
    /** Section assumption. */
    const SectionType sectionType;

    /**
     * @brief Data and state associated with a quadrature point.
     */
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

      std::unique_ptr< MarmotMaterialC0PenaltyGradientPlasticity > material;

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
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdX( dNdXiSized::Zero() ),
          B( BSized::Zero() ){};
    };

    /// Quadrature points owned by the element.
    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct element with ID, quadrature rule and section assumption.
     */
    C0PenaltyGradientPlasticityElement( int                                         elementID,
                                        FiniteElement::Quadrature::IntegrationTypes integrationType,
                                        SectionType                                 sectionType );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

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

    void computeBodyForce( double*       P,
                           double*       K,
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

    void computeConsistentInertia( double* M );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( stateName == "sdv" ) {
        MarmotJournal::warningToMSG( MakeString()
                                     << __PRETTY_FUNCTION__
                                     << " on 'sdv' is discouraged and deprecated, please use precise state name" );
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

    std::vector< double >                getCoordinatesAtCenter();
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();
    int                                  getNumberOfQuadraturePoints();
  };

  // ---------------------------------------------------------------------------
  //  Implementation
  // ---------------------------------------------------------------------------

  template < int nDim, int nNodes >
  C0PenaltyGradientPlasticityElement< nDim, nNodes >::C0PenaltyGradientPlasticityElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : elementProperties( nullptr, 0 ), elLabel( elementID ), sectionType( sectionType )
  {
    for ( const auto& qpInfo :
          FiniteElement::Quadrature::getGaussPointInfo( ParentGeometryElement::shape, integrationType ) ) {
      qps.emplace_back( qpInfo.xi, qpInfo.weight );
    }
  }

  template < int nDim, int nNodes >
  int C0PenaltyGradientPlasticityElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties )
      Map< const VectorXd >( elementPropertiesInfo.elementProperties, elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< MarmotMaterialC0PenaltyGradientPlasticity >(
        MarmotLibrary::MarmotMaterialC0PenaltyGradientPlasticityFactory::createMaterial( section.materialName,
                                                                                         section.materialProperties,
                                                                                         section.nMaterialProperties,
                                                                                         elLabel ) );
    }
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > C0PenaltyGradientPlasticityElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        nodeFields[i].push_back( "nonlocal plastic strain" );
      }
    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > C0PenaltyGradientPlasticityElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      // Displacement DOFs first (grouped by node)
      for ( int i = 0; i < nNodes; i++ )
        for ( int j = 0; j < nDim; j++ )
          permutationPattern.push_back( i * ( nDim + 1 ) + j );
      // Then nonlocal plastic strain DOFs (one per node)
      for ( int i = 0; i < nNodes; i++ )
        permutationPattern.push_back( i * ( nDim + 1 ) + nDim );
    }
    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::initializeYourself()
  {
    for ( QuadraturePoint& qp : qps ) {
      const auto          dNdXi = ParentGeometryElement::dNdXi( qp.xi );
      const JacobianSized J     = ParentGeometryElement::Jacobian( dNdXi );
      const JacobianSized JInv  = J.inverse();
      qp.detJ                   = J.determinant();
      qp.N                      = ParentGeometryElement::N( qp.xi );
      qp.dNdX                   = ParentGeometryElement::dNdX( dNdXi, JInv );
      qp.B                      = ParentGeometryElement::B( qp.dNdX );

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

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::computeYourself( const double* QTotal_,
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

    // Extract displacement and nonlocal plastic strain DOFs
    const Ref< const USizedVector > dQU( dQ.head( sizeDoFU ) );
    const Ref< const KSizedVector > dQK( dQ.tail( sizeDoFK ) );
    const Ref< const KSizedVector > qK( QTotal.tail( sizeDoFK ) );

    // Sub-stiffness matrices
    Ref< Matrix< double, sizeDoFU, sizeDoFU > > kUU( Ke.topLeftCorner( sizeDoFU, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFU, sizeDoFK > > kUK( Ke.topRightCorner( sizeDoFU, sizeDoFK ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFU > > kKU( Ke.bottomLeftCorner( sizeDoFK, sizeDoFU ) );
    Ref< Matrix< double, sizeDoFK, sizeDoFK > > kKK( Ke.bottomRightCorner( sizeDoFK, sizeDoFK ) );

    // Right-hand side
    Ref< USizedVector > fU( Pe.head( sizeDoFU ) );
    Ref< KSizedVector > fK( Pe.tail( sizeDoFK ) );

    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    // Element properties: penalty parameter and length scale
    double penalty, lengthScale;
    if constexpr ( nDim == 2 ) {
      penalty     = elementProperties[1];
      lengthScale = elementProperties[2];
    }
    else if constexpr ( nDim == 3 ) {
      penalty     = elementProperties[0];
      lengthScale = elementProperties[1];
    }
    else {
      // 1D
      penalty     = elementProperties[1];
      lengthScale = elementProperties[2];
    }

    const double l2 = lengthScale * lengthScale;

    using response  = typename MarmotMaterialC0PenaltyGradientPlasticity::response;
    using tangents  = typename MarmotMaterialC0PenaltyGradientPlasticity::tangents;
    using increment = typename MarmotMaterialC0PenaltyGradientPlasticity::increment;

    for ( size_t i = 0; i < this->qps.size(); i++ ) {

      QuadraturePoint& qp = this->qps[i];

      const BSized& B = qp.B;
      const NSized& N = qp.N;

      Voigt dE = B * dQU;

      // Interpolate nonlocal plastic strain at Gauss point
      double kappaNL  = N * qK;
      double dKappaNL = N * dQK;

      response  res;
      tangents  tan;
      increment inc;
      try {
        if constexpr ( nDim == 2 ) {
          Vector6d dE6  = ContinuumMechanics::VoigtNotation::planeVoigtToVoigt( dE );
          res.stress    = qp.managedStateVars->stress;
          res.stateVars = qp.managedStateVars->materialStateVars.data();
          inc           = { dE6, kappaNL, dKappaNL, time[1], dT };
          CSized C      = CSized::Zero();
          Voigt  S      = Voigt::Zero();

          if ( sectionType == SectionType::PlaneStress ) {
            qp.material->computePlaneStress( res, tan, inc );
            S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            C = ContinuumMechanics::PlaneStress::getPlaneStressTangent( tan.dStressDDStrain );
          }
          else if ( sectionType == SectionType::PlaneStrain ) {
            qp.material->computeStress( res, tan, inc );
            S = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( res.stress );
            C = ContinuumMechanics::PlaneStrain::getPlaneStrainTangent( tan.dStressDDStrain );
          }
          else {
            throw std::invalid_argument( "Invalid section type for 2D element, expected PlaneStress or PlaneStrain" );
          }

          // Momentum balance: fU = -∫ B^T σ dΩ
          fU -= B.transpose() * S * qp.J0xW;
          // kUU = ∫ B^T C B dΩ
          kUU += B.transpose() * C * B * qp.J0xW;

          // Gradient plasticity equation:
          // fK = -∫ [ N^T β(κ̄ - κ) + l² ∇N^T ∇κ̄ ] dΩ
          fK -= ( N.transpose() * penalty * ( kappaNL - res.kappaLocal ) + l2 * qp.dNdX.transpose() * qp.dNdX * qK ) *
                qp.J0xW;

          const auto dSdK            = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt( tan.dStressDDKappaNL );
          const auto dKappaLocal_dDE = ContinuumMechanics::VoigtNotation::voigtToPlaneVoigt(
            Eigen::Matrix< double, 6, 1 >( tan.dKappaLocalDDStrain.transpose() ) );

          // kUK = ∫ B^T (dσ/dκ̄) N dΩ
          kUK += B.transpose() * dSdK * N * qp.J0xW;

          // kKU = -∫ N^T β (dκ/dε) B dΩ
          kKU -= N.transpose() * penalty * dKappaLocal_dDE.transpose() * B * qp.J0xW;

          // kKK = ∫ [ N^T β (1 - dκ/dκ̄) N + l² ∇N^T ∇N ] dΩ
          kKK += ( N.transpose() * penalty * ( 1.0 - tan.dKappaLocalDDKappaNL ) * N +
                   l2 * qp.dNdX.transpose() * qp.dNdX ) *
                 qp.J0xW;
        }

        else if constexpr ( nDim == 3 ) {
          if ( sectionType != Solid ) {
            throw std::invalid_argument( "Invalid section type for 3D element! Must be Solid" );
          }

          res.stress    = qp.managedStateVars->stress;
          res.stateVars = qp.managedStateVars->materialStateVars.data();
          inc           = { dE, kappaNL, dKappaNL, time[1], dT };
          qp.material->computeStress( res, tan, inc );

          // Momentum balance
          fU -= B.transpose() * res.stress * qp.J0xW;
          kUU += B.transpose() * tan.dStressDDStrain * B * qp.J0xW;

          // Gradient plasticity equation
          fK -= ( N.transpose() * penalty * ( kappaNL - res.kappaLocal ) + l2 * qp.dNdX.transpose() * qp.dNdX * qK ) *
                qp.J0xW;

          // kUK = ∫ B^T (dσ/dκ̄) N dΩ
          kUK += B.transpose() * tan.dStressDDKappaNL * N * qp.J0xW;

          // kKU = -∫ N^T β (dκ/dε) B dΩ
          kKU -= N.transpose() * penalty * tan.dKappaLocalDDStrain * B * qp.J0xW;

          // kKK = ∫ [ N^T β (1 - dκ/dκ̄) N + l² ∇N^T ∇N ] dΩ
          kKK += ( N.transpose() * penalty * ( 1.0 - tan.dKappaLocalDDKappaNL ) * N +
                   l2 * qp.dNdX.transpose() * qp.dNdX ) *
                 qp.J0xW;
        }
      }
      catch ( StressUpdateFailed& e ) {
        pNewDT = 0.25;
        return;
      }
      qp.managedStateVars->stress = res.stress;
      qp.managedStateVars->strain += make3DVoigt< ParentGeometryElement::voigtSize >( dE );
    }
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Map< KeSizedMatrix > Me( M );
    Me.setZero();

    for ( const auto& qp : qps ) {
      const auto   N_  = ParentGeometryElement::NB( ParentGeometryElement::N( qp.xi ) );
      const double rho = qp.material->getDensity( qp.managedStateVars->materialStateVars.data() );
      Me.topLeftCorner( sizeDoFU, sizeDoFU ) += N_.transpose() * N_ * qp.detJ * qp.J0xW * rho;
    }
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
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

      FiniteElement::BoundaryElement boundaryEl( ParentGeometryElement::shape,
                                                 elementFace,
                                                 nDim,
                                                 ParentGeometryElement::coordinates );

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

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                             double*       K,
                                                                             const double* load,
                                                                             const double* QTotal,
                                                                             const double* time,
                                                                             double        dT )
  {
    Map< RhsSized >                        Pe( P );
    Ref< USizedVector >                    fU( Pe.head( sizeDoFU ) );
    Map< const Matrix< double, nDim, 1 > > f( load );

    for ( const auto& qp : qps ) {
      const auto N_ = ParentGeometryElement::NB( ParentGeometryElement::N( qp.xi ) );
      fU += N_.transpose() * f * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void C0PenaltyGradientPlasticityElement< nDim, nNodes >::setInitialConditions( StateTypes    state,
                                                                                 const double* values )
  {
    if constexpr ( nDim > 1 ) {
      switch ( state ) {
      case MarmotElement::GeostaticStress: {
        for ( QuadraturePoint& qp : qps ) {

          XiSized coordAtGauss = ParentGeometryElement::NB( ParentGeometryElement::N( qp.xi ) ) *
                                 ParentGeometryElement::coordinates;

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
      default: break;
      }
    }
  }

  template < int nDim, int nNodes >
  std::vector< double > C0PenaltyGradientPlasticityElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    coordsMap = ParentGeometryElement::NB( ParentGeometryElement::N( centerXi ) ) * ParentGeometryElement::coordinates;
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > C0PenaltyGradientPlasticityElement< nDim,
                                                                           nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      coordsMap = ParentGeometryElement::NB( ParentGeometryElement::N( qp.xi ) ) * ParentGeometryElement::coordinates;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int C0PenaltyGradientPlasticityElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

} // namespace Marmot::Elements
