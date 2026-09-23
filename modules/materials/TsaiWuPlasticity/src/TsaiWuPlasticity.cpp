#include "Marmot/TsaiWuPlasticity.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/TsaiWuPlasticityConstants.h"

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  TsaiWuPlasticityModel::TsaiWuPlasticityModel( const double* materialProperties,
                                                const int     nMaterialProperties,
                                                const int     materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      E1( this->materialProperties[0] ),
      E2( this->materialProperties[1] ),
      E3( this->materialProperties[2] ),
      nu12( this->materialProperties[3] ),
      nu23( this->materialProperties[4] ),
      nu13( this->materialProperties[5] ),
      G12( this->materialProperties[6] ),
      G23( this->materialProperties[7] ),
      G13( this->materialProperties[8] ),
      T1( this->materialProperties[9] ),
      C1( this->materialProperties[10] ),
      T2( this->materialProperties[11] ),
      C2( this->materialProperties[12] ),
      T3( this->materialProperties[13] ),
      C3( this->materialProperties[14] ),
      S12( this->materialProperties[15] ),
      S23( this->materialProperties[16] ),
      S13( this->materialProperties[17] ),
      HLin( this->materialProperties[18] ),
      deltaYieldStress( this->materialProperties[19] ),
      delta( this->materialProperties[20] ),
      direction1( { this->materialProperties[21], this->materialProperties[22], this->materialProperties[23] } ),
      direction2( { this->materialProperties[24], this->materialProperties[25], this->materialProperties[26] } )
  {
    initializeStateLayout();

    double F11 = 1.0 / ( T1 * C1 );
    double F22 = 1.0 / ( T2 * C2 );
    double F33 = 1.0 / ( T3 * C3 );
    double F44 = 1.0 / ( S23 * S23 );
    double F55 = 1.0 / ( S13 * S13 );
    double F66 = 1.0 / ( S12 * S12 );

    double F1 = 1.0 / T1 - 1.0 / C1;
    double F2 = 1.0 / T2 - 1.0 / C2;
    double F3 = 1.0 / T3 - 1.0 / C3;

    double F12 = -0.5 * std::sqrt( F11 * F22 );
    double F13 = -0.5 * std::sqrt( F11 * F33 );
    double F23 = -0.5 * std::sqrt( F22 * F33 );

    q      = Eigen::Matrix< double, 6, 1 >::Zero();
    q( 0 ) = F1;
    q( 1 ) = F2;
    q( 2 ) = F3;

    P         = Eigen::Matrix< double, 6, 6 >::Zero();
    P( 0, 0 ) = F11;
    P( 1, 1 ) = F22;
    P( 2, 2 ) = F33;
    P( 3, 3 ) = F66;
    P( 4, 4 ) = F44;
    P( 5, 5 ) = F55;

    P( 0, 1 ) = P( 1, 0 ) = F12;
    P( 0, 2 ) = P( 2, 0 ) = F13;
    P( 1, 2 ) = P( 2, 1 ) = F23;

    // Ellipticity checks as requested
    if ( F11 * F22 - F12 * F12 <= 0.0 ) {
      throw std::runtime_error( "TsaiWuPlasticity: Yield surface is not elliptic in 1-2 plane." );
    }
    if ( F22 * F33 - F23 * F23 <= 0.0 ) {
      throw std::runtime_error( "TsaiWuPlasticity: Yield surface is not elliptic in 2-3 plane." );
    }
    if ( F11 * F33 - F13 * F13 <= 0.0 ) {
      throw std::runtime_error( "TsaiWuPlasticity: Yield surface is not elliptic in 1-3 plane." );
    }

    /* double detP3x3 = P.block< 3, 3 >( 0, 0 ).determinant(); */
    /* if ( detP3x3 < 0.0 ) { */
    /*   std::cerr << "detP3x3 = " << detP3x3 << std::endl; */
    /*   std::cerr << "F11 = " << F11 << ", F22 = " << F22 << ", F33 = " << F33 << std::endl; */
    /*   std::cerr << "F12 = " << F12 << ", F13 = " << F13 << ", F23 = " << F23 << std::endl; */
    /*   std::cerr << "P: \n" << P.block< 3, 3 >( 0, 0 ) << std::endl; */
    /*   throw std::runtime_error( */
    /*     "TsaiWuPlasticity: Yield surface is not elliptic in 3D principal stress space (determinant < 0)." ); */
    /* } */
    Cel = ContinuumMechanics::Elasticity::Orthotropic::stiffnessTensor( E1, E2, E3, nu12, nu23, nu13, G12, G23, G13 );

    localCoordinateSystem = Marmot::Math::orthonormalCoordinateSystem( direction1, direction2 );

    CelGlobal = ContinuumMechanics::VoigtNotation::Transformations::
      transformStiffnessToGlobalSystem( Cel, localCoordinateSystem );
  }

  double TsaiWuPlasticityModel::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 28 ) {
      throw std::runtime_error(
        "TsaiWuPlasticityModel::getDensity: Not enough material properties given! nMaterialProperties < 28" );
    }
    return this->materialProperties[27];
  }

  void TsaiWuPlasticityModel::computeStress( state3D&        state,
                                             Matrix6d&       dStress_dStrain,
                                             const Vector6d& dStrain,
                                             const timeInfo& timeInfo ) const
  {
    mVector6d  nomStress( state.stress.data() );
    mMatrix6d  C( dStress_dStrain.data() );
    const auto dEGlobal = dStrain;

    using namespace Marmot::ContinuumMechanics::VoigtNotation;
    Vector6d dELocal = Transformations::transformStrainToLocalSystem( dEGlobal, localCoordinateSystem );

    if ( dEGlobal.isZero( 1e-14 ) ) {
      C = CelGlobal;
      return;
    }

    double&  kappa     = stateLayout.getAs< double& >( state.stateVars, "kappa" );
    Vector6d SLocalOld = Transformations::transformStressToLocalSystem( nomStress, localCoordinateSystem );

    auto fy = [&]( double kappa_ ) {
      return 1.0 + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };

    auto dfy_ddKappa = [&]( double kappa_ ) { return HLin + deltaYieldStress * delta * std::exp( -delta * kappa_ ); };

    auto compute_rho = [&]( const Vector6d& S ) {
      double qs  = q.dot( S );
      double sps = S.dot( P * S );
      return 0.5 * ( qs + std::sqrt( qs * qs + 4.0 * sps ) );
    };

    auto compute_n = [&]( const Vector6d& S ) {
      double qs      = q.dot( S );
      double sps     = S.dot( P * S );
      double deltaSq = qs * qs + 4.0 * sps;
      if ( deltaSq < 1e-16 )
        return Vector6d( 0.5 * q );
      return Vector6d( 0.5 * ( q + ( q * qs + 4.0 * P * S ) / std::sqrt( deltaSq ) ) );
    };

    auto compute_H = [&]( const Vector6d& S ) {
      double qs      = q.dot( S );
      double sps     = S.dot( P * S );
      double deltaSq = qs * qs + 4.0 * sps;
      if ( deltaSq < 1e-16 )
        return Matrix6d( ( 0.5 * q ) * q.transpose() );
      double   rootDelta = std::sqrt( deltaSq );
      Vector6d v         = q * qs + 4.0 * P * S;
      return Matrix6d(
        0.5 * ( ( q * q.transpose() + 4.0 * P ) / rootDelta - ( v * v.transpose() ) / ( rootDelta * deltaSq ) ) );
    };

    Vector6d trialStress = SLocalOld + Cel * dELocal;
    double   rhoTrial    = compute_rho( trialStress );

    if ( rhoTrial - fy( kappa ) >= 0.0 ) {
      // plastic step
      Vector6d S       = trialStress;
      double   dLambda = 0.0;
      int      counter = 0;

      Matrix6d JssInv;

      while ( true ) {
        double   rho     = compute_rho( S );
        Vector6d n       = compute_n( S );
        Matrix6d H       = compute_H( S );
        double   fy_val  = fy( kappa + dLambda );
        double   dfy_val = dfy_ddKappa( kappa + dLambda );

        Vector6d Rs = S - trialStress + dLambda * ( Cel * n );
        double   Rl = rho - fy_val;

        if ( std::abs( Rl ) < TsaiWuPlasticityConstants::innerNewtonTol &&
             Rs.norm() < TsaiWuPlasticityConstants::innerNewtonTol ) {
          break;
        }

        if ( counter >= TsaiWuPlasticityConstants::nMaxInnerNewtonCycles ) {
          throw std::runtime_error( "return mapping failed to converge in TsaiWuPlasticityModel::computeStress" );
        }

        Matrix6d    Jss = Matrix6d::Identity() + dLambda * Cel * H;
        Vector6d    Jsl = Cel * n;
        RowVector6d Jls = n.transpose();
        double      Jll = -dfy_val;

        JssInv            = Jss.colPivHouseholderQr().solve( Matrix6d::Identity() );
        double   denom    = ( Jls * JssInv * Jsl ).value() - Jll;
        double   ddLambda = ( Rl - ( Jls * ( JssInv * Rs ) ).value() ) / denom;
        Vector6d dS       = -JssInv * ( Rs + Jsl * ddLambda );

        S += dS;
        dLambda += ddLambda;
        counter++;
      }

      kappa += dLambda;
      Vector6d deltaStressLocal = S - SLocalOld;
      nomStress += Transformations::transformStressToGlobalSystem( deltaStressLocal, localCoordinateSystem );

      Matrix6d Calg    = JssInv * Cel;
      double   dfy_val = dfy_ddKappa( kappa );
      Vector6d n       = compute_n( S );
      double   bottom  = ( n.transpose() * Calg * n ).value() + dfy_val;
      Matrix6d C_loc   = Calg - ( Calg * n * n.transpose() * Calg ) / bottom;

      C = Transformations::transformStiffnessToGlobalSystem( C_loc, localCoordinateSystem );
    }
    else {
      // elastic step
      Vector6d deltaStressLocal = Cel * dELocal;
      nomStress += Transformations::transformStressToGlobalSystem( deltaStressLocal, localCoordinateSystem );
      C = CelGlobal;
    }
  }

} // namespace Marmot::Materials
