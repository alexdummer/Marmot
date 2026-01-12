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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Fastor;

  /**
   * @class MarmotMaterialFiniteStrainSubstepped
   * @brief A decorator that applies time substepping with analytical tangent accumulation.
   * * @tparam BaseMaterialType The concrete material class to wrap.
   * * Properties: [nSubsteps, Prop1, Prop2, ...]
   */
  template < typename BaseMaterialType >
  class MarmotMaterialFiniteStrainSubstepped : public MarmotMaterialFiniteStrain {
  protected:
    std::unique_ptr< BaseMaterialType > baseMaterial;
    int                                 nSubsteps;

  public:
    MarmotMaterialFiniteStrainSubstepped( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ )
    {
      // 1. Parse Substepping parameters
      if ( nMaterialProperties_ < 1 ) {
        throw std::invalid_argument( "MarmotMaterialFiniteStrainSubstepped requires at least 1 property (nSubsteps)" );
      }
      this->nSubsteps = static_cast< int >( matProperties_[0] );
      if ( this->nSubsteps < 1 )
        this->nSubsteps = 1;

      // 2. Instantiate Base Material (shift properties by 1)
      baseMaterial = std::make_unique< BaseMaterialType >( matProperties_ + 1,
                                                           nMaterialProperties_ - 1,
                                                           materialNumber_ );

      // 3. Initialize Layout
      initializeStateLayout();
    }

    virtual ~MarmotMaterialFiniteStrainSubstepped() = default;

    void initializeStateLayout() override
    {

      // Add our state variable: Deformation Gradient at start of GLOBAL step (F_n)
      this->stateLayout.add( "Substepping_F_n", 9 );

      // Add Base Material state variables
      this->stateLayout.add( "materialstate", baseMaterial->getNumberOfRequiredStateVars() );

      this->stateLayout.finalize();
    }

    void initializeYourself( double* stateVars, int nStateVars ) override
    {
      // 1. Initialize Base Material part
      int baseVarsCount = baseMaterial->getNumberOfRequiredStateVars();

      baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "materialstate" ), baseVarsCount );

      // 2. Initialize F_n to Identity
      FastorStandardTensors::Tensor33d&
        Fn = this->stateLayout.getAs< FastorStandardTensors::Tensor33d& >( stateVars, "Substepping_F_n" );
      Fn.eye();
    }

    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override
    {
      using namespace Eigen;
      using namespace FastorStandardTensors;

      // print state vars
      /* std::cout << "State Vars: " << VectorXd::Map(response.stateVars, this->stateLayout.totalSize()) .transpose() <<
       * std::endl; */

      // 1. Recover F_n (State at t_n)
      Tensor33d&      Fn_ref = this->stateLayout.getAs< Tensor33d& >( response.stateVars, "Substepping_F_n" );
      const Tensor33d Fn     = Fn_ref;
      const Tensor33d Fn1    = deformation.F;

      // 2. Setup Analytical Accumulation
      int nBaseState = baseMaterial->getNumberOfRequiredStateVars();

      // J_accum = d(State_i)/d(F_global). Size: nState x 9
      MatrixXd J_accum = MatrixXd::Zero( nBaseState, 9 );

      // Substep delta T
      double dt_sub = timeIncrement.dT / static_cast< double >( nSubsteps );
      double t_curr = timeIncrement.time;

      // Temporary data structures
      ConstitutiveResponse< 3 > subResponse;
      subResponse.stateVars = this->stateLayout.getPtr( response.stateVars, "materialstate" );

      AlgorithmicModuli< 3 >                         subTangents;
      MarmotMaterialFiniteStrain::StateSensitivities subSensitivities;
      Deformation< 3 >                               subDef;

      // 3. Substepping Loop
      for ( int i = 1; i <= nSubsteps; ++i ) {
        double xi = static_cast< double >( i ) / static_cast< double >( nSubsteps );

        // Linear Interpolation: F(xi) = (1-xi)*F_n + xi*F_n1
        subDef.F = ( 1.0 - xi ) * Fn + xi * Fn1;

        TimeIncrement subTime = { t_curr, dt_sub };

        // Call Base Material with Sensitivities
        // This updates stateVars to step i and fills subTangents/subSensitivities
        baseMaterial->computeStressWithSensitivities( subResponse, subTangents, subSensitivities, subDef, subTime );

        // Derivative of interpolation w.r.t global F_{n+1}
        // F_sub = F_n + xi * (F_{n+1} - F_n)
        // d(F_sub)/d(F_{n+1}) = xi * I
        // Matrix form (9x9) is xi * Identity
        double dFsub_dFglobal = xi;

        /* std::cout << "Substep " << i << "/" << nSubsteps << ", xi=" << xi << ", F_sub=\n" << subDef.F << std::endl;
         */
        /* std::cout << "  tau_sub=\n" << subResponse.tau << std::endl; */
        /* std::cout << "  subTangents.dTau_dF=\n"; */
        /* Map<Matrix<double, 9, 9>> C_alg_sub(subTangents.dTau_dF.data()); */
        /* std::cout << C_alg_sub << std::endl; */

        /* std::cout << "  subSensitivities.dState_dF=\n" << subSensitivities.dState_dF << std::endl; */
        /* std::cout << "  subSensitivities.dState_dStateOld=\n" << subSensitivities.dState_dStateOld << std::endl; */

        if ( i < nSubsteps ) {
          // Update J_accum for the next step
          // J_new = dState/dF_sub * (xi*I) + dState/dStateOld * J_old
          if ( nBaseState > 0 ) {
            J_accum = subSensitivities.dState_dF * dFsub_dFglobal + subSensitivities.dState_dStateOld * J_accum;
          }
        }
        else {
          // Final Step: Compute Global Tangent
          // C_global = C_sub * (xi*I) + dStress/dStateOld * J_old

          // 1. Term 1: Local Tangent contribution
          // Map Fastor 4th order tensor to Eigen 9x9 matrix
          // Note: Fastor layout is Row-Major, Eigen Map treats it linearly.
          // Since dTau_dF from numerical diff was generated using linear layout, this is consistent.
          Map< Matrix< double, 9, 9 > > C_alg_sub( subTangents.dTau_dF.data() );

          MatrixXd Term1 = C_alg_sub * dFsub_dFglobal;

          // 2. Term 2: History contribution
          MatrixXd Term2 = MatrixXd::Zero( 9, 9 );
          if ( nBaseState > 0 ) {
            Term2 = subSensitivities.dStress_dStateOld * J_accum;
          }

          // 3. Combine Terms to get C_global
          // Note the transpose of Term2 due to index ordering differences
          // between Fastor and Eigen conventions.
          MatrixXd C_global = Term1 + Term2.transpose();

          // Copy back to output tensor
          std::memcpy( tangents.dTau_dF.data(), C_global.data(), 81 * sizeof( double ) );
        }

        t_curr += dt_sub;
      }

      // 4. Output results
      response.tau                  = subResponse.tau;
      response.elasticEnergyDensity = subResponse.elasticEnergyDensity;

      // 5. Update F_n in the state vector for the next convergence check or time step
      // We update the stored "Substepping_F_n" to the current F
      Fn_ref = deformation.F;
    }
  };

} // namespace Marmot::Materials
