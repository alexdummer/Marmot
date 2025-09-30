/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
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
#include <Eigen/Core>
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <string>

namespace Marmot::Materials {

  class ViscoelasticNeoHooke : public MarmotMaterialFiniteStrain {
  public:
    using MarmotMaterialFiniteStrain::MarmotMaterialFiniteStrain;

    ViscoelasticNeoHooke( const double* materialProperties, int nMaterialProperties, int materialLabel );

    void computeStress( ConstitutiveResponse< 3 >&,
                        AlgorithmicModuli< 3 >&,
                        const Deformation< 3 >&,
                        const TimeIncrement& );


    class ViscoelasticNeoHookeStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "S0", .length = 9 },
        { .name = "maxwellStateVars", .length = 0 },
      } );

      Fastor::TensorMap< double, 3, 3 > S0;
      double* maxwellStateVars;

    ViscoelasticNeoHookeStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), 
        S0( &find("S0")  ), 
        maxwellStateVars( &find( "maxwellStateVars" ) ){};
    };
    std::unique_ptr< ViscoelasticNeoHookeStateVarManager > stateVars;

    int getNumberOfRequiredStateVars() { return ViscoelasticNeoHookeStateVarManager::layout.nRequiredStateVars + nMaxwellElements * 9; };

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& result );
  protected:
    const double K; // bulk modulus
    const double G; // shear modulus
    
    const size_t nMaxwellElements;

    Eigen::VectorXd tau;
    Eigen::VectorXd gamma;
  };

} // namespace Marmot::Materials
