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

#include "Marmot/C0PenaltyGradientPlasticityElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  // 2D Plane Stress elements
  const static bool PGCPS4_isRegistered = MarmotElementFactory::
    registerElement( "PGCPS4",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 4 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 4 >::PlaneStress >() );

  const static bool PGCPS8_isRegistered = MarmotElementFactory::
    registerElement( "PGCPS8",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 8 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 8 >::PlaneStress >() );

  const static bool PGCPS8R_isRegistered = MarmotElementFactory::
    registerElement( "PGCPS8R",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 8 >,
                                          ReducedIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 8 >::PlaneStress >() );

  // 2D Plane Strain elements
  const static bool PGCPE4_isRegistered = MarmotElementFactory::
    registerElement( "PGCPE4",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 4 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 4 >::PlaneStrain >() );

  const static bool PGCPE8_isRegistered = MarmotElementFactory::
    registerElement( "PGCPE8",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 8 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 8 >::PlaneStrain >() );

  const static bool PGCPE8R_isRegistered = MarmotElementFactory::
    registerElement( "PGCPE8R",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 2, 8 >,
                                          ReducedIntegration,
                                          C0PenaltyGradientPlasticityElement< 2, 8 >::PlaneStrain >() );

  // 3D Solid elements
  const static bool PGC3D8_isRegistered = MarmotElementFactory::
    registerElement( "PGC3D8",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 3, 8 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 3, 8 >::Solid >() );

  const static bool PGC3D20_isRegistered = MarmotElementFactory::
    registerElement( "PGC3D20",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 3, 20 >,
                                          FullIntegration,
                                          C0PenaltyGradientPlasticityElement< 3, 20 >::Solid >() );

  const static bool PGC3D20R_isRegistered = MarmotElementFactory::
    registerElement( "PGC3D20R",
                     makeFactoryFunction< C0PenaltyGradientPlasticityElement< 3, 20 >,
                                          ReducedIntegration,
                                          C0PenaltyGradientPlasticityElement< 3, 20 >::Solid >() );

} // namespace Marmot::Elements::Registration
