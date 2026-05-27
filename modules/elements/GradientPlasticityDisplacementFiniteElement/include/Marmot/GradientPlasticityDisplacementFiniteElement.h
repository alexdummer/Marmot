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
#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotMaterialGradientPlasticityHypoElasticFactory.h"

namespace Marmot::Elements {

  template < int nDim, int nNodes, int nNonLocalNodes = nNodes >
  class GradientPlasticityDisplacementFiniteElement
    : public GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, 1, nNonLocalNodes > {
  public:
    using BaseType    = GeneralGradientEnhancedDisplacementFiniteElement< nDim, nNodes, 1, nNonLocalNodes >;
    using SectionType = typename BaseType::SectionType;

    GradientPlasticityDisplacementFiniteElement( int                                         elementID,
                                                 FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                 SectionType                                 sectionType )
      : BaseType( elementID, integrationType, sectionType )
    {
    }

    void assignProperty( const MarmotMaterialSection& section ) override
    {
      for ( auto& qp : this->qps ) {
        qp.material = std::unique_ptr< MarmotMaterialGeneralGradientEnhancedHypoElastic< 1 > >(
          MarmotLibrary::MarmotMaterialGradientPlasticityHypoElasticFactory::createMaterial(
            section.materialName, section.materialProperties, section.nMaterialProperties, this->elLabel ) );
      }
    }
  };

} // namespace Marmot::Elements
