#include "Marmot/GradientPlasticityDisplacementFiniteElement.h"
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

  const static bool GPS4_isRegistered = MarmotElementFactory::
    registerElement( "GPS4",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool GPS8_isRegistered = MarmotElementFactory::
    registerElement( "GPS8",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool GPE4_isRegistered = MarmotElementFactory::
    registerElement( "GPE4",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool GPE8_isRegistered = MarmotElementFactory::
    registerElement( "GPE8",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool GPS8R_isRegistered = MarmotElementFactory::
    registerElement( "GPS8R",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool GPE8R_isRegistered = MarmotElementFactory::
    registerElement( "GPE8R",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool GP3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GP3D8",
                     makeFactoryFunction< GradientPlasticityDisplacementFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          GradientPlasticityDisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool GP3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GP3D20",
                     makeFactoryFunction<
                       GradientPlasticityDisplacementFiniteElement< 3, 20 >,
                       FullIntegration,
                       GradientPlasticityDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool GP3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GP3D20R",
                     makeFactoryFunction<
                       GradientPlasticityDisplacementFiniteElement< 3, 20 >,
                       ReducedIntegration,
                       GradientPlasticityDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

} // namespace Marmot::Elements::Registration
