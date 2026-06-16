#include "Marmot/C0GradientPlasticityFiniteElement.h"
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

  const static bool C0GPPE4_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE4",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool C0GPPE8R_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8R",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C0GPPE9R_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE9R",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 9 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 9 >::PlaneStrain >() );

  const static bool C0GPPE8_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C0GPPE9_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE9",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 9 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 9 >::PlaneStrain >() );

  const static bool C0GPPE8U9L9GR_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8U9L9GR",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8, 9, 9 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8, 9, 9 >::PlaneStrain >() );

  const static bool C0GPPE8U9L4GR_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8U9L4GR",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8, 9, 4 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8, 9, 4 >::PlaneStrain >() );
} // namespace Marmot::Elements::Registration