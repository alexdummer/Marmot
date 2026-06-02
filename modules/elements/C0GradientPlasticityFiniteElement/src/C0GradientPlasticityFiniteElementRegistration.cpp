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

  const static bool C0GPPS4_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPS4",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool C0GPPE4_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE4",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool C0GPPS8R_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPS8R",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool C0GPPE8R_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8R",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C0GPPE8_isRegistered = MarmotElementFactory::
    registerElement( "C0GPPE8",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 2, 8 >::PlaneStrain >() );

  const static bool C0GP3D8_isRegistered = MarmotElementFactory::
    registerElement( "C0GP3D8",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          C0GradientPlasticityFiniteElement< 3, 8 >::Solid >() );

  const static bool C0GP3D20R_isRegistered = MarmotElementFactory::
    registerElement( "C0GP3D20R",
                     makeFactoryFunction< C0GradientPlasticityFiniteElement< 3, 20 >,
                                          ReducedIntegration,
                                          C0GradientPlasticityFiniteElement< 3, 20 >::Solid >() );

} // namespace Marmot::Elements::Registration