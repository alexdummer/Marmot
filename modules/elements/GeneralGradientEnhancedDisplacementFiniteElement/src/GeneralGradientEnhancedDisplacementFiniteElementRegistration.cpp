#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

#undef CONCAT

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool GCPS4_isRegistered = MarmotElementFactory::
    registerElement( "GCPS4",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool GCPS8_isRegistered = MarmotElementFactory::
    registerElement( "GCPS8",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >() );

  const static bool GCPE4_isRegistered = MarmotElementFactory::
    registerElement( "GCPE4",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool GCPE8_isRegistered = MarmotElementFactory::
    registerElement( "GCPE8",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          FullIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );
  const static bool GCPS8R_isRegistered = MarmotElementFactory::
    registerElement( "GCPS8R",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >() );
  const static bool GCPE8R_isRegistered = MarmotElementFactory::
    registerElement( "GCPE8R",
                     makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >,
                                          ReducedIntegration,
                                          GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >() );

  /* const static bool GGC3D4_isRegistered = MarmotLibrary::MarmotElementFactory:: */
  /*   registerElement( "GGC3D4", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGC3D4, */
  /*                    []( int elementID ) -> MarmotElement* { */
  /*                      return new GeneralGradientEnhancedDisplacementFiniteElement< */
  /*                        3, */
  /*                        4 >( elementID, */
  /*                             Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, */
  /*                             GeneralGradientEnhancedDisplacementFiniteElement< 3, 4 >::SectionType::Solid ); */
  /*                    } ); */
  /* const static bool GGC3D10_isRegistered = MarmotLibrary::MarmotElementFactory:: */
  /*   registerElement( "GGC3D10", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGC3D10, */
  /*                    []( int elementID ) -> MarmotElement* { */
  /*                      return new GeneralGradientEnhancedDisplacementFiniteElement< */
  /*                        3, */
  /*                        10 >( elementID, */
  /*                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, */
  /*                              GeneralGradientEnhancedDisplacementFiniteElement< 3, 10 >::SectionType::Solid ); */
  /*                    } ); */

  const static bool GC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8 >::SectionType::Solid >() );

  const static bool GC3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool GC3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >,
                       ReducedIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid >() );

  const static bool G2GCPS4_isRegistered = MarmotElementFactory::
    registerElement( "G2GCPS4",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 4, 2 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 4, 2 >::PlaneStress >() );

  const static bool G2GCPS8_isRegistered = MarmotElementFactory::
    registerElement( "G2GCPS8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >::PlaneStress >() );

  const static bool G2GCPS8R_isRegistered = MarmotElementFactory::
    registerElement( "G2GCPS8R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >,
                       ReducedIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 2 >::PlaneStress >() );

  const static bool G2GC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G2GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >::SectionType::Solid >() );

  // const static bool G2GC3D8R_isRegistered = MarmotLibrary::MarmotElementFactory::
  //   registerElement( "G2GC3D8R",
  //                    GeneralGradientEnhancedDisplacementFiniteElementCode::G2GC3D8R,
  //                    []( int elementID ) -> MarmotElement* {
  //                      return new GeneralGradientEnhancedDisplacementFiniteElement<
  //                        3,
  //                        8,
  //                        2 >( elementID,
  //                             Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
  //                             GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >::SectionType::Solid );
  //                    } );

  const static bool G2GC3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G2GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >::SectionType::Solid >() );

  const static bool G2GC3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G2GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >,
                       ReducedIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 2 >::SectionType::Solid >() );

  const static bool G6GCPS8_isRegistered = MarmotElementFactory::
    registerElement( "G6GCPS8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6 >::PlaneStress >() );

  const static bool G6GC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G6GC3D8",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 6 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 6 >::SectionType::Solid >() );

  const static bool G6GC3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G6GC3D20",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >::SectionType::Solid >() );

  const static bool G6GC3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G6GC3D20R",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >,
                       ReducedIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6 >::SectionType::Solid >() );

  const static bool G6GCPS84_isRegistered = MarmotElementFactory::
    registerElement( "G6GCPS84",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6, 4 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 2, 8, 6, 4 >::PlaneStress >() );

  const static bool G6GC3D20R8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G6GC3D20RM",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >,
                       ReducedIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >::SectionType::Solid >() );

  const static bool G6GC3D208_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G6GC3D20M",
                     makeFactoryFunction<
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >,
                       FullIntegration,
                       GeneralGradientEnhancedDisplacementFiniteElement< 3, 20, 6, 8 >::SectionType::Solid >() );

} // namespace Marmot::Elements::Registration
