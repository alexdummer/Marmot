#include "Marmot/GeneralGradientEnhancedDisplacementFiniteElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace Marmot::Elements::Registration {

#define CONCAT( a, b ) a##b

  enum GeneralGradientEnhancedDisplacementFiniteElementCode {

    /* TAG explanation
     * XXXX
     * ||||___    4: type of element
     * |||____    3: active fields
     * ||_____    2: number of nodes
     * |______    1: number of nodes
     *
     * active fields:   1: displacement + nonlocal damage,
     *
     * type of element: 1: 1D full integration,
     *                  2: 2D full integration, plane stress
     *                  3: 3D full integration,
     *                  4: 1D red. integration,
     *                  5: 2D red. integration, plane stress
     *                  6: 3D red. integration
     *                  7: 2D full integration, plane strain
     *                  8: 2D red. integration, plane strain
     * */

    GGCPS4  = CONCAT( 1193, 412 ),
    GGCPS8  = CONCAT( 1193, 812 ),
    GGCPS8R = CONCAT( 1193, 815 ),

    // Plane Strain
    GGCPE4  = CONCAT( 1193, 417 ),
    GGCPE4R = CONCAT( 1193, 418 ),
    GGCPE8  = CONCAT( 1193, 817 ),
    GGCPE8R = CONCAT( 1193, 818 ),

    // Solid
    GGC3D4   = CONCAT( 1193, 413 ),
    GGC3D10  = CONCAT( 1193, 1013 ),
    GGC3D8   = CONCAT( 1193, 813 ),
    GGC3D8R  = CONCAT( 1193, 816 ),
    GGC3D20  = CONCAT( 1193, 2013 ),
    GGC3D20R = CONCAT( 1193, 2016 ),
    // solid with two damage fields
    G2GC3D8 = CONCAT( 1193, 823 ),

  };

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

  /* const static bool GGCPS4_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPS4", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPS4, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >, */
  /*                                         FullIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStress >()
   * ); */

  /* const static bool GGCPS8_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPS8", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPS8, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >, */
  /*                                         FullIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >()
   * ); */

  /* const static bool GGCPE4_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPE4", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPE4, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >, */
  /*                                         FullIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 4 >::PlaneStrain >()
   * ); */

  /* const static bool GGCPE8_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPE8", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPE8, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >, */
  /*                                         FullIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >()
   * ); */
  /* const static bool GGCPS8R_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPS8R", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPS8R, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >, */
  /*                                         ReducedIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStress >()
   * ); */
  /* const static bool GGCPE8R_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "GGCPE8R", */
  /*                    GeneralGradientEnhancedDisplacementFiniteElementCode::GGCPE8R, */
  /*                    makeFactoryFunction< GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >, */
  /*                                         ReducedIntegration, */
  /*                                         GeneralGradientEnhancedDisplacementFiniteElement< 2, 8 >::PlaneStrain >()
   * ); */

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

  const static bool GGC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GGC3D8",
                     GeneralGradientEnhancedDisplacementFiniteElementCode::GGC3D8,
                     []( int elementID ) -> MarmotElement* {
                       return new GeneralGradientEnhancedDisplacementFiniteElement<
                         3,
                         8 >( elementID,
                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                              GeneralGradientEnhancedDisplacementFiniteElement< 3, 8 >::SectionType::Solid );
                     } );

  const static bool GGC3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GGC3D20",
                     GeneralGradientEnhancedDisplacementFiniteElementCode::GGC3D20,
                     []( int elementID ) -> MarmotElement* {
                       return new GeneralGradientEnhancedDisplacementFiniteElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                               GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid );
                     } );

  const static bool GGC3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "GGC3D20R",
                     GeneralGradientEnhancedDisplacementFiniteElementCode::GGC3D20R,
                     []( int elementID ) -> MarmotElement* {
                       return new GeneralGradientEnhancedDisplacementFiniteElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                               GeneralGradientEnhancedDisplacementFiniteElement< 3, 20 >::SectionType::Solid );
                     } );

  const static bool G2GC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "G2GC3D8",
                     GeneralGradientEnhancedDisplacementFiniteElementCode::G2GC3D8,
                     []( int elementID ) -> MarmotElement* {
                       return new GeneralGradientEnhancedDisplacementFiniteElement<
                         3,
                         8,
                         2 >( elementID,
                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                              GeneralGradientEnhancedDisplacementFiniteElement< 3, 8, 2 >::SectionType::Solid );
                     } );

} // namespace Marmot::Elements::Registration
