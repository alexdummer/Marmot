#include "Marmot/PhaseFieldDisplacementFiniteElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace Marmot::Elements::Registration {

#define CONCAT( a, b ) a##b

  enum PhaseFieldDisplacementFiniteElementCode {

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

    PFCPS4  = CONCAT( 1193, 412 ),
    PFCPS8  = CONCAT( 1193, 812 ),
    PFCPS8R = CONCAT( 1193, 815 ),

    // Plane Strain
    PFCPE4  = CONCAT( 1193, 417 ),
    PFCPE4R = CONCAT( 1193, 418 ),
    PFCPE8  = CONCAT( 1193, 817 ),
    PFCPE8R = CONCAT( 1193, 818 ),

    // Solid
    PFC3D4   = CONCAT( 1193, 413 ),
    PFC3D10  = CONCAT( 1193, 1013 ),
    PFC3D8   = CONCAT( 1193, 813 ),
    PFC3D8R  = CONCAT( 1193, 816 ),
    PFC3D20  = CONCAT( 1193, 2013 ),
    PFC3D20R = CONCAT( 1193, 2016 ),

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

  /* const static bool PFCPS4_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPS4", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPS4, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 4 >, */
  /*                                         FullIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 4 >::PlaneStress >() ); */

  /* const static bool PFCPS8_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPS8", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPS8, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 8 >, */
  /*                                         FullIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 8 >::PlaneStress >() ); */

  /* const static bool PFCPE4_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPE4", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPE4, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 4 >, */
  /*                                         FullIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 4 >::PlaneStrain >() ); */

  /* const static bool PFCPE8_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPE8", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPE8, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 8 >, */
  /*                                         FullIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 8 >::PlaneStrain >() ); */
  /* const static bool PFCPS8R_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPS8R", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPS8R, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 8 >, */
  /*                                         ReducedIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 8 >::PlaneStress >() ); */
  /* const static bool PFCPE8R_isRegistered = MarmotElementFactory:: */
  /*   registerElement( "PFCPE8R", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFCPE8R, */
  /*                    makeFactoryFunction< PhaseFieldDisplacementFiniteElement< 2, 8 >, */
  /*                                         ReducedIntegration, */
  /*                                         PhaseFieldDisplacementFiniteElement< 2, 8 >::PlaneStrain >() ); */

  /* const static bool PFC3D4_isRegistered = MarmotLibrary::MarmotElementFactory:: */
  /*   registerElement( "PFC3D4", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFC3D4, */
  /*                    []( int elementID ) -> MarmotElement* { */
  /*                      return new PhaseFieldDisplacementFiniteElement< */
  /*                        3, */
  /*                        4 >( elementID, */
  /*                             Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, */
  /*                             PhaseFieldDisplacementFiniteElement< 3, 4 >::SectionType::Solid ); */
  /*                    } ); */
  /* const static bool PFC3D10_isRegistered = MarmotLibrary::MarmotElementFactory:: */
  /*   registerElement( "PFC3D10", */
  /*                    PhaseFieldDisplacementFiniteElementCode::PFC3D10, */
  /*                    []( int elementID ) -> MarmotElement* { */
  /*                      return new PhaseFieldDisplacementFiniteElement< */
  /*                        3, */
  /*                        10 >( elementID, */
  /*                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, */
  /*                              PhaseFieldDisplacementFiniteElement< 3, 10 >::SectionType::Solid ); */
  /*                    } ); */

  const static bool PFC3D8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "PFC3D8",
                     PhaseFieldDisplacementFiniteElementCode::PFC3D8,
                     []( int elementID ) -> MarmotElement* {
                       return new PhaseFieldDisplacementFiniteElement<
                         3,
                         8 >( elementID,
                              Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                              PhaseFieldDisplacementFiniteElement< 3, 8 >::SectionType::Solid );
                     } );

  const static bool PFC3D20_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "PFC3D20",
                     PhaseFieldDisplacementFiniteElementCode::PFC3D20,
                     []( int elementID ) -> MarmotElement* {
                       return new PhaseFieldDisplacementFiniteElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                               PhaseFieldDisplacementFiniteElement< 3, 20 >::SectionType::Solid );
                     } );

  const static bool PFC3D20R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "PFC3D20R",
                     PhaseFieldDisplacementFiniteElementCode::PFC3D20R,
                     []( int elementID ) -> MarmotElement* {
                       return new PhaseFieldDisplacementFiniteElement<
                         3,
                         20 >( elementID,
                               Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                               PhaseFieldDisplacementFiniteElement< 3, 20 >::SectionType::Solid );
                     } );

} // namespace Marmot::Elements::Registration
