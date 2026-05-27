#include "Marmot/GradientPlasticityDisplacementFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTesting.h"

using namespace Marmot;
using namespace Marmot::Elements;
using namespace Marmot::Testing;

void testBasicProperties()
{
  using ElemType = GradientPlasticityDisplacementFiniteElement< 2, 4 >;
  auto element   = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStrain );

  throwExceptionOnFailure( element->getNNodes() == 4, "wrong number of nodes" );
  throwExceptionOnFailure( element->getNSpatialDimensions() == 2, "wrong number of dimensions" );
  throwExceptionOnFailure( element->getNDofPerElement() == 12, "wrong number of dofs" );
}

void testAssignGradientPlasticityMaterial()
{
  using ElemType = GradientPlasticityDisplacementFiniteElement< 2, 4 >;
  auto element   = std::make_unique< ElemType >( 1,
                                               FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               ElemType::SectionType::PlaneStrain );

  const std::vector< double > coords = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
  element->assignNodeCoordinates( coords.data() );
  const std::vector< double > elementProps = { 1.0 };
  element->assignProperty( ElementProperties( elementProps.data(), elementProps.size() ) );

  const std::vector< double > materialProperties = { 210000.0, 0.3, 250.0, -2500.0, 0.05, 7850.0, 0.01 };
  element->assignProperty(
    MarmotMaterialSection( "GRADIENTVONMISES", materialProperties.data(), materialProperties.size() ) );
  element->initializeYourself();

  auto stateVars = std::vector< double >( element->getNumberOfRequiredStateVars(), 0.0 );
  element->assignStateVars( stateVars.data(), stateVars.size() );

  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 4, "unexpected number of quadrature points" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = { testBasicProperties, testAssignGradientPlasticityMaterial };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
