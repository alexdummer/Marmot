#include "Marmot/MarmotMaterialGradientPlasticityHypoElasticFactory.h"
#include "Marmot/MarmotJournal.h"

using namespace MarmotLibrary;

template <>
MarmotMaterialGradientPlasticityHypoElasticFactory< 1 >::MaterialFactoryMap& MarmotMaterialGradientPlasticityHypoElasticFactory<
  1 >::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}

template <>
MarmotMaterialGradientPlasticityHypoElastic< 1 >* MarmotMaterialGradientPlasticityHypoElasticFactory<
  1 >::createMaterial( const std::string& materialName,
                       const double*      materialProperties,
                       int                nMaterialProperties,
                       int                materialNumber )
{
  auto& map = materialFactoryFunctionByName();
  auto  it  = map.find( materialName );
  if ( it == map.end() ) {
    std::string reg = "Registered materials are: ";
    for ( const auto& pair : map ) {
      reg += pair.first + ", ";
    }
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__ << " Material " + materialName + " not registered!" + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}