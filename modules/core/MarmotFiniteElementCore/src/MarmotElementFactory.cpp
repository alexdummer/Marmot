#include "Marmot/MarmotElementFactory.h"

using namespace Marmot::Registration;

MarmotElementFactory::ElementFactoryMap& MarmotElementFactory::elementFactoryFunctionByName()
{
  static ElementFactoryMap map;
  return map;
}
