#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/ViscoelasticNeoHooke.h"

int main()
{

  using namespace Marmot::Materials;
  std::array< double, 6 > materialProperties_ = { 3500, 1500, 0.5, 1, 0.2, 10 };
  const double            nMaterialProperties = materialProperties_.size();
  const int               elLabel             = 1;

  ViscoelasticNeoHooke mat = ViscoelasticNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  const int nStateVars = mat.getNumberOfRequiredStateVars();
  mat.assignStateVars( new double[nStateVars], nStateVars );

  ViscoelasticNeoHooke::Deformation< 3 > def;
  ViscoelasticNeoHooke::TimeIncrement    timeInc = { 0, 0.001 };

  ViscoelasticNeoHooke::ConstitutiveResponse< 3 > response;
  ViscoelasticNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F( 0, 0 ) += 1e-4;
  /* def.F( 1, 1 ) += 1e-4; */
  /* def.F( 2, 2 ) += 1e-4; */

  mat.computeStress( response, tangent, def, timeInc );
  std::cout << "Kirchhoff Stress: \n" << response.tau << std::endl;

  mat.computeStress( response, tangent, def, { 0.1, 10 } );
  std::cout << "Kirchhoff Stress: \n" << response.tau << std::endl;

  return 0;
}
