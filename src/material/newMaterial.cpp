
#include <jem/util/Properties.h>

#include "HookeMaterial.h"
#include "OrthotropicMaterial.h"
#include "DamageMaterial.h"

using namespace jem;


//-----------------------------------------------------------------------
//   newMaterial
//-----------------------------------------------------------------------


Ref<Material>         newMaterial

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  Properties     matProps = props.getProps ( name );
  Properties     matConf  = conf.makeProps ( name );

  Ref<Material>  mat;
  String         type;
  int            dim;

  matProps.get ( type, "type" );
  matConf .set ( "type", type );

  matProps.get ( dim, "dim"   );
  matConf .set ( "dim", dim   );
  
  if      ( type == "Hooke" )
  {
    mat = newInstance<HookeMaterial> ( dim, globdat );
  }
  else if ( type == "OrthoHooke" )
  {
    mat = newInstance<OrthotropicMaterial> ( dim, globdat );
  }
  else if ( type == "Damage" )
  {
    mat = newInstance<DamageMaterial> ( dim, globdat );
  }
  else
  {
    matProps.propertyError (
      name,
      "invalid material type: " + type
    );
  }

  return mat;
}
