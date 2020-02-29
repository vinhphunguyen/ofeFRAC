
#include <jem/util/Properties.h>

#include "XuNeedlemanMaterial.h"
#include "TuronCohesiveMaterial.h"
#include "ContactMaterial.h"

using namespace jem;


//-----------------------------------------------------------------------
//   newCohesiveMaterial
//-----------------------------------------------------------------------


Ref<CohesiveMaterial>   newCohesiveMaterial

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  Properties     matProps = props.getProps ( name );
  Properties     matConf  = conf.makeProps ( name );

  Ref<CohesiveMaterial>  mat;
  String                 type;
  int                    dim;

  matProps.get ( type, "type" );
  matConf .set ( "type", type );

  matProps.get ( dim, "dim"   );
  matConf .set ( "dim", dim   );
  
  if ( type == "XuNeedleman" )
  {
    mat = newInstance<XuNeedlemanMaterial> ( dim, globdat );
  }  
  else if ( type == "TuronMixed" )
  {
    mat = newInstance<TuronCohesiveMaterial> ( dim, globdat );
  }
  else if ( type == "Contact" )
  {
    mat = newInstance<ContactMaterial> ( dim, globdat );
  }
  else
  {
    matProps.propertyError (
      name,
      "invalid cohesive material type: " + type
    );
  }

  return mat;
}
