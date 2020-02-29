
#include <jem/util/Properties.h>

#include "SBDiffusivity.h"
#include "ConstantDiffusivity.h"

using namespace jem;


//-----------------------------------------------------------------------
//   newDiffusitivity
//-----------------------------------------------------------------------


Ref<Diffusivity>         newDiffusivity

  ( const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  Ref<Diffusivity>    mat;
  String              type;

  props.get ( type, "type" );
  conf .set ( "type", type );

  if            ( type == "StressBased" )
  {
    mat = newInstance<SBDiffusivity>   ( conf, props, globdat );
  }
  else if       ( type == "Constant" )
  {
    mat = newInstance<ConstantDiffusivity>   ( conf, props, globdat );
  }
  else
  {
    props.propertyError (
      props.getName(),
      "invalid diffusitivity type: " + type
    );
  }

  return mat;
}
