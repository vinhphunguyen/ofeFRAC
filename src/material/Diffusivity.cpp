
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/util/Properties.h>

#include "Diffusivity.h"

using namespace jem;
using jem::util::Properties;

// ------------------------------------------------
//  Constructor
// ------------------------------------------------

Diffusivity::Diffusivity 

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  props.get ( c0_, "c0" );
  conf .set ( "c0", c0_ );
}

// ------------------------------------------------
//  Destructor
// ------------------------------------------------

Diffusivity::~Diffusivity()
{}
