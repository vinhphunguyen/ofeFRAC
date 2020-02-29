
#include <jive/util/DofSpace.h>


using jive::IntVector;
using jive::util::DofSpace;


//-----------------------------------------------------------------------
//   getDispDofs
//-----------------------------------------------------------------------

// Returns an integer vector with the displacement DOF indices.

IntVector getDispDofs ( const DofSpace& dofs )
{
  const char*  DISP_NAMES[3] = { "dx", "dy", "dz" };

  const int    itemCount = dofs.itemCount ();

  IntVector    iitems    ( itemCount );
  IntVector    jtypes    ( 3 );

  IntVector    idofs;

  int          n = 0;

  // Determine the displacement DOF type indices.

  for ( int i = 0; i < 3; i++ )
  {
    int  j = dofs.findType ( DISP_NAMES[i] );

    if ( j >= 0 )
    {
      jtypes[n++] = j;
    }
  }

  jtypes.reshape ( n );

  for ( int i = 0; i < itemCount; i++ )
  {
    iitems[i] = i;
  }

  idofs.resize ( n * itemCount );

  // To make this function truly general, the displacement DOFs do not
  // have to be defined for all items.

  n = dofs.collectDofIndices ( idofs, iitems, jtypes );

  idofs.reshape ( n );

  return idofs;
}
