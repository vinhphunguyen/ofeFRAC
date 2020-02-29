
#include "Material.h"

Material::Material 

  ( const int          rank,
    const Properties&  globdat )

  : rank_ ( rank )

{}

Material::~Material()
{}

void   Material::commit()
{}


void   Material::configure ( const Properties& props,
                             const Properties& globdat )
{}

void   Material::getConfig ( const Properties& props,
                             const Properties& globdat ) const
{}

void Material::allocPoints ( int count )
{}

double Material::checkLocalisation (       Vector& normal,
                                    const Vector& stress,
                                    const Matrix& tangent,
                                         int     ip ) const
{ return 0.;}

void  Material::            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he )
{}
  
void  Material::            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint, 
      double                he )
{}

Vector Material::           giveStress

    ( int                    ip ) const
{
    Vector v(3);
    return v;
}

