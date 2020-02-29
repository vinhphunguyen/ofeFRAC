#ifndef IMPORT_H
#define IMPORT_H

#include <jive/Array.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/matmul.h>

#include "util/utilities.h"
#include "material/DamageMaterial.h"
#include "material/Material.h"

namespace jem
{
  namespace util
  {
    class Properties;
  }
}

namespace jive
{
  namespace util
  {
    class XTable;
    class XDofSpace;
  }

  namespace algebra
  {
    class   MatrixBuilder;

    typedef MatrixBuilder  MBuilder;
  }

  namespace model
  {
    class Model;
  }

  namespace fem
  {
    class NodeSet;
    class ElementSet;
    class BezierElement;
    class InternalElement;
    class ElementGroup;
  }

  namespace geom
  {
    class   InternalShape;
    typedef InternalShape IShape;
  }
}

using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::numeric::Function;
using jem::util::Flex;
using jem::util::Properties;
using jive::Vector;
using jive::IdxVector;
using jive::Matrix;
using jive::Cubix;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::algebra::MBuilder;
using jive::model::Model;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::fem::BezierElement;

typedef MatmulChain<double,3>   MChain3;
typedef MatmulChain<double,1>   MChain1;


#endif
