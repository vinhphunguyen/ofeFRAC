
#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Properties.h>
#include <jive/util/Assignable.h>
#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/geom/Line.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/ElementGroup.h>


using namespace jem;

using jem::util::Properties;
using jive::Vector;
using jive::StringVector;
using jive::Matrix;
using jive::IdxVector;
using jive::util::Assignable;
using jive::util::XDofSpace;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::geom::BShape;

//=======================================================================
//   typedefs
//=======================================================================

typedef ElementSet    ElemSet;
typedef ElementGroup  ElemGroup;

//=======================================================================
//   class LineLoadModel
//=======================================================================

/*
 * This class implements a line load model to compute the external force 
 * vector due to traction (N/m) perpendicular to the boundary on which it is
 * applied.
 */

class LineLoadModel : public Model
{
 public:

  typedef Model             Super;
  typedef LineLoadModel     Self;

  static const char*        TYPE_NAME;
  static const char*        LOAD_PROP;
  static const char*        SHAPE_PROP;


                            LineLoadModel ();

                            LineLoadModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~LineLoadModel ();


 private:

  void                      evalLoad_

    ( const Vector&           vec,
      double                  scale,
      const Properties&       globdat )      const;


 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;
  Ref<BShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;

  double                    load_;
};


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LineLoadModel::TYPE_NAME  = "LineLoad";
const char*  LineLoadModel::LOAD_PROP  = "load";
const char*  LineLoadModel::SHAPE_PROP = "shape";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


LineLoadModel::LineLoadModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::geom::BShapeFactory;

  Properties    myConf  = conf .makeProps ( myName_ );
  Properties    myProps = props.findProps ( myName_ );

  const String  context = getContext ();

  String        shapeType;


  egroup_ = ElemGroup::get     ( myConf, myProps, globdat, context );
  elems_  = egroup_.getElements ();
  nodes_  = elems_.getNodes    ();

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  const int rank = nodes_.rank ();

  dofTypes_.resize( rank );
  for ( int i = 0; i < rank; i++ )
  {
    dofTypes_[i] = i ;
  }

  shape_ = BShapeFactory::newInstance
    (joinNames (myName_, SHAPE_PROP ), conf, props);
  
  myProps.get ( load_, LOAD_PROP );

  elems_.checkSomeElements ( context,
                             egroup_.getIndices(),
			     shape_->nodeCount() );
}


LineLoadModel::~LineLoadModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LineLoadModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  if ( ! props.contains( myName_ ) )
  {
    return;
  }

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LineLoadModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( LOAD_PROP, load_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LineLoadModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jem::Exception;
  using jive::model::Actions;
  using jive::model::ActionParams;


  if ( action == Actions::GET_EXT_VECTOR )
  {
    Vector  fext;
    double  scale;

    scale = 1.0;

    params.find ( scale, ActionParams::SCALE_FACTOR );
    params.get  ( fext,  ActionParams::EXT_VECTOR   );
    evalLoad_   ( fext,  scale, globdat );

    //System::out() << fext << "\n";

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   evalLoad_
//-----------------------------------------------------------------------


void LineLoadModel::evalLoad_

  ( const Vector&      fext,
    double             scale,
    const Properties&  globdat ) const

{
  using jem::select;
  using jem::numeric::norm2;

  IdxVector    ielems    = egroup_.getIndices ();
  Matrix       sfuncs    = shape_->getShapeFunctions ();

  const idx_t  rank      = nodes_.rank ();
  const idx_t  elemCount = ielems.size ();
  const idx_t  nodeCount = shape_->nodeCount   ();
  const idx_t  ipCount   = shape_->ipointCount ();
  const idx_t  loadCount = dofTypes_.size ();
  const idx_t  dofCount  = nodeCount * loadCount;

  Vector       fe        ( dofCount );
  IdxVector    idofs     ( dofCount );
  IdxVector    inodes    ( nodeCount );
  Vector       weights   ( ipCount );
  Matrix       coords    ( rank, nodeCount );
  Matrix       normals   ( rank, ipCount   );

  // Loop over the boundary elements (can be line elements or surface
  // elements)

  for ( idx_t ie = 0; ie < elemCount; ie++ )
  {
    int  ielem = ielems[ie];

    elems_.getElemNodes   ( inodes, ielem );
    nodes_.getSomeCoords  ( coords, inodes );
    dofs_->getDofIndices  ( idofs,  inodes, dofTypes_ );

    shape_->getNormals ( normals, weights, coords );
    
    //System::out() << coords << "\n";
    //System::out() << ie << "\n";


    // fe = int_{transpose(N) * traction}
    // fe stored in column vector = {f1x, f1y, f2x, f2y, ..., fnx, fny}^T (2D)
    // loop on integration points
    //  loop on directions: x, y and z

    fe = 0.0;

    for ( idx_t ip = 0; ip < ipCount; ip++ )
    {
      for ( idx_t j = 0; j < loadCount; j++ )
      {
        double  s = scale * weights[ip] * load_ * normals(j,ip);

        //System::out () << normals(j,ip) << "\n";

	    fe[slice(j,END,loadCount)] += s * sfuncs(ALL,ip);
      }
    }

    // assemble fe into global external force vector

    select ( fext, idofs ) += fe;
  }

  //System::out() << "fext: " << fext << "\n";
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newLineLoadModel
//-----------------------------------------------------------------------


Ref<Model>            newLineLoadModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<LineLoadModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareLineLoadModel
//-----------------------------------------------------------------------


void declareLineLoadModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( LineLoadModel::TYPE_NAME, newLineLoadModel );
}
