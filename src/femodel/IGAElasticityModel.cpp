
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/util/Properties.h>

#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/Database.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/Printer.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/InternalElement.h>
#include <jive/fem/BezierElement.h>
#include <jive/fem/ElementGroup.h>
#include <jive/geom/StdSquare.h>
#include <jive/geom/StdCube.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/StdShapeFactory.h>
#include <jive/geom/StdBezierShape.h>
#include <jive/geom/ParametricArea.h>
#include <jive/geom/ParametricVolume.h>
#include <jive/geom/Geometries.h>
#include <jive/util/Assignable.h>

#include "IGAElasticityModel.h"


using jive::geom::IShape;
using jem::io::FileWriter;
using jive::fem::IElement;
using jem::io::PrintWriter;
using jive::util::Printer;

//=======================================================================
//   class IGAElasticityModel
//=======================================================================

/**
  * This class implements the isogeometric solid mechanics model. 
  * Using Bezier extractor for NURBS/Tsplines.
  *
  * The following Bernstein basis are currently available in my jive
  *
  * Bezier9B:    bi-quadratic Bezier elements
  * Bezier16B:   bi-cubic Bezier elements
  * Bezier25B:   bi-quartic Bezier elements
  * Bezier3x2B:  quadratic-linear Bezier elements
  * Bezier4x2B:  cubic-linear Bezier elements
  * Bezier4x3B:  cubic-quadratic Bezier elements
  * Bezier5x2B:  quartic-linear Bezier elements
  * Bezier5x3B:  quartic-quadratic Bezier elements
  * Bezier6x3B:  quintic-quadratic Bezier elements
  *
  * 3D elements include:
  *
  * Bezier27B:     tri-quadratic Bezier elements
  * Bezier3x2x2B:  quadratic-linear-linear Bezier elements
  * Bezier3x3x2B:  quadratic-quadratic-linear Bezier elements
  *
  * Note: The code assumes the mesh is homogeneous in the sense that there is
  * only one Bernstein basis (p,q). Unstructured Tsplines in which (p,q) varies
  * from element to element is not supported yet. However this can be implemented with ease.
  *
  * 2/3D: verified correctly compared to MIGFEM.
  *
  * This model applies for an element group. There might be more than 1 material in the model.
  * 6 December 2016: using InternalElement class now, so the code works with both Lagrange
  *   elements and NURBS/T-splines elements.
  *
  *
  * VP Nguyen, Cardiff University, Wales, UK.
  * March 2013.
  */

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  IGAElasticityModel::TYPE_NAME     = "IGAElasticity";
const char*  IGAElasticityModel::MATERIAL_PROP = "material";
const char*  IGAElasticityModel::SHAPE_PROP    = "shape";
const char*  IGAElasticityModel::DOF_NAMES[3]  = {"dx","dy","dz"};


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


IGAElasticityModel::IGAElasticityModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Model ( name )

{
  using jive::util::Database;
  using jive::geom::StdShape;
  using jive::geom::StdSquare;
  using jive::geom::StdCube;
  using jive::geom::StdBezierShape;
  using jive::geom::ParametricArea;
  using jive::geom::ParametricVolume;
  using jive::geom::StdShapeFactory;
  using jive::geom::Geometries;
  using jive::util::Assignable;
  using jive::fem::ElementGroup;
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  String         context = getContext ();

  ElementGroup egroup = ElementGroup::get  ( myConf, myProps, globdat, context );
  ElementSet   elems  = egroup.getElements ( );
  NodeSet      nodes  = elems.getNodes     ( );

  elemCount_ = egroup.size( );

  ielems_.resize ( elemCount_ );
  ielems_ = egroup.getIndices ();

  // create dofspace object

  rank_       = nodes.rank ();
  strCount_   = STRAIN_COUNTS[rank_];

  dofs_ = XDofSpace::get ( nodes.getData(), globdat );

  dofTypes_.resize( rank_ );

  for ( int i = 0; i < rank_; i++ )
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i]);
  }

  dofs_->addDofs ( IdxVector ( iarray( nodes.size() ) ), dofTypes_);

  cons_  = Constraints::get ( dofs_, globdat );

  String elemType;
  myProps.get ( elemType , "elemType" );
    
  Ref<IShape>   iShape;

  if ( elemType == "Lagrange" )
  {
    iShape  = IShapeFactory::newInstance(
                 joinNames (myName_, SHAPE_PROP ),
                 conf,
                 props );
    element_ = newInstance<IElement> ( elems, iShape);
  }
  else
  {

    // first create Bernstein shape object 
    // then create the Bezier shape object from the Bernstein object

    Properties     shapeProps = myProps.getProps  ( SHAPE_PROP );
    Properties     shapeConf  = myConf.makeProps  ( SHAPE_PROP );

    String         geom, extDbName;
    IdxVector      gauss;
    
    shapeProps.get  ( gauss,       "gauss"     );

    Ref<StdShape> shape;
    
    if      ( rank_ == 2 )
    {
        shape = StdShapeFactory::newInstance (Geometries::SQUARE, shapeConf, shapeProps );
        iShape = newInstance<ParametricArea> (
                             "shape",
                             StdSquare::getGaussScheme   ( gauss[0], gauss[1] ),
                             newInstance<StdBezierShape> ( shape ) );
    }
    else if ( rank_ == 3 )
    {
        shape = StdShapeFactory::newInstance (Geometries::CUBE, shapeConf, shapeProps );
        iShape = newInstance<ParametricVolume> (
                             "shape",
                             StdCube::getGaussScheme   ( gauss[0], gauss[1], gauss[2] ),
                             newInstance<StdBezierShape> ( shape ) );
    }

    element_ = newInstance<BezierElement> (
      elems,
      Database::get ( "C", elems.getData(), globdat, context ),
      iShape
    );
  }

  // ------------------------------------------------------
  //       create material objects
  // ------------------------------------------------------

  StringVector mats;

  myProps.get ( mats ,MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  int numMats = mats.size ( );
  materials_.  resize  ( numMats     );
  elemMatMap_. resize  ( elemCount_  );

  for ( int iMat = 0; iMat < numMats; iMat++ )
  {
    String matName         = mats[iMat]; 
    Ref<Material> tmpMat   = newMaterial ( matName, myConf, myProps, globdat );
    Properties    matProps = myProps.findProps ( matName );
    Properties    matConf  = myConf.makeProps  ( matName );

    tmpMat->configure ( matProps, globdat );
    tmpMat->getConfig ( matConf,  globdat );

    // get the elements of this material

    ElementGroup  egroup    = ElementGroup::get ( matConf, matProps, globdat, context );
    IdxVector     eIndices  = egroup.getIndices ( ); 

    int elemCount = eIndices.size ();

    for ( int ielem = 0; ielem < elemCount; ielem++ )
    {
      elemMatMap_[eIndices[ielem]] = iMat;
    }
    materials_[iMat] = tmpMat;
  }

  globdat.set ("dd",elemMatMap_);
  
  // correct function to compute B matrix

  getShapeGrads_  = getShapeGradsFunc        ( rank_ );

  //System::out() <<  "IGAElasticityModel:  initialization done.\n";
}


IGAElasticityModel::~IGAElasticityModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool IGAElasticityModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

 /* if ( action == Actions::GET_CONSTRAINTS )
  {
    cons_->printTo ( Printer::get() );
    Printer::get() << "\n\n";
    Printer::flush ();
  }*/

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MBuilder>  mbuilder;

    Vector         state;
    Vector         fint;

    StateVector::get ( state,    dofs_, globdat );
    params      .get ( fint,     ActionParams::INT_VECTOR );
    params      .get ( mbuilder, ActionParams::MATRIX0 );

    calcMatrix_ ( *mbuilder, fint, state );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   calcMatrix_
//-----------------------------------------------------------------------


void IGAElasticityModel::calcMatrix_

  ( MBuilder&      mb,
    const Vector&  fint,
    const Vector&  state ) const

{
  using jem::numeric::matmul;

  IElement&  elem      = *element_;

  const idx_t     nodeCount = elem.nodeCount;
  const idx_t     dofCount  = nodeCount*rank_;
        idx_t     ielem;

  //System::out() << nodeCount << "\n";
  //System::out() << elemCount_ << "\n";

  Matrix          elMat     ( dofCount, dofCount ); // element stiffness matrix
  Ref<Material>   eMat;
  Vector          elVec     ( dofCount );           // element internal force vector
  Vector          elState   ( dofCount );           // element displacement vector
  IdxVector       idofs     ( dofCount );           // element global positions
  idx_t           iMat;
  Matrix          stiff     ( strCount_, strCount_ );
  Vector          stress    ( strCount_ );
  Vector          strain    ( strCount_ );
  Matrix          B         ( strCount_, dofCount  );
  Matrix          Bt        = B.transpose ();

  MChain3        mc3;

  elem.selectPoints ( BezierElement::IPOINTS );

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elem.loadData        ( ielem );
    elem.getShapeGrads   ();
    dofs_->getDofIndices ( idofs, elem.inodes, dofTypes_ );

    // get correct material of current element

    iMat   = elemMatMap_[ielem]; 
    eMat   = materials_[iMat];

    elMat = 0.0;

    for ( idx_t ip = 0; ip < elem.pntCount; ip++ )
    {
      Matrix  grads = elem.shapeGrads[ip]; // get grads of shape functions at ip
      getShapeGrads_ ( B, grads );         // using grads to build B matrix
      eMat->update ( stress, stiff, strain, 0 );
      elMat += elem.pntWeights[ip] * mc3.matmul ( Bt, stiff, B );
    }
    System::out() << elMat << "\n";

    elState = state[idofs];
    elVec   = matmul ( elMat, elState );
    mb.addBlock ( idofs, idofs, elMat );
    fint[idofs] += elVec;
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool IGAElasticityModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat ) const

{
//  using jive::model::ActionParams;
//  using jive::model::StateVector;
//
//  Ref<XTable>  table;
//  Vector       state;
//  String       name;
//
//
//  params.get ( name, ActionParams::TABLE_NAME );
//
//  if ( name == "flux" )
//  {
//    params.get ( table, ActionParams::TABLE );
//
//    if ( table->getRowItems() == element_->elems.getData() )
//    {
//      StateVector::get ( state, dofs_, globdat );
//      calcFlux_        ( *table, state );
//
//      return true;
//    }
//  }
//
  return false;
}


//-----------------------------------------------------------------------
//   calcFlux_
//-----------------------------------------------------------------------
/*

void IGAElasticityModel::calcFlux_

  ( XTable&        table,
    const Vector&  state ) const

{
  using jem::numeric::matmul;

  BezierElement&  elem      = *element_;

  const idx_t     elemCount = elem.elems.size ();

  Vector          elState   ( elem.nodeCount );
  IdxVector       idofs     ( elem.nodeCount );

  IdxMatrix       jcols;
  Matrix          flux;


  elem.selectPoints ( BezierElement::VERTICES );

  jcols.resize ( elem.rank, elem.pntCount );
  flux .resize ( elem.rank, elem.pntCount );

  for ( idx_t j = 0; j < elem.pntCount; j++ )
  {
    for ( idx_t i = 0; i < elem.rank; i++ )
    {
      jcols(i,j) = table.addColumn (
        String::format ( "f(%d,%d)", i, j )
      );
    }
  }

  for ( idx_t ielem = 0; ielem < elemCount; ielem++ )
  {
    elem.loadData        ( ielem );
    elem.getShapeGrads   ();
    dofs_->getDofIndices ( idofs, elem.inodes, dofTypes_ );

    elState = state[idofs];

    for ( idx_t i = 0; i < elem.pntCount; i++ )
    {
      flux[i] = matmul ( elem.shapeGrads[i], elState );
    }

    table.setRowValues ( ielem, flatten( jcols ), flatten( flux ) );
  }
}
*/


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newESolidModel
//-----------------------------------------------------------------------


static Ref<Model>     newIGAElasticityModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<IGAElasticityModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareIGAElasticityModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( IGAElasticityModel::TYPE_NAME, & newIGAElasticityModel );
}

