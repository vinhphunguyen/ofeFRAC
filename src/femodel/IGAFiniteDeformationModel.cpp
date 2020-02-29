
#include <jem/io/FileWriter.h>

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/util/Properties.h>
#include <jive/util/XTable.h>
#include <jive/util/Database.h>
#include <jive/util/XDofSpace.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/BezierElement.h>
#include <jive/fem/ElementGroup.h>
#include <jive/geom/StdSquare.h>
#include <jive/geom/StdCube.h>
#include <jive/geom/StdShapeFactory.h>
#include <jive/geom/StdBezierShape.h>
#include <jive/geom/ParametricArea.h>
#include <jive/geom/ParametricVolume.h>
#include <jive/geom/Geometries.h>
#include <jive/util/Assignable.h>

#include "IGAFiniteDeformationModel.h"


//=======================================================================
//   class IGAFiniteDeformationModel
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
  *
  *
  * VP Nguyen, Cardiff University, Wales, UK.
  * March 2013.
  */

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  IGAFiniteDeformationModel::TYPE_NAME     = "IGAFiniteDeformation";
const char*  IGAFiniteDeformationModel::MATERIAL_PROP = "material";
const char*  IGAFiniteDeformationModel::SHAPE_PROP    = "shape";
const char*  IGAFiniteDeformationModel::DOF_NAMES[3]  = {"ux","uy","uz"};


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


IGAFiniteDeformationModel::IGAFiniteDeformationModel

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

  // first create Bernstein shape object 
  // then create the Bezier shape object from the Bernstein object

  Properties     shapeProps = myProps.getProps  ( SHAPE_PROP );
  Properties     shapeConf  = myConf.makeProps  ( SHAPE_PROP );

  String         geom, extDbName;
  IdxVector      gauss;
  
  shapeProps.get  ( gauss,       "gauss"     );
  shapeProps.get  ( extDbName,   "extDbName" );

  Ref<StdShape> shape;
  Ref<IShape>   iShape;
  
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
    Database::get ( extDbName, elems.getData(), globdat, context ),
    iShape
  );

  // ------------------------------------------------------
  //       create material objects
  // ------------------------------------------------------

  StringVector mats;

  myProps.get ( mats ,MATERIAL_PROP );
  myConf. set ( MATERIAL_PROP, mats );

  int numMats = mats.size ( );
  materials_. resize  ( numMats  );

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
  
  // correct function to compute B matrix

  getShapeGrads_     = getShapeGradsFunc ( rank_ );
  getBMatrixLinFunc_ = getBMatrixLinFunc ( rank_ );
  
  //System::out() <<  "IGAElasticityModel:  initialization done.\n";
}


IGAFiniteDeformationModel::~IGAFiniteDeformationModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool IGAFiniteDeformationModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;


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


void IGAFiniteDeformationModel::calcMatrix_

  ( MBuilder&      mb,
    const Vector&  fint,
    const Vector&  state ) const

{
  using jem::numeric::matmul;

  BezierElement&  elem      = *element_;

  const idx_t     nodeCount = elem.nodeCount;
  const idx_t     dofCount  = nodeCount*rank_;
        idx_t     ielem;

  //System::out() << nodeCount << "\n";
  //System::out() << elemCount_ << "\n";

  Matrix          elMat     ( dofCount, dofCount ); // element stiffness matrix
  Ref<Material>   eMat;
  Vector          elVec     ( dofCount );           // element internal force vector
  Vector          elDisp    ( dofCount );           // element displacement vector
  IdxVector       idofs     ( dofCount );           // element global positions
  idx_t           iMat;
  Matrix          stiff     ( strCount_, strCount_ );
  Vector          stress    ( strCount_ );
  Vector          strain    ( strCount_ );
  Matrix          B         ( strCount_, dofCount  );
  Matrix          Bt        = B.transpose ();
  Matrix          F          ( rank_, rank_ ); // deformation gradient matrix
  
  MChain3        mc3;
  MChain1        mc1;

  elem.selectPoints ( BezierElement::IPOINTS );

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    ielem = ielems_[ie];
    elem.loadData        ( ielem );
    elem.getShapeGrads   ();
    dofs_->getDofIndices ( idofs, elem.inodes, dofTypes_ );

    // get correct material of current element

    iMat    = elemMatMap_[ielem]; 
    eMat    = materials_[iMat];
    elDisp  = state[idofs];

    elMat = 0.0;
    elVec = 0.0;

    for ( idx_t ip = 0; ip < elem.pntCount; ip++ )
    {
      Matrix  grads = elem.shapeGrads[ip]; // get grads of shape functions at ip

      evalDeformationGradient ( F, elDisp, grads);
      getGreenLagrangeStrain  ( strain, F   );
      getBMatrixLinFunc_      ( B, F, grads );
      eMat->update ( stress, stiff, strain, 0 );
      //System::out()<<B<<"\n";
      //System::out()<<stress<<"\n";
      elMat += elem.pntWeights[ip] * mc3.matmul ( Bt, stiff, B );
      elVec += elem.pntWeights[ip] * mc1.matmul ( Bt, stress   );

      addElemMatLargeD ( elMat, stress, grads,  elem.pntWeights[ip]  );
    }

    mb.addBlock ( idofs, idofs, elMat );
    fint[idofs] += elVec;
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool IGAFiniteDeformationModel::getTable_

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


//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newESolidModel
//-----------------------------------------------------------------------


static Ref<Model>     newIGAFiniteDeformationModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<IGAFiniteDeformationModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSolidModel
//-----------------------------------------------------------------------


void declareIGAFiniteDeformationModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( IGAFiniteDeformationModel::TYPE_NAME, & newIGAFiniteDeformationModel );
}

