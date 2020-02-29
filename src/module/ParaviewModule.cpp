/**
 *  This class is used to write FE results to VTU files that can be processed by Paraview
 *  It writes one PVD file that contains a collection of VTU files of which each
 *  VTU file is associated with one time step.
 *  Note that the FE mesh is repeated in every VTU files.
 *  When using this module, you should shutdown the simulation by typing "exit"
 *  (not by Ctrl-C) otherwise the PVD file is not created.
 *  This pvd is now kind of obsolete as Paraview automatically groups vtu files together.
 *
 *
 *  VP Nguyen, Cardiff University, UK. March 2013.
 */

#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/base/limits.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/DenseTable.h>
#include <jive/util/Printer.h>
#include <jive/model/StateVector.h>
#include <jive/model/Model.h>
#include <jive/model/Actions.h>
#include <jive/implict/SolverInfo.h>
#include <jive/app/ModuleFactory.h>
#include <jive/SparseMatrix.h>

#include "ParaviewModule.h"
#include "util/utilities.h"

using namespace jem;
using jem::io::FileWriter;
using jive::Matrix;
using jive::IdxVector;
using jive::Vector;
using jem::io::endl;
using jive::util::XDofSpace;
using jive::model::StateVector;
using jive::util::Globdat;
using jive::util::Printer;
using jive::util::DenseTable;
using jive::SparseMatrix;

//=======================================================================
//   class ParaviewModule
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  ParaviewModule::FILE_PROP     = "fileName";
const char*  ParaviewModule::DATA_PROP     = "data";
const char*  ParaviewModule::INTERVAL_PROP = "interval";
const char*  ParaviewModule::DOF_PROP      = "dofs";
const char*  ParaviewModule::TABLE_PROP    = "table";



//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


ParaviewModule::ParaviewModule ( const String& name ) :
  Super ( name )
{}

ParaviewModule::~ParaviewModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status ParaviewModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
   using jive::mp::Globdat;

   Properties  myProps = props.getProps ( myName_ );
   Properties  myConf  = conf.makeProps ( myName_ );

   const String context = getContext();

   egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );
   elems_  = egroup_.getElements ( );
   nodes_  = elems_.getNodes     ( );
   dofs_   = XDofSpace::get ( nodes_.getData(), globdat );
  
   nodeCount_    = nodes_.size     ( );
   elemCount_    = egroup_.size    ( );
   rank_         = nodes_.rank     ( );
   
   // all dof types in the mesh
   // for a hydromechanical model with fracturing fluid we have
   // 2D: dx, dy, dp, dpf
   // 3D: dx,dy,dz, dp, dpf
   dofTypeCount_ = dofs_->typeCount( ); // all dof types in the mesh
  
   StringVector typeNames; // = dofs_->getTypeNames ();
   myProps.get  ( typeNames,   DOF_PROP );
   
   dofTypes_.resize ( typeNames.size() );

   for (idx_t i = 0; i < typeNames.size(); i++ )
   {
     dofTypes_[i] = dofs_->getTypeIndex ( typeNames[i] );
   }

   ielems_.resize ( elemCount_ );
   ielems_ = egroup_.getIndices ();
   noNodePerElem_= elems_.getElemNodeCount(ielems_[0]);

   mpx_  = Globdat::getMPContext ( globdat );

   // default values for Q4, Q8 and Q9
   numVertexesPerCell_ = 4;
   vtkCellCode_        = 9;

   // triangle elements T3 and T^6
   if ( (rank_ == 2) & (noNodePerElem_%3 ==0) )
   {
     numVertexesPerCell_ = 3;
     vtkCellCode_        = 5;
   }

   // 3D elements: only 8-node hexahedron elements
   // and T4 elements
   if ( (rank_ == 3) & (noNodePerElem_ ==8) )
   {
     numVertexesPerCell_ = 8;
     vtkCellCode_        = 12;
   }
   if ( (rank_ == 3) & (noNodePerElem_ ==4) )
   {
     numVertexesPerCell_ = 4;
     vtkCellCode_        = 10;
   }

   isTable_ = false;

   return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

Module::Status ParaviewModule::run ( const Properties& globdat )
{
  using jem::Ref;
  using jem::System;
  using jive::Vector;
  using jive::util::ItemSet;
  using jive::model::Model;
  using jive::model::ActionParams;
  using jive::model::Actions;

  const String   context = getContext ();
        int      step;

  globdat.get ( step, Globdat::TIME_STEP );

  if ( ( step % interval_ ) != 0 ) return OK;

  Ref<PrintWriter> vtuFile   = newInstance<PrintWriter> ( 
                                   newInstance<FileWriter> 
		                    ( fileName_ + String(step)  + ".vtu" ) );

  //if ( mpx_->myRank () == 0 )
  //{
  writeMeshVTKFile_   ( vtuFile );
  writeStateVTKFile_  ( vtuFile, globdat, context );

  // write data retrieved from models
  // Actually this module is asking the model (all models defined in the input
  // file) to compute the tables that are corresponding with the dataNames_.

  const idx_t dataCount = dataNames_.size();

  Properties      params  ("actionParams");
  Vector          weights ( nodeCount_ );

  for ( idx_t i = 0; i < dataCount; i++ )
  {
    String          tname (dataNames_[i]);

    weights = 0.0;

    Ref<Model>      model  = Model::get ( globdat, getContext() );
    Ref<XTable>     xtable = newInstance<DenseTable>  ( tname, nodes_.getData() );

    params.set ( ActionParams::TABLE,         xtable  );
    params.set ( ActionParams::TABLE_NAME,    tname   );
    params.set ( ActionParams::TABLE_WEIGHTS, weights );

    bool updated = model->takeAction ( Actions::GET_TABLE, params, globdat );

    weights = where ( abs( weights ) < Limits<double>::TINY_VALUE, 1.0, 1.0 / weights );

    xtable->scaleRows    ( weights );
    
    if (updated) 
    {
      writeTableVTKFile_   ( vtuFile, xtable, tname );
    }
  }
  
  if ( isTable_ )
  {
    System::out() << "OK OK OK\n";
    Ref<XTable>  phi = XTable::get ( tabName_, nodes_.getData(), globdat );
   phi->printTo ( Printer::get() );
    writeTableVTKFile_   ( vtuFile, phi, "phi" );
  }

  writeClosingVTKFile_   ( vtuFile );

  vtuFile->flush();
  vtuFile->close();
 // }

  return OK;
}

//-----------------------------------------------------------------------
//   writeMeshVTKFile_
//-----------------------------------------------------------------------
    

void ParaviewModule::writeMeshVTKFile_ ( Ref<PrintWriter> vtuFile )
{
  const String str1 = String::format("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"> \n"
		                             "<UnstructuredGrid> \n"
	                                 "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\"> \n", 
				                     nodeCount_, elemCount_ );

  *vtuFile << str1;

  // get node coordinates

  *vtuFile << "<Points> \n";
  *vtuFile << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" >\n";

  Matrix coords(rank_,nodeCount_);

  nodes_.getCoords ( coords );

  if ( rank_ == 2 )
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuFile << coords(0,i) << " " << coords(1,i) << " " << 0. << '\n';
    }
  }
  else
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuFile << coords(0,i) << " " << coords(1,i) << " " << coords(2,i) << '\n';
    }
  }

  *vtuFile << "</DataArray>\n </Points>\n";

  // write element connectivity

  IdxVector inodes(noNodePerElem_);

  *vtuFile << "<Cells> \n";
  *vtuFile << "<DataArray  type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  { 
    int  ielem = ielems_[ie]; 
    elems_.getElemNodes  ( inodes, ielem );

    // differentiate between first-order and high-order elements
    // for high order elements, only write the corner nodes which are
    // sufficient for visualisation.

    if (noNodePerElem_ == numVertexesPerCell_ )
    {
      for ( idx_t in = 0; in < noNodePerElem_; in++ ){
        *vtuFile << inodes[in] << " "; 
      }
    }
    else
    {
      for ( idx_t in = 0; in < numVertexesPerCell_; in++ ){
        *vtuFile << inodes[in*2] << " "; 
      }
    }

    *vtuFile << endl;
  }

  *vtuFile << "</DataArray>\n";

  // write cell offset

  *vtuFile << "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n";

  int offset = 0;
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    offset += numVertexesPerCell_;
    *vtuFile <<  offset << endl;
  }

  *vtuFile << "</DataArray>\n";

  // Print cell types
  
  *vtuFile << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\"> \n";
  
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
     *vtuFile <<  vtkCellCode_ << endl;
  }

  *vtuFile << "</DataArray> \n </Cells> \n";

}

//-----------------------------------------------------------------------
//   writeStateVTKFile_
//-----------------------------------------------------------------------
/*
 * Write state (displacements, velocities, accelerations, temperatures..) to VTK file. 
 */
    

void ParaviewModule::writeStateVTKFile_ ( Ref<PrintWriter> vtuFile,
                                           const Properties& globdat,
                                           const String&     context )
{
   *vtuFile << "<PointData  Vectors=\"sigma\"> \n";
   
   // print displacement field

   Vector  state;
   StateVector::get ( state, dofs_, globdat );

   int maxDofCount = rank_ == 2 ? rank_+1 : rank_;

   const String str1 = String::format("<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , maxDofCount);

   //*vtuFile << " <DataArray  type=\"Float64\"  Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\"> \n";
   *vtuFile << str1;

   IdxVector  idofs    ( rank_       );
   Vector     nodeDisp ( maxDofCount );

   for (idx_t i = 0; i < nodeCount_; ++i)
   {
      nodeDisp = 0.;
      try{
        dofs_->getDofIndices ( idofs, i, dofTypes_[slice(BEGIN, rank_)] );

        nodeDisp[slice(BEGIN,rank_)] = select ( state, idofs );
      }
      catch ( const jem::IllegalInputException& ex )
      {
      }

      for (idx_t j = 0; j< maxDofCount; ++j)
      {
         *vtuFile << String::format( "%12.6e   ", nodeDisp[j] );
      }
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";

   if ( dofTypeCount_ > rank_ )
   {
   const String str2 = String::format("<DataArray type=\"Float64\" Name=\"scalar\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , 1 );

   *vtuFile << str2;

   idx_t    idof;
   double   val;

   for (idx_t i = 0; i < nodeCount_; ++i)
   {
      val = 0.;
      try{
        idof = dofs_->getDofIndex ( i, dofTypes_[rank_] );
        val  = state[idof];
      }
      catch ( const jem::IllegalInputException& ex )
      {
      }
      *vtuFile << String::format( "%12.6e   ", val );
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";
   }
}

//-----------------------------------------------------------------------
//   writeClosingVTKFile_
//-----------------------------------------------------------------------
    
void ParaviewModule::writeClosingVTKFile_ ( Ref<PrintWriter> vtuFile )
{
   *vtuFile << "</PointData> \n";
   // end of VTK file

   *vtuFile << "</Piece> \n";
   *vtuFile << "</UnstructuredGrid> \n";
   *vtuFile << "</VTKFile> \n";
}

//-----------------------------------------------------------------------
//   writeTableVTKFile_
//-----------------------------------------------------------------------
    

void ParaviewModule::writeTableVTKFile_ 

     ( Ref<PrintWriter> vtuFile,
       Ref<XTable>      table,
       const String&    tname )
{
   using jive::StringVector;

   StringVector     colNames = table->getColumnNames  ();
   idx_t            colCount = colNames.size          ();
   idx_t            rowCount = table->rowCount        ();

   Matrix           coord; // matrix representation of the table

   //table->printTo ( Printer::get() );

   // convert a table to a matrix (function in utilities.cpp)
   
   matrixFromTable ( coord, *table, colNames );

   //SparseMatrix coord = table->toMatrix();
   
   *vtuFile << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n",
                             tname, colCount);

   for (idx_t i = 0; i < rowCount; ++i)
   {
      for (idx_t j = 0; j< colCount; ++j)
      *vtuFile << String::format( "%12.8f   ", coord(i,j) );
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void ParaviewModule::shutdown ( const Properties& globdat )
{
   int     step;

   globdat.get ( step, Globdat::TIME_STEP );

   Ref<PrintWriter> pvdFile   = newInstance<PrintWriter> (newInstance<FileWriter> 
		                    ( fileName_ + ".pvd" ));
   *pvdFile << "<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n";
   *pvdFile << "<Collection>\n";

   for (idx_t i = 1; i <= step; ++i)
   {
      if ( ( i % interval_ ) != 0 ) continue;
      
      String fileName = fileName_ + String(i)  + ".vtu";
      *pvdFile << "<DataSet file='"+ fileName + "' groups='' part='0' timestep='"+String(i)+"'/>\n";
   }

   *pvdFile << "</Collection>\n";
   *pvdFile << "</VTKFile>\n";
  
   pvdFile->flush();
   pvdFile->close();

  //writePVDFile_       ( step );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ParaviewModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  interval_ = 1;

  Properties  myProps = props.findProps ( myName_ );

  myProps.get  ( fileName_,     FILE_PROP     );
  myProps.get  ( dataNames_,    DATA_PROP     );
  myProps.find ( interval_,     INTERVAL_PROP );
  
  if ( myProps.find ( tabName_, TABLE_PROP ) )
  {
    isTable_ = true;
  }
    
  //pvdFile_    = newInstance<PrintWriter> (newInstance<FileWriter> ( fileName_+".pvd" ));
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ParaviewModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );
  
  myConf.set  ( FILE_PROP, fileName_     );
  myConf.set  ( DATA_PROP, dataNames_    );
  myConf.set  ( INTERVAL_PROP, interval_ );
  myConf.set  ( TABLE_PROP, isTable_     );
}

Ref<Module>  ParaviewModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )
{
  return newInstance<Self> ( name );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareParaviewModule
//-----------------------------------------------------------------------

void declareParaviewModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "Paraview", & ParaviewModule::makeNew );
}
