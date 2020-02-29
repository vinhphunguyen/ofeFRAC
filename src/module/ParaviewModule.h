/**
 *  This class is used to write FE results to VTU files that can be processed by Paraview
 *  It writes one PVD file that contains a collection of VTU files of which each
 *  VTU file is associated with one time step.
 * 
 *  VP Nguyen, Cardiff University, UK. March 2013.
 */

#ifndef PARAVIEW_MODULE_H
#define PARAVIEW_MODULE_H

#include <jive/app/Module.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Assignable.h>
#include <jem/mp/utilities.h>
#include <jive/mp/Globdat.h>



namespace jem
{
  namespace io
  {
    class PrintWriter;	  
  }
}

namespace jive
{
  namespace util
  {
    class XTable;	  
    class XDofSpace;
  }
}

using namespace jem;
using jem::String;
using jem::io::PrintWriter;
using jem::util::Properties;
using jem::mp::Context;

using jive::app::Module;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::util::Assignable;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::IdxVector;
using jive::StringVector;

typedef ElementSet              ElemSet;
typedef ElementGroup            ElemGroup;


//-----------------------------------------------------------------------
//   class ParaviewModule
//-----------------------------------------------------------------------


class ParaviewModule : public Module
{
 public:

  typedef ParaviewModule    Self;
  typedef Module            Super;

  static const char*        FILE_PROP;
  static const char*        DATA_PROP;
  static const char*        INTERVAL_PROP;
  static const char*        DOF_PROP;
  static const char*        TABLE_PROP;


  explicit                  ParaviewModule

    ( const String&           name = "paraview" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )        const;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       globdat,
      const Properties&       props );

 protected:

  virtual                  ~ParaviewModule   ();

 private:

  void                     writeMeshVTKFile_    ( Ref<PrintWriter> vtuFile );
  void                     writeClosingVTKFile_ ( Ref<PrintWriter> vtuFile );

  void                     writeStateVTKFile_  
    ( Ref<PrintWriter>       vtuFile,
      const Properties&      globdat,
      const String&          context );

  void                     writeTableVTKFile_    
     
   ( Ref<PrintWriter> vtuFile,
     Ref<XTable>      table,
     const String&    tname );

 private:

  String                   fileName_;
  Ref<PrintWriter>         pvdFile_;

  Assignable<ElemGroup>    egroup_;
  Assignable<ElemSet>      elems_;
  Assignable<NodeSet>      nodes_;
  Ref<XDofSpace>           dofs_;
  Ref<Context>             mpx_;

  idx_t                    dofTypeCount_;
  idx_t                    numVertexesPerCell_;
  idx_t                    vtkCellCode_;
  IdxVector                dofTypes_;

  int                      elemCount_;
  int                      nodeCount_;
  int                      noNodePerElem_;
  int                      rank_;
  int                      interval_; // write only for intervals

  IdxVector                inodes_;
  IdxVector                ielems_;
  StringVector             dataNames_;

  String                   tabName_; // write Table already stored in globdat
  bool                     isTable_; // e.g., from Staggered Phase field calculations
};


#endif
