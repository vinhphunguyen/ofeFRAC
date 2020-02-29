/*
 *  Copyright (C) 2012 TU Delft. All rights reserved.
 *  
 *  Object that gathers data in a table for output purposes.
 *  All models compute all data.
 *
 *  Author: F.P. van der Meer
 *  Date: November 2012
 *
 */

#include <jem/base/array/select.h>
#include <jem/base/array/operators.h>
#include <jem/util/StringUtils.h>

#include "utilities.h"
#include "TbFiller.h"

using jem::io::endl;
using jem::util::StringUtils;

typedef  ArrayBuffer<idx_t>    IdxArrayB;

//=======================================================================
//   class TbFiller::NamedSlice
//=======================================================================

TbFiller::NamedSlice::NamedSlice () 
{}

TbFiller::NamedSlice::NamedSlice

  ( const Slice   slice,
    const String& name ) : slice_(slice), name_(name)
{}

Slice TbFiller::NamedSlice::slice() const
{
  return slice_;
}

String TbFiller::NamedSlice::name() const
{
  return name_;
}

//=======================================================================
//   class TbFiller
//=======================================================================

// ------------------------------------------------------------------------
//   static members
// ------------------------------------------------------------------------

IdxVector    TbFiller::permi_         = IdxVector();
IdxVector    TbFiller::permj1_        = IdxVector();
IdxVector    TbFiller::permj2_        = IdxVector();

// const char*  TbFiller::X_NAMES[3]     = {"x","y","z"};

//-----------------------------------------------------------------------
//   constructor,destructor
//-----------------------------------------------------------------------

TbFiller::TbFiller 

    ( const idx_t       rank ) : rank_ ( rank )

{
  strCount_ = STRAIN_COUNTS[rank_];
  ntype_    = 0;
}

TbFiller::~TbFiller ()
{}

//-----------------------------------------------------------------------
//   announce
//-----------------------------------------------------------------------

Slice TbFiller::announce
  
  ( const String&  fullName )

{
  // typical usage:
  // Slice      iistrain   = tbFiller.announce ( "strain.tensor" );
  // Slice      iistress   = tbFiller.announce ( "stress.tensor" );
  // Slice      iidamage   = tbFiller.announce ( "damage"        );

  StringVector parts  = StringUtils::split ( fullName, '.' );
  idx_t        ncomp  = 0;
  idx_t        before = ntype_;
  String       name;

  if ( parts.size() > 1 )
  {
           name = parts[0];  // -> 'strain'
    String type = parts[1];  // -> 'tensor'

    if      ( type == "tensor" )
    {
      ncomp = strCount_;
      announceTensor_ ( name );
    }
    else if ( type == "vector" )
    {
      ncomp = rank_;
      announceVector_ ( name );
    }
    else if ( type == "diag" )
    {
      ncomp = 3;
      announceDiag_ ( name );
    }
    else
    {
      System::warn() << "unknonwn data type in " << fullName << endl;
      ncomp = 1;
    }
  }
  else
  {
    typeNames_.pushBack ( fullName );
    ncomp = 1;
  }

  ntype_ += ncomp;

  SliceFromTo  sl  ( before, ntype_ );

  if ( parts.size() > 1 )
  {
    superNames_.pushBack ( NamedSlice(sl,name) );
  }

  return sl;
}

//-----------------------------------------------------------------------
//   announce
//-----------------------------------------------------------------------

Slice TbFiller::announce
  
  ( const StringVector&  names )

{
  idx_t before = ntype_;

  for ( idx_t i = 0; i < names.size(); ++ i )
  {
    announce ( names[i] );
  }

  return SliceFromTo ( before, ntype_ );
}

//-----------------------------------------------------------------------
//   setFilter
//-----------------------------------------------------------------------

void TbFiller::setFilter
  
  ( const String&        filter )

{
  // decide which data is to be written to file
  // xoutput = 
  // {
  //   table = " stress | eqveps | strain_xy"; -> 'filter'
  // };

  StringVector strs = StringUtils::split ( filter, '|' );

  write_ . resize ( ntype_ );
  write_ = false;

  // loop over 'table' -> stress, eqveps, strain_xy
   
  for ( idx_t istr = 0; istr < strs.size(); ++istr )
  {
    String str = strs[istr].stripWhite(); //System::out() << typeNames_ <<"\n";
    bool found = false;

    // first try whether it refers to an defined single type

    for ( idx_t itype = 0; itype < ntype_; ++itype )
    {
      if ( str.equals ( typeNames_[itype] ) )
      {
        write_[itype] = found = true;
        break;
      }
    }

    // then try shorthand for full stress/strain vectors etc

    for ( idx_t isup = 0; !found && isup < superNames_.size(); ++isup )
    {
      if ( str.equals ( superNames_[isup].name() ) )
      {
        write_[superNames_[isup].slice()] = found = true;
      }
    }

    if ( !found )
    {
//       // debugging information, printed for every model: too annoying
//       System::out() << "Specified output `" << str << "` is unknown. It " <<
//         "will be ignored.\n" << 
//         "Possible candidates are: \n" << typeNames_ << endl;
    }
  }
}

//-----------------------------------------------------------------------
//   prepareTable
//-----------------------------------------------------------------------

void TbFiller::prepareTable
  
  (       IdxVector&  i2table,
          IdxVector&  jcols,
    const Ref<XTable> table ) const

{
  ArrayBuffer<idx_t> ibuf;
  ArrayBuffer<idx_t> jbuf;

  for ( idx_t itype = 0; itype < ntype_; ++itype )
  {
    if ( write_[itype] )
    {
      ibuf.pushBack ( itype );
      jbuf.pushBack ( table->addColumn ( typeNames_[itype] ) );
    }
  }

  i2table.ref ( ibuf.toArray() );
  jcols  .ref ( jbuf.toArray() );
}


//-----------------------------------------------------------------------
//   announceTensor_
//-----------------------------------------------------------------------

void TbFiller::announceTensor_

  ( const String&  name ) 

{
  if      ( rank_ == 1 )
  {
    typeNames_.pushBack ( name + "_xx" );
  }
  else if ( rank_ == 2 )
  {
    typeNames_.pushBack ( name + "_xx" );
    typeNames_.pushBack ( name + "_yy" );
    typeNames_.pushBack ( name + "_zz" );
    typeNames_.pushBack ( name + "_xy" );
  }
  else
  {
    typeNames_.pushBack ( name + "_xx" );
    typeNames_.pushBack ( name + "_yy" );
    typeNames_.pushBack ( name + "_zz" );
    typeNames_.pushBack ( name + "_xy" );
    typeNames_.pushBack ( name + "_yz" );
    typeNames_.pushBack ( name + "_xz" );
  }
}

//-----------------------------------------------------------------------
//   announceVector_
//-----------------------------------------------------------------------

void TbFiller::announceVector_

  ( const String&  name ) 

{
  if      ( rank_ == 1 )
  {
    typeNames_.pushBack ( name + "_x" );
  }
  else if ( rank_ == 2 )
  {
    typeNames_.pushBack ( name + "_x" );
    typeNames_.pushBack ( name + "_y" );
  }
  else
  {
    typeNames_.pushBack ( name + "_x" );
    typeNames_.pushBack ( name + "_y" );
    typeNames_.pushBack ( name + "_z" );
  }
}

//-----------------------------------------------------------------------
//   announceDiag_
//-----------------------------------------------------------------------

void TbFiller::announceDiag_

  ( const String&  name ) 

{
  if      ( rank_ == 1 )
  {
    typeNames_.pushBack ( name + "_11" );
  }
  else
  {
    typeNames_.pushBack ( name + "_11" );
    typeNames_.pushBack ( name + "_22" );
    typeNames_.pushBack ( name + "_33" );
  }
}

//-----------------------------------------------------------------------
//   permTri6
//-----------------------------------------------------------------------

void TbFiller::permTri6

  ( const Vector& weights,
    const Matrix& values )

{
  if ( permi_.size() != 3 )
  {
    permi_ .resize ( 3 );
    permj1_.resize ( 3 );
    permj2_.resize ( 3 );

    permi_ [0] = 0; permi_ [1] = 2; permi_ [2] = 4;
    permj1_[0] = 1; permj1_[1] = 3; permj1_[2] = 5;
    permj2_[0] = 5; permj2_[1] = 1; permj2_[2] = 3;
  }

  weights[permi_] = weights[permj1_] + weights[permj2_];
  
  select ( values, permi_, ALL ) = select ( values, permj1_, ALL )
                                 + select ( values, permj2_, ALL );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   print (for NamedSlice)
//-----------------------------------------------------------------------

void print

  (       jem::io::Writer&      out, 
    const TbFiller::NamedSlice& ns)

{
  print ( out, ns.name(), " [", ns.slice().first(), "," );
  print ( out, ns.slice().last(), ")" );
}

