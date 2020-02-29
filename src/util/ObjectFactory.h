/*
 *
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *
 *  This class implements a rather simple but generic object
 *  factory, following a column in gamedev.net.
 *
 *  Usage:
 *
 *   1. Define a factory for a concrete super class: 
 *
 *         ObjectFactory<Material,int>  MaterialFactory;
 *
 *   2. Register concrete classes to this factory
 *
 *         MaterialFactory.register<HookeMaterial> (1);
 *
 *   3. When need to create an instance of HookeMaterial, then
 *
 *         Material* HookeMaterial = MaterialFactory.create(1);
 *
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 August 2008
 *
 */


#ifndef OBJECT_FACTORY_H
#define OBJECT_FACTORY_H


#include <map>

using std::map;

// -----------------------------------------------------------------------
//   function CreateObject
// -----------------------------------------------------------------------

template <class BaseClass,
          class ConcreteClass >

BaseClass* CreateObject ()
{
  return new ConcreteClass ();
};

// -----------------------------------------------------------------------
//   class ObjectFactory
// -----------------------------------------------------------------------

template <class    BaseClass,
          typename Identifier >

class ObjectFactory
{
  public:

    // -------------------------------------------------------------------
    //   some typedefs
    // -------------------------------------------------------------------

    typedef  BaseClass* *()  createObjectFunc;

    typedef  map<Identifier, createObjectFunc>::const_iterator cit;
    typedef  map<Identifier, createObjectFunc>::iterator       it;


    // -------------------------------------------------------------------
    //     register
    // -------------------------------------------------------------------

    template <class ConcreteClass>

    bool   register ( Identifier id )
    {
      cit iter = id2CreateMethodMap_.find ( id );

      if ( iter != id2CreateMethodMap_.end () )
      {
        id2CreateMethodMap_[id] = &CreateObject<BaseClass,ConcreteClass> ( );

        return true;
      }
    }

    // -------------------------------------------------------------------
    //     create (id)
    // -------------------------------------------------------------------

    BaseClass* create ( Identifier id )
    {
      cit iter = id2CreateMethodMap_.find ( id );

      if ( iter != id2CreateMethodMap_.end () )
      {
        return id2CreateMethodMap_[id].second ( );
      }
    }

  protected:

  private:

    map<Identifier, createObjectFunc>   id2CreateMethodMap_;
};

#endif

