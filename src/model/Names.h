#ifndef NAMES_H
#define NAMES_H


//-----------------------------------------------------------------------
//   class Names
//-----------------------------------------------------------------------


class Names
{
 public:

  static const char*    DOFS[4];

};


//-----------------------------------------------------------------------
//   class PropertyNames
//-----------------------------------------------------------------------


class PropertyNames
{
 public:

  static const char*    SOLVER;

  static const char*    MODEL;
  static const char*    ARC_FUNC;
  static const char*    MIN_INCR;
  static const char*    MAX_INCR;
  static const char*    OPT_ITER;
  static const char*    RED_STEP;
  static const char*    LOAD_INCR;
  static const char*    LOAD_SCALE;
  static const char*    WEIGHT_TABLE;
  static const char*    WEIGHT_TRACK;
  static const char*    WEIGHT_THRES;

};


typedef PropertyNames   PropNames;



#endif

