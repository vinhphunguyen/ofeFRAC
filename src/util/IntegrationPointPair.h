#ifndef INTEGRATION_POINT_PAIR_H
#define INTEGRATION_POINT_PAIR_H

//=======================================================================
//   class IPointPair
//=======================================================================

// The helper class IPointPair stores information about a pair of
// integration points that are geometrically close to each other.

class IPointPair
{
 public:

  inline                    IPointPair  ();

 public:

  int                       ipoint;
  int                       jpoint;
  double                    dist;
  double                    weight1;
  double                    weight2;

};


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


inline IPointPair::IPointPair ()
{
  ipoint  = -1;
  jpoint  = -1;
  dist    =  0.0;
  weight1 =  0.0;
  weight2 =  0.0;
}


#endif

