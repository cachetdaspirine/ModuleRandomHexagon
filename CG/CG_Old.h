#ifndef CG_H
#define CG_H

#include <blitz/array.h>
#include "Ham.h"
#include "nr3.h"

class Fiber;
class Lattice;
class Parameter;

using namespace blitz;

class CG
{
 public:
  CG(Parameter* parameter);
  void Evolv(Lattice* lattice, Fiber* fiber);
  
 private:
  int Z,Lx,Ly,Lz;
  double eps,gamma;
  
  
}


#endif CG_H
