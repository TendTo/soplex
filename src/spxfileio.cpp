/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxfileio.cpp,v 1.4 2003/01/15 17:26:07 bzfkocht Exp $"

//#define DEBUGGING 1

#include <assert.h>

#include "spxdefines.h"
#include "spxsolver.h"
#include "spxfileio.h"

namespace soplex
{
bool SPxSolver::readBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
{
   METHOD( "SPxSolver::readBasisFile()" );

   spxifstream file(filename);

   if (!file)
      return false;
 
   return readBasis(file, rowNames, colNames);
}

bool SPxSolver::writeBasisFile(
   const char*    filename, 
   const NameSet& rowNames,
   const NameSet& colNames)
{
   METHOD( "SPxSolver::writeBasisFile()" );
   std::ofstream file(filename);

   if (!file)
      return false;
 
   writeBasis(file, rowNames, colNames);

   return true;
}


bool SPxSolver::readFile( 
   const char* filename, 
   NameSet*    rowNames,
   NameSet*    colNames, 
   DIdxSet*    intVars)
{
   METHOD( "SPxSolver::readFile()" );

   spxifstream file(filename);

   if (!file)
      return false;

   return read(file, rowNames, colNames, intVars);
}

bool SPxLP::readFile( 
   const char* filename, 
   NameSet*    rowNames,
   NameSet*    colNames, 
   DIdxSet*    intVars)
{
   METHOD( "SPxLP::readFile()" );

   spxifstream file(filename);

   if (!file)
      return false;

   return read(file, rowNames, colNames, intVars);
}

void SPxSolver::dumpFile(const char* filename) const
{
   METHOD( "SPxSolver::dumpFile()" );
   std::ofstream file(filename);

   if (file.good())
      file << *this;
}

} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
