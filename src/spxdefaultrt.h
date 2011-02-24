/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: spxdefaultrt.h,v 1.21 2011/02/24 13:53:39 bzfgleix Exp $"

/**@file  spxdefaultrt.h
 * @brief Textbook ratio test for SoPlex.
 */
#ifndef _SPXDEFAULTRT_H_
#define _SPXDEFAULTRT_H_


#include <assert.h>

#include "spxdefines.h"
#include "spxratiotester.h"

namespace soplex
{

/**@brief   Textbook ratio test for SoPlex.
   @ingroup Algo
   
   Class SPxDefaultRT provides an implementation of the textbook ratio test
   as a derived class of SPxRatioTester. This class is not intended for
   reliably solving LPs (even though it does the job for ``numerically simple''
   LPs). Instead, it should serve as a demonstration of how to write ratio
   tester classes.

   See SPxRatioTester for a class documentation.
*/
class SPxDefaultRT : public SPxRatioTester
{
public:

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// default constructor
   SPxDefaultRT() 
      : SPxRatioTester("Default")
   {}
   /// copy constructor
   SPxDefaultRT(const SPxDefaultRT& old) 
      : SPxRatioTester(old)
   {}
   /// assignment operator
   SPxDefaultRT& operator=( const SPxDefaultRT& rhs)
   {
      if(this != &rhs)
      {
         SPxRatioTester::operator=(rhs);
      }

      return *this;
   }
   /// destructor
   virtual ~SPxDefaultRT()
   {}
   /// clone function for polymorphism
   inline virtual SPxRatioTester* clone() const
   {
      return new SPxDefaultRT(*this);
   }
   //@}

   //-------------------------------------
   /**@name Select enter/leave */
   //@{
   /// 
   virtual int selectLeave(Real& val, SPxId);
   ///
   virtual SPxId selectEnter(Real& val, int);
};

} // namespace soplex
#endif // _SPXDEFAULTRT_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
