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
#pragma ident "@(#) $Id: spxfileio.h,v 1.1 2003/01/15 17:26:07 bzfkocht Exp $"

/**@file  spxfileio.h
 * @brief declaration of types for file output
 *
 * This is to make the use of compressed input files transparent
 * when programming.
 */
#ifndef _SPXFILEIO_H_
#define _SPXFILEIO_H_

#include <iostream>
#include <fstream>

/*-----------------------------------------------------------------------------
 * compressed file support
 *-----------------------------------------------------------------------------
 */
#ifdef WITH_ZLIB
#include "gzstream.h"
#endif // WITH_GSZSTREAM

namespace soplex
{
#ifdef WITH_ZLIB
   typedef gzstream::igzstream spxifstream;
#else
   typedef std::ifstream spxifstream;
#endif // WITH_ZLIB

} // namespace soplex
#endif // _SPXFILEIO_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
