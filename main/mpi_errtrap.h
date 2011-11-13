/*
 *
 * VERSION : $Revision: 1.2 $
 * DATE    : $Date: 2002/11/20 15:55:08 $
 * FILE    : mpi_errtrap.h
 * DESCR.  : Declare errhandler-type
 * AUTHOR  : J. Thyssen - 980709
 *
 */

#ifndef _MPI_ERRTRAP_H
#define _MPI_ERRTRAP_H

#if defined (VAR_MPI)

#include <mpi.h>
typedef MPI_Errhandler errhandler_type;

#else /* VAR_MPI */

typedef int errhandler_type;

#endif /* VAR_MPI */

#endif /* _MPI_ERRTRAP_H */
