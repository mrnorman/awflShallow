
#ifndef _IO_H_
#define _IO_H_

#include "types.h"
#include "Array.h"

void ncwrap( int ierr , int line );
void output_init( str_par &par, str_dom &dom, str_dyn &dyn, str_stat &stat );

#endif
