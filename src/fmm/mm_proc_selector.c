/* Here we exploit the capability of C to define a 
 * pointer to a procedure that is selected just once
 * but called several times.
 * Hence we avoid repeated IF tests and floating flags that
 * F90 would demand.  */ 

/* Bug fix: dec. 2002.
   renamed mm_selected_t_contractor_ to emm_selected_t_contractor_
   to avoid Portland compiler failure, and Intel Compiler warnings.
   reason not yet determined!! */ 

/* NOTICE: problems with this C interface encountered when compiling 
           with "checkbounds -C" option with some compilers!  */

#include <stdio.h>

typedef void (*arg1)(void*);
typedef void (*arg2)(void*, void*);
typedef void (*arg3)(void*, void*, void*);
typedef void (*arg4)(void*, void*, void*, void*);
typedef void (*arg5)(void*, void*, void*, void*, void*);

static arg1 T_contractor;
static arg1 W_contractor;
static arg2 T_buffer;
static arg2 W_buffer;
static arg4 stored_T_tester;
static arg3 stored_test_pq;
static arg4 T_mould1;
static arg3 T_mould2;
static arg3 T_mould3;
static arg5 T_mould4;
static arg4 stored_WS_para;

/*---------------------------------------------------------------------------*/

void mm_store_t_contractor_(arg1 T_contractor_in)
{
  T_contractor = T_contractor_in;
}
void mm_selected_t_contractor_(void*a)
{
  T_contractor(a);
}

void mm_store_t_buffer_(arg2 T_buffer_in)
{
  T_buffer = T_buffer_in;
}
void mm_selected_t_buffer_(void*a, void*b)
{
  T_buffer(a,b);
}

/*---------------------------------------------------------------------------*/

void mm_store_w_contractor_(arg1 W_contractor_in)
{
  W_contractor = W_contractor_in;
}
void mm_selected_w_contractor_(void*a)
{
  W_contractor(a);
}

void mm_store_w_buffer_(arg2 W_buffer_in)
{
  W_buffer = W_buffer_in;
}
void mm_selected_w_buffer_(void*a, void*b)
{
  W_buffer(a,b);
}

/*---------------------------------------------------------------------------*/

void mm_store_t_tester_(arg4 tester_in)
{
  stored_T_tester = tester_in;
}
void mm_test_and_buffer_t_pair_(void*a, void*b, void*c, void*d)
{
  stored_T_tester(a,b,c,d);
}

void mm_store_test_(arg3 test_pq_in)
{
  stored_test_pq = test_pq_in;
}
void mm_included_pair_(void*a, void*b, void*c)
{
  stored_test_pq(a,b,c);
}

/*---------------------------------------------------------------------------*/

void mm_store_t_pair_mould1_(arg4 T_mould_in)
{
  T_mould1 = T_mould_in;
}
void mm_store_t_pair_mould2_(arg3 T_mould_in)
{
  T_mould2 = T_mould_in;
}
void mm_store_t_pair_mould3_(arg3 T_mould_in)
{
  T_mould3 = T_mould_in;
}
void mm_store_t_pair_mould4_(arg5 T_mould_in)
{
  T_mould4 = T_mould_in;
}
void mm_stored_t_pair_mould_(void*LHS, void*RHS, void*id, void*wt, void*T_pair) 
{
  T_mould1(LHS,RHS,id,T_pair);
  T_mould2(LHS,id,T_pair);
  T_mould3(RHS,id,T_pair);
  T_mould4(LHS,RHS,id,wt,T_pair);
}

/*---------------------------------------------------------------------------*/

void mm_store_ws_para_(arg4 WS_para_in)
{
  stored_WS_para = WS_para_in;
}
void mm_get_ws_para_(void*a, void*b, void*c, void*d)
{
  stored_WS_para(a,b,c,d);
}

/*---------------------------------------------------------------------------*/

