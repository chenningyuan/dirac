#ifndef DISPATCH_H
#define DISPATCH_H

// Some macros to help with connecting dynamic code to template code.

// Build a table with function pointers to function<0> .. function<N-1>
// fpointer is the typedef'd type of the function pointer, i.e. 
// typedef void (*fpointer)(int)
#define DISPATCH_TABLE(table,fpointer,function,N) \
template< int i > fpointer *function##dispatch_init(fpointer function_table[]) \
{ \
  function_table[i] = function<i>; \
  return function##dispatch_init< i - 1 >(function_table); \
} \
template<> fpointer *function##dispatch_init<-1>(fpointer function_table[]) \
{  \
  return function_table; \
} \
const fpointer *table = function##dispatch_init<N>(new fpointer[N])

// Macro to generate lists f<0>, f<1>, ..
#define FUNLIST(f,N) FUNLIST_OFF(f,0,N)
#define FUNLIST_OFF(f,OFF,N) FUNLIST_OFF##N(f,OFF)
#define FUNLIST_OFF0(f,OFF) f<OFF>
#define FUNLIST_OFF1(f,OFF) f<OFF>, f<OFF+1>
#define FUNLIST_OFF2(f,OFF) f<OFF>, FUNLIST_OFF1(f,OFF + 1)
#define FUNLIST_OFF3(f,OFF) f<OFF>, FUNLIST_OFF2(f,OFF + 1)
#define FUNLIST_OFF4(f,OFF) f<OFF>, FUNLIST_OFF3(f,OFF + 1)
#define FUNLIST_OFF5(f,OFF) f<OFF>, FUNLIST_OFF4(f,OFF + 1)
#define FUNLIST_OFF6(f,OFF) f<OFF>, FUNLIST_OFF5(f,OFF + 1)
#define FUNLIST_OFF7(f,OFF) f<OFF>, FUNLIST_OFF6(f,OFF + 1)
#define FUNLIST_OFF8(f,OFF) f<OFF>, FUNLIST_OFF7(f,OFF + 1)
#define FUNLIST_OFF9(f,OFF) f<OFF>, FUNLIST_OFF8(f,OFF + 1)
#define FUNLIST_OFF10(f,OFF) f<OFF>, FUNLIST_OFF9(f,OFF + 1)

#endif
