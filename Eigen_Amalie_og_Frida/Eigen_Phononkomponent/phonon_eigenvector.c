/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: phonon_eigenvector.instr (template_body_simple)
 * Date:       Tue Mar  5 06:38:13 2024
 * File:       ./phonon_eigenvector.c
 * CFLAGS=
 */

#define MCCODE_STRING ""
#define FLAVOR        "mcstas"
#define FLAVOR_UPPER  "MCSTAS"

#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED

#include <string.h>

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];
#define MCCODE_BASE_TYPES

#ifdef OPENACC
#undef MC_TRACE_ENABLED
#endif

#ifndef MC_NUSERVAR
#define MC_NUSERVAR 10
#endif

/* Particle JUMP control logic */
struct particle_logic_struct {
int dummy;
};

struct _struct_particle {
  double x,y,z; /* position [m] */
  double vx,vy,vz; /* velocity [m/s] */
  double sx,sy,sz; /* spin [0-1] */
  int mcgravitation; /* gravity-state */
  void *mcMagnet;    /* precession-state */
  int allow_backprop; /* allow backprop */
  unsigned long randstate[7];
  double t, p;    /* time, event weight */
  long long _uid;  /* event ID */
  long _index;     /* component index where to send this event */
  long _absorbed;  /* flag set to TRUE when this event is to be removed/ignored */
  long _scattered; /* flag set to TRUE when this event has interacted with the last component instance */
  long _restore;   /* set to true if neutron event must be restored */
  long flag_nocoordschange;   /* set to true if particle is jumping */
  struct particle_logic_struct _logic;
};
typedef struct _struct_particle _class_particle;

_class_particle _particle_global_randnbuse_var;
_class_particle* _particle = &_particle_global_randnbuse_var;

#pragma acc routine
_class_particle mcgenstate(void);
#pragma acc routine
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
			   double t, double sx, double sy, double sz, double p, int mcgravitation, void *mcMagnet, int mcallowbackprop);

extern int mcgravitation;      /* flag to enable gravitation */
#pragma acc declare create ( mcgravitation )
int mcallowbackprop;        
#pragma acc declare create ( mcallowbackprop )

_class_particle mcgenstate(void) {
  _class_particle particle = mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, mcgravitation, NULL, mcallowbackprop);
  return(particle);
}
/*Generated user variable handlers:*/

#pragma acc routine
double particle_getvar(_class_particle *p, char *name, int *suc);

#ifdef OPENACC
#pragma acc routine
int str_comp(char *str1, char *str2);
#endif

double particle_getvar(_class_particle *p, char *name, int *suc){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int s=1;
  double rval=0;
  if(!str_comp("x",name)){rval=p->x;s=0;}
  if(!str_comp("y",name)){rval=p->y;s=0;}
  if(!str_comp("z",name)){rval=p->z;s=0;}
  if(!str_comp("vx",name)){rval=p->vx;s=0;}
  if(!str_comp("vy",name)){rval=p->vy;s=0;}
  if(!str_comp("vz",name)){rval=p->vz;s=0;}
  if(!str_comp("sx",name)){rval=p->sx;s=0;}
  if(!str_comp("sy",name)){rval=p->sy;s=0;}
  if(!str_comp("sz",name)){rval=p->sz;s=0;}
  if(!str_comp("t",name)){rval=p->t;s=0;}
  if(!str_comp("p",name)){rval=p->p;s=0;}
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
void* particle_getvar_void(_class_particle *p, char *name, int *suc);

#ifdef OPENACC
#pragma acc routine
int str_comp(char *str1, char *str2);
#endif

void* particle_getvar_void(_class_particle *p, char *name, int *suc){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int s=1;
  void* rval=0;
  if(!str_comp("x",name)) {rval=(void*)&(p->x); s=0;}
  if(!str_comp("y",name)) {rval=(void*)&(p->y); s=0;}
  if(!str_comp("z",name)) {rval=(void*)&(p->z); s=0;}
  if(!str_comp("vx",name)){rval=(void*)&(p->vx);s=0;}
  if(!str_comp("vy",name)){rval=(void*)&(p->vy);s=0;}
  if(!str_comp("vz",name)){rval=(void*)&(p->vz);s=0;}
  if(!str_comp("sx",name)){rval=(void*)&(p->sx);s=0;}
  if(!str_comp("sy",name)){rval=(void*)&(p->sy);s=0;}
  if(!str_comp("sz",name)){rval=(void*)&(p->sz);s=0;}
  if(!str_comp("t",name)) {rval=(void*)&(p->t); s=0;}
  if(!str_comp("p",name)) {rval=(void*)&(p->p); s=0;}
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
int particle_setvar_void(_class_particle *, char *, void*);

int particle_setvar_void(_class_particle *p, char *name, void* value){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int rval=1;
  if(!str_comp("x",name)) {memcpy(&(p->x),  value, sizeof(double)); rval=0;}
  if(!str_comp("y",name)) {memcpy(&(p->y),  value, sizeof(double)); rval=0;}
  if(!str_comp("z",name)) {memcpy(&(p->z),  value, sizeof(double)); rval=0;}
  if(!str_comp("vx",name)){memcpy(&(p->vx), value, sizeof(double)); rval=0;}
  if(!str_comp("vy",name)){memcpy(&(p->vy), value, sizeof(double)); rval=0;}
  if(!str_comp("vz",name)){memcpy(&(p->vz), value, sizeof(double)); rval=0;}
  if(!str_comp("sx",name)){memcpy(&(p->sx), value, sizeof(double)); rval=0;}
  if(!str_comp("sy",name)){memcpy(&(p->sy), value, sizeof(double)); rval=0;}
  if(!str_comp("sz",name)){memcpy(&(p->sz), value, sizeof(double)); rval=0;}
  if(!str_comp("p",name)) {memcpy(&(p->p),  value, sizeof(double)); rval=0;}
  if(!str_comp("t",name)) {memcpy(&(p->t),  value, sizeof(double)); rval=0;}
  return rval;
}

#pragma acc routine
int particle_setvar_void_array(_class_particle *, char *, void*, int);

int particle_setvar_void_array(_class_particle *p, char *name, void* value, int elements){
#ifndef OPENACC
#define str_comp strcmp
#endif
  int rval=1;
  return rval;
}

#pragma acc routine
void particle_restore(_class_particle *p, _class_particle *p0);

void particle_restore(_class_particle *p, _class_particle *p0) {
  p->x  = p0->x;  p->y  = p0->y;  p->z  = p0->z;
  p->vx = p0->vx; p->vy = p0->vy; p->vz = p0->vz;
  p->sx = p0->sx; p->sy = p0->sy; p->sz = p0->sz;
  p->t = p0->t;  p->p  = p0->p;
  p->_absorbed=0; p->_restore=0;
}

#pragma acc routine
double particle_getuservar_byid(_class_particle *p, int id, int *suc){
  int s=1;
  double rval=0;
  switch(id){
  }
  if (suc!=0x0) {*suc=s;}
  return rval;
}

#pragma acc routine
void particle_uservar_init(_class_particle *p){
}

#define MC_EMBEDDED_RUNTIME
/* embedding file "mccode-r.h" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 3.2
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int numipar;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <float.h>
#include <inttypes.h>
#include <stdint.h>
#ifdef OPENACC
#include <openacc.h>
#ifndef GCCOFFLOAD
#include <accelmath.h>
#else
#include <math.h>
#endif
#pragma acc routine
int noprintf();
#pragma acc routine
size_t str_len(const char *s);
#else
#include <math.h>
#endif

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#  define mcstatic
#else
#  define mcstatic
#endif

#ifdef __dest_os
#  if (__dest_os == __mac_os)
#    define MAC
#  endif
#endif

#ifdef __FreeBSD__
#  define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#  define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#  include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#  ifdef WIN32
#    define MC_PATHSEP_C '\\'
#    define MC_PATHSEP_S "\\"
#  else  /* !WIN32 */
#    define MC_PATHSEP_C '/'
#    define MC_PATHSEP_S "/"
#  endif /* !WIN32 */
#endif /* MC_PATHSEP_C */



/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#  define MCCODE_STRING "McStas 3.2 - Dec. 12, 2022"
#endif

#ifndef MCCODE_DATE
#  define MCCODE_DATE "Dec. 12, 2022"
#endif

#ifndef MCCODE_VERSION
#  define MCCODE_VERSION "3.2"
#endif

#ifndef MCCODE_NAME
#  define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#  define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_PARTICLE_CODE
#  define MCCODE_PARTICLE_CODE 2112
#endif

#ifndef MCCODE_LIBENV
#  define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#  define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifdef MAC
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#if (USE_MPI == 0)
#  undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifdef OPENACC  /* default is to disable signals with PGI/OpenACC */
#  ifndef NOSIGNALS
#    define NOSIGNALS 1
#  endif
#endif

#ifndef OPENACC
#  ifndef USE_OFF  /* default is to enable OFF when not using PGI/OpenACC */
#    define USE_OFF
#  endif
#  ifndef CPUFUNNEL  /* allow to enable FUNNEL-mode on CPU */
#  ifdef FUNNEL      /* by default disable FUNNEL-mode when not using PGI/OpenACC */
#    undef FUNNEL
#  endif
#  endif
#endif

#if (NOSIGNALS == 0)
#  undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_int,
    instr_type_string, instr_type_char,
    instr_type_vector, instr_type_double
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

#ifndef MCCODE_BASE_TYPES
typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];
#endif

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[];         /* list of instrument parameters */
extern int    numipar;                                    /* number of instrument parameters */
extern char   instrument_name[], instrument_source[]; /* instrument name and filename */
extern char  *instrument_exe;                           /* executable path = argv[0] or NULL */
extern char   instrument_code[];                        /* contains the initial 'instr' file */

#ifndef MC_ANCIENT_COMPATIBILITY
extern int traceenabled, defaultmain;
#endif
#endif


/* Useful macros ============================================================ */


/* SECTION: Dynamic Arrays */
typedef int* IArray1d;
IArray1d create_iarr1d(int n);
void destroy_iarr1d(IArray1d a);

typedef int** IArray2d;
IArray2d create_iarr2d(int nx, int ny);
void destroy_iarr2d(IArray2d a);

typedef int*** IArray3d;
IArray3d create_iarr3d(int nx, int ny, int nz);
void destroy_iarr3d(IArray3d a);

typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);


/* MPI stuff */
#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#ifndef OPENACC    /* ... but only if we are not also running on GPU */
#undef NOSIGNALS
#endif
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */


#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
#pragma acc routine
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount-1 */

/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
char  *mcsig_message;
#define SIG_MESSAGE(msg) mcsig_message=(char *)(msg);
#else
#define SIG_MESSAGE(...)
#endif /* !NOSIGNALS */


/* Useful macros and constants ============================================== */


#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif


#  ifndef M_E
#    define M_E        2.71828182845904523536  // e
#  endif
#  ifndef M_LOG2E
#    define M_LOG2E    1.44269504088896340736  //  log2(e)
#  endif
#  ifndef M_LOG10E
#    define M_LOG10E   0.434294481903251827651 //  log10(e)
#  endif
#  ifndef M_LN2
#    define M_LN2      0.693147180559945309417 //  ln(2)
#  endif
#  ifndef M_LN10
#    define M_LN10     2.30258509299404568402  //  ln(10)
#  endif
#  ifndef M_PI
#    define M_PI       3.14159265358979323846  //  pi
#  endif
#  ifndef PI
#    define PI       M_PI                      //  pi - also used in some places
#  endif
#  ifndef M_PI_2
#    define M_PI_2     1.57079632679489661923  //  pi/2
#  endif
#  ifndef M_PI_4
#    define M_PI_4     0.785398163397448309616 //  pi/4
#  endif
#  ifndef M_1_PI
#    define M_1_PI     0.318309886183790671538 //  1/pi
#  endif
#  ifndef M_2_PI
#    define M_2_PI     0.636619772367581343076 //  2/pi
#  endif
#  ifndef M_2_SQRTPI
#    define M_2_SQRTPI 1.12837916709551257390  //  2/sqrt(pi)
#  endif
#  ifndef M_SQRT2
#    define M_SQRT2    1.41421356237309504880  //  sqrt(2)
#  endif
#  ifndef M_SQRT1_2
#    define M_SQRT1_2  0.707106781186547524401 //  1/sqrt(2)
#  endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


#define UNSET nan("0x6E6F74736574")
int nans_match(double, double);
int is_unset(double);
int is_valid(double);
int is_set(double);
int all_unset(int n, ...);
int all_set(int n, ...);
int any_unset(int n, ...);
int any_set(int n, ...);


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) (instrument->_position_absolute[index])
#define POS_R_COMP_INDEX(index) (instrument->_position_relative[index])

/* setting parameters based COMP_GETPAR (returned as pointer)         */
/* compname must be given as a string, type and par are symbols.      */
#define COMP_GETPAR3(type, compname, par) \
    &( ((_class_ ## type ##_parameters *) _getvar_parameters(compname))->par )
/* the body of this function depends on component instances, and is cogen'd */
void* _getvar_parameters(char* compname);

int _getcomp_index(char* compname);

/* Note: The two-stage approach to COMP_GETPAR is NOT redundant; without it,
* after #define C sample, COMP_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of COMP_GETPAR requires that we use sometimes bare names...
* NOTE: This can ONLY be used in instrument descriptions, not components.
*/
#define COMP_GETPAR2(comp, par) (_ ## comp ## _var._parameters.par)
#define COMP_GETPAR(comp, par) COMP_GETPAR2(comp,par)

#define INSTRUMENT_GETPAR(par) (_instrument_var._parameters.par)

/* Current component name, index, position and orientation */
/* These macros work because, using class-based functions, "comp" is usually
*  the local variable of the active/current component. */
#define INDEX_CURRENT_COMP (_comp->_index)
#define NAME_CURRENT_COMP (_comp->_name)
#define TYPE_CURRENT_COMP (_comp->_type)
#define POS_A_CURRENT_COMP (_comp->_position_absolute)
#define POS_R_CURRENT_COMP (_comp->_position_relative)
#define ROT_A_CURRENT_COMP (_comp->_rotation_absolute)
#define ROT_R_CURRENT_COMP (_comp->_rotation_relative)

#define NAME_INSTRUMENT (instrument->_name)


/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define DEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", instrument_name, instrument_source); }
#define DEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  printf("Component %30s AT (%g,%g,%g)\n", name, c.x, c.y, c.z); \
  }
#define DEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define DEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define DEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define DEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define DEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define DEBUG_INSTR()
#define DEBUG_COMPONENT(name,c,t)
#define DEBUG_INSTR_END()
#define DEBUG_ENTER()
#define DEBUG_COMP(c)
#define DEBUG_LEAVE()
#define DEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);


/* random number generation. ================================================ */

/* available random number generators */
#define _RNG_ALG_MT         1
#define _RNG_ALG_KISS       2

/* selection of random number generator */
#ifndef RNG_ALG
#  define RNG_ALG  _RNG_ALG_KISS
#endif


#if RNG_ALG == _RNG_ALG_MT  // MT (currently not functional for gpu)
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define randstate_t unsigned long // this could be anything
#  define RANDSTATE_LEN 1
#  define srandom(seed) mt_srandom_empty()
#  define random() mt_random()
#  define _random() mt_random()
#elif RNG_ALG == _RNG_ALG_KISS  // KISS
#  ifndef ULONG_MAX
#    define ULONG_MAX ((unsigned long)0xffffffffffffffffUL)
#  endif
#  define MC_RAND_MAX ULONG_MAX
#  define randstate_t unsigned long
#  define RANDSTATE_LEN 7
#  define srandom(seed) kiss_srandom(_particle->randstate, seed)
#  define random() kiss_random(_particle->randstate)
#  define _random() kiss_random(state)
#endif

#pragma acc routine
double _randnorm2(randstate_t* state);


// component writers interface
#define randnorm() _randnorm2(_particle->randstate) // NOTE: can not use _randnorm on gpu
#define rand01() _rand01(_particle->randstate)
#define randpm1() _randpm1(_particle->randstate)
#define rand0max(p1) _rand0max(p1, _particle->randstate)
#define randminmax(p1, p2) _randminmax(p1, p2, _particle->randstate)
#define randtriangle() _randtriangle(_particle->randstate)

// Mersenne Twister rng
unsigned long mt_random(void);
void mt_srandom (unsigned long x);
void mt_srandom_empty();

// KISS rng
#pragma acc routine
unsigned long *kiss_srandom(unsigned long state[7], unsigned long seed);
#pragma acc routine
unsigned long kiss_random(unsigned long state[7]);

// Scrambler / hash function
#pragma acc routine seq
randstate_t _hash(randstate_t x);

// internal RNG (transforms) interface
#pragma acc routine
double _rand01(randstate_t* state);
#pragma acc routine
double _randpm1(randstate_t* state);
#pragma acc routine
double _rand0max(double max, randstate_t* state);
#pragma acc routine
double _randminmax(double min, double max, randstate_t* state);
#pragma acc routine
double _randtriangle(randstate_t* state);


#ifdef USE_OPENCL
#include "opencl-lib.h"
#include "opencl-lib.c"
#endif

#ifndef DANSE
int init(void);
int raytrace(_class_particle*);
int save(FILE *);
int finally(void);
int display(void);
#endif


/* GPU related algorithms =================================================== */

/*
*  Divide-and-conquer strategy for parallel sort absorbed last.
*/
#ifdef FUNNEL
long sort_absorb_last(_class_particle* particles, _class_particle* pbuffer, long len, long buffer_len, long flag_split, long* multiplier);
#endif
long sort_absorb_last_serial(_class_particle* particles, long len);


/* simple vector algebra ==================================================== */


#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
#pragma acc routine seq
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#pragma acc routine seq
mcstatic void norm_func(double *x, double *y, double *z);
#define NORM(x,y,z)	norm_func(&x, &y, &z)

#pragma acc routine seq
void normal_vec(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

#pragma acc routine
Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
#pragma acc routine
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
#pragma acc routine
Coords coords_add(Coords a, Coords b);
#pragma acc routine
Coords coords_sub(Coords a, Coords b);
#pragma acc routine
Coords coords_neg(Coords a);
#pragma acc routine
Coords coords_scale(Coords b, double scale);
#pragma acc routine
double coords_sp(Coords a, Coords b);
#pragma acc routine
Coords coords_xp(Coords b, Coords c);
#pragma acc routine
double coords_len(Coords a);
#pragma acc routine seq
void   coords_print(Coords a);
#pragma acc routine seq
mcstatic void coords_norm(Coords* c);

#pragma acc routine seq
void rot_set_rotation(Rotation t, double phx, double phy, double phz);
#pragma acc routine seq
int  rot_test_identity(Rotation t);
#pragma acc routine seq
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
#pragma acc routine seq
void rot_copy(Rotation dest, Rotation src);
#pragma acc routine seq
void rot_transpose(Rotation src, Rotation dst);
#pragma acc routine seq
Coords rot_apply(Rotation t, Coords a);

#pragma acc routine seq
void mccoordschange(Coords a, Rotation t, _class_particle *particle);
#pragma acc routine seq
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters
is no longer equal */

_class_particle mcgenstate(void);

// trajectory/shape intersection routines
#pragma acc routine seq
int inside_rectangle(double, double, double, double);
#pragma acc routine seq
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
      double vx, double vy, double vz, double dx, double dy, double dz);
#pragma acc routine seq
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r, double h);
#pragma acc routine seq
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
      double vx, double vy, double vz, double r);
// second order equation roots
#pragma acc routine seq
int solve_2nd_order(double *t1, double *t2,
      double A,  double B,  double C);

// random vector generation to shape
// defines silently introducing _particle as the last argument
#define randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, radius) \
  _randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, radius, _particle)
#define randvec_target_rect_angular(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A) \
  _randvec_target_rect_angular(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, _particle)
#define randvec_target_rect_real(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order) \
  _randvec_target_rect_real(xo, yo, zo, solid_angle, xi, yi, zi, height, width, A, lx, ly, lz, order, _particle)
// defines forwarding to "inner" functions
#define randvec_target_sphere randvec_target_circle
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9) \
  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
// headers for randvec
#pragma acc routine seq
void _randvec_target_circle(double *xo, double *yo, double *zo,
  double *solid_angle, double xi, double yi, double zi, double radius,
  _class_particle* _particle);
#pragma acc routine seq
void _randvec_target_rect_angular(double *xo, double *yo, double *zo,
  double *solid_angle, double xi, double yi, double zi, double height,
  double width, Rotation A,
  _class_particle* _particle);
#pragma acc routine seq
void _randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
  double xi, double yi, double zi, double height, double width, Rotation A,
  double lx, double ly, double lz, int order,
  _class_particle* _particle);


// this is the main()
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif


/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *dirname             = NULL;      /* name of output directory */
static   char *siminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * siminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *siminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle);

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

/* embedding file "mcstas-r.h" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char instrument_name[], instrument_source[];
*   int traceenabled, defaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  instrument.counter_AbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundent with mcstas.h */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER0 do {DEBUG_SCATTER(); SCATTERED++;} while(0)
#define SCATTER SCATTER0

#define JUMPTOCOMP(comp) mcneutron->_index = INDEX_COMP(comp);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  } while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(x, y, z, t, vx, vy, vz, \
               &sx, &sy, &sz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    x += vx*(dt); \
    y += vy*(dt); \
    z += vz*(dt); \
    t += (dt); \
    if (isnan(p) || isinf(p)) { ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { ABSORB; }\
    if (mcMagnet) /*printf("Spin precession gravity\n")*/; \
    x  += vx*(dt) + (Ax)*(dt)*(dt)/2; \
    y  += vy*(dt) + (Ay)*(dt)*(dt)/2; \
    z  += vz*(dt) + (Az)*(dt)*(dt)/2; \
    vx += (Ax)*(dt); \
    vy += (Ay)*(dt); \
    vz += (Az)*(dt); \
    t  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; ABSORB; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -vz, -z); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); z=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(vz == 0) { ABSORB; }; \
    mc_dt = -z/vz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    z = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -vx, -x); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); x=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(vx == 0) { ABSORB; }; \
    mc_dt = -x/vx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    x = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -vy, -y); \
    if (mc_ret) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); y=0;}\
    else if (mcallowbackprop == 0 && mc_dt < 0) { ABSORB; }; } \
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(vy == 0) { ABSORB; }; \
    mc_dt = -y/vy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { ABSORB; }; \
    mcPROP_DT(mc_dt); \
    y = 0; \
    DISALLOW_BACKPROP; \
  } while(0)


#ifdef DEBUG

#define DEBUG_STATE() if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define DEBUG_SCATTER() if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define DEBUG_STATE()
#define DEBUG_SCATTER()

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

/* embedding file "mccode-r.c" */

/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int traceenabled = 0;
int defaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
#pragma acc declare create ( mcseed )
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* OpenACC-related segmentation parameters: */
int vecsize = 128;
int numgangs = 7813;
long gpu_innerloop = 2147483647;

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
#ifdef MCDEFAULT_NCOUNT
mcstatic unsigned long long int mcncount             = MCDEFAULT_NCOUNT;
#else
mcstatic unsigned long long int mcncount             = 1000000;
#endif
#pragma acc declare create ( mcncount )
mcstatic unsigned long long int mcrun_num            = 0;
#pragma acc declare create ( mcrun_num )
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

#ifndef NX_COMPRESION
#define NX_COMPRESSION NX_COMP_NONE
#endif

/* String nullification on GPU and other replacements */
#ifdef OPENACC
int noprintf() {
  return 0;
}

int str_comp(char *str1, char *str2) {
  while (*str1 && *str1 == *str2) {
    str1++;
    str2++;
  }
  return (*str1 - *str2);
}

size_t str_len(const char *s)
{
  size_t len = 0;
  if(s != NULL)
  {
    while(*s != '\0')
    {
      ++len;
      ++s;
    }
  }
  return len;
}

#endif

/* SECTION: Predefine (component) parameters ================================= */

int nans_match(double a, double b){
  return (*(uint64_t*)&a == *(uint64_t*)&b);
}
int is_unset(double x){
  return nans_match(x, UNSET);
}
int is_set(double x){
  return !nans_match(x, UNSET);
}
int is_valid(double x){
  return !isnan(x)||is_unset(x);
}
int all_unset(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=1;
  for (int i=0; i<n; ++i) if(is_set(va_arg(ptr, double))) ret=0;
  va_end(ptr);
  return ret;
}
int all_set(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=1;
  for (int i=0; i<n; ++i) if(is_unset(va_arg(ptr, double))) ret=0;
  va_end(ptr);
  return ret;
}
int any_unset(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=0;
  for (int i=0; i<n; ++i) if(is_unset(va_arg(ptr, double))) ret=1;
  va_end(ptr);
  return ret;
}
int any_set(int n, ...){
  va_list ptr;
  va_start(ptr, n);
  int ret=0;
  for (int i=0; i<n; ++i) if(is_set(va_arg(ptr, double))) ret=1;
  va_end(ptr);
  return ret;
}


/* SECTION: Dynamic Arrays ================================================== */
IArray1d create_iarr1d(int n){
  IArray1d arr2d;
  arr2d = calloc(n, sizeof(int));
  return arr2d;
}
void destroy_iarr1d(IArray1d a){
  free(a);
}

IArray2d create_iarr2d(int nx, int ny){
  IArray2d arr2d;
  arr2d = calloc(nx, sizeof(int *));

  int *p1;
  p1 = calloc(nx*ny, sizeof(int));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_iarr2d(IArray2d a){
  free(a[0]);
  free(a);
}

IArray3d create_iarr3d(int nx, int ny, int nz){
  IArray3d arr3d;
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(int **));

  // d2
  int **p1;
  p1 = calloc(nx*ny, sizeof(int *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  int *p2;
  p2 = calloc(nx*ny*nz, sizeof(int));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}

void destroy_iarr3d(IArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}

void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}

void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;

  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}

void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}


/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

#ifndef STRACPY
/* this is a replacement to strncpy, but ensures that the copy ends with NULL */
/* http://stracpy.blogspot.fr/2011/04/stracpy-strncpy-replacement.html */
#define STRACPY
char *stracpy(char *destination, const char *source, size_t amount)
{
        if (!destination || !source || !amount) return(NULL);
        while(amount--)
          if((*destination++ = *source++) == '\0') break;
        *destination = '\0';
        return destination;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=dirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = dirname ? strlen(dirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, dirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within dirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);

  mem  = mcfull_file(name, ext); /* create dirname/name.ext */

  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;

  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);

    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+");

  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n",
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: detector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m)
    return(detector);

  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (detector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];

  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m)
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin);
      else x = 0;
      if (detector.n && detector.p)
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin);
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw)
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */

      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }

      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {

    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH,
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }

  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n",
      detector.component, strlen(detector.filename) ? detector.filename : "");

  return(detector);

} /* mcdetector_statistics */

/*******************************************************************************
* detector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=siminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR detector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", instrument_name, instrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, instrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", instrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }


  return(detector);
} /* detector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: siminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < numipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, dirname, MC_PATHSEP_C, siminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, instrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);

  fprintf(f, "%sTrace_enabled: %s\n", pre, traceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, defaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre,
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: siminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre,
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, instrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, dirname ? dirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < numipar; i++) {
      if (mcinputtable[i].par){
	/* Parameters with a default value */
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
	  (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
	  fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        /* ... and those without */
	}else{
	  fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
  fflush(f);
} /* mcruninfo_out */

/*******************************************************************************
* siminfo_out:    wrapper to fprintf(siminfo_file)
*******************************************************************************/
void siminfo_out(char *format, ...)
{
  va_list ap;

  if(siminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(siminfo_file, format, ap);
    va_end(ap);
  }
} /* siminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;

  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" :
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f,
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" :
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre,
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);

  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  /* Write data set information to simulation description file. */
  MPI_MASTER(
    siminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n",
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    siminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", siminfo_file, detector);
    siminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);

      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);

}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        siminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", siminfo_file, detector);
        siminfo_out("end data\n");

        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
        fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      }
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2,
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0,
          outfile, detector.istransposed);
      }
      fclose(outfile);

      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1,
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=CHAR_BUF_LENGTH; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;

  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);

  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);

    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }

  return(ret);

} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: siminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, instrument_name);

  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK)
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {

    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s",
      dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name,
      instrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   instrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < numipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }

      nxprintattr(f, "name",          instrument_name);
      nxprintf   (f, "name",          instrument_name);
      nxprintattr(f, "Source",        instrument_source);

      nxprintattr(f, "Trace_enabled", traceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  defaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );

      /* add instrument source code when available */
      buffer = mcinfo_readfile(instrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", instrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", instrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)",
          instrument_source, instrument_name);

      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",instrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", instrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        dirname && strlen(dirname) ? dirname : ".", MC_PATHSEP_S, siminfo_name);

      nxprintf   (f, "name",      "%s",     siminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",instrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", dirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif

      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < numipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */

      NXclosegroup(f); /* simulation */
    } /* NXsimulation */

    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];
  if (!f || !detector.m || mcdisable_output_files) return;

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */

    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" :
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" :
                 "xylimits", detector.limits);
      nxprintattr(f, "variables",
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");

      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */

} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[CHAR_BUF_LENGTH];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMPRESSION, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);

    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{

  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int fulldims[3]={detector.m,detector.n,detector.p};
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;

  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);

  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) fulldims[0] = NX_UNLIMITED;

  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, fulldims, NX_COMPRESSION, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */

  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",
        strlen(detector.filename) ? detector.filename : detector.component);
    else
      printf("Append:   \"%s\"\n",
	     strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }

  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);

  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[CHAR_BUF_LENGTH];

  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);

  strcpy_valid(data_name,
    detector.filename && strlen(detector.filename) ?
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (siminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {

      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar,
          1, detector.m, detector.xmin, detector.xmax);

        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar,
          2, detector.n, detector.ymin, detector.ymax);

        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar,
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */

      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);

      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */

  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may.
*   Then the number of recv/send is not constant along nodes, and simulation stalls.
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );

  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];

	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );

  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );

#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* siminfo_init:   open SIM and write header
*******************************************************************************/
FILE *siminfo_init(FILE *f)
{
  int exists=0;

  /* check format */
  if (!mcformat || !strlen(mcformat)
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE")
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }

  /* open the SIM file if not defined yet */
  if (siminfo_file || mcdisable_output_files)
    return (siminfo_file);

#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  siminfo_file = mcnew_file(siminfo_name, "h5", &exists);
    if(!siminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      siminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(siminfo_file); /* points to nxhandle */
  }
#endif

  /* write main description file (only MASTER) */
  MPI_MASTER(

  siminfo_file = mcnew_file(siminfo_name, "sim", &exists);
  if(!siminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    siminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    siminfo_out("%s simulation description file for %s.\n",
      MCCODE_NAME, instrument_name);
    siminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    siminfo_out("Program: %s\n\n", MCCODE_STRING);

    siminfo_out("begin instrument: %s\n", instrument_name);
    mcinfo_out(   "  ", siminfo_file);
    siminfo_out("end instrument\n");

    siminfo_out("\nbegin simulation: %s\n", dirname);
    mcruninfo_out("  ", siminfo_file);
    siminfo_out("end simulation\n");

  }
  return (siminfo_file);

  ); /* MPI_MASTER */

} /* siminfo_init */

/*******************************************************************************
*   siminfo_close:  close SIM
*******************************************************************************/
void siminfo_close()
{
#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
  if(siminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else {
#endif
      fclose(siminfo_file);
#ifdef USE_NEXUS
    }
#endif
#ifdef USE_MPI
  }
#endif
    siminfo_file = NULL;
  }
} /* siminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));

} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*
*   t:    title
*   xl:   x-label
*   yl:   y-label
*   xvar: measured variable length
*   x1:   x axus min
*   x2:   x axis max
*   n:    1d data vector lenght
*   p0:   pntr to start of data block#0
*   p1:   pntr to start of data block#1
*   p2:   pntr to start of data block#2
*   f:    filename
*
*   Not included in the macro, and here forwarded to detector_import:
*   c:    ?
*   posa: ?
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = detector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));

} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   Special case for list: master creates file first, then slaves append their
*   blocks without header-
*
*   t:    title
*   xl:   x-label
*   yl:   y-label
*   x1:   x axus min
*   x2:   x axis max
*   y1:   y axis min
*   y2:   y axis max
*   m:    dim 1 (x) size
*   n:    dim 2 (y) size
*   p0:   pntr to start of data block#0
*   p1:   pntr to start of data block#1
*   p2:   pntr to start of data block#2
*   f:    filename
*
*   Not included in the macro, and here forwarded to detector_import:
*   c:    ?
*   posa: ?
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];

  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = detector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));

} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;

  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,labs(m),1,labs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);

  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    dirname = dir;
  else
    dirname = dir+strlen("file://");


#ifdef USE_MPI
  if(mpi_node_rank == mpi_node_root) {
#endif
    if(mkdir(dirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
#ifdef USE_MPI
    }
#endif

  /* remove trailing PATHSEP (if any) */
  while (strlen(dirname) && dirname[strlen(dirname) - 1] == MC_PATHSEP_C)
    dirname[strlen(dirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", instrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", dirname ? dirname : ".");
  mcruninfo_out("  ", stdout);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays */
/* Within the TRACE scope we are now using _particle->uid directly */
unsigned long long int mcget_run_num() // shuld be (_class_particle* _particle) somehow
{
  /* This function only remains for the few cases outside TRACE where we need to know
     the number of simulated particles */
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(!mcseed) {
  //  srandom(mcseed);
  //} else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*PI/24),y,z+r*cos(i*2*PI/24),
                    x+r*sin((i+1)*2*PI/24),y,z+r*cos((i+1)*2*PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {
  #ifndef OPENACC
  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  #endif
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/* coords_test_zero: check if zero vector*/
int coords_test_zero(Coords a){
  return ( a.x==0 && a.y==0 && a.z==0 );
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}


/* SECTION: GPU algorithms ================================================== */


/*
*  Divide-and-conquer strategy for parallelizing this task: Sort absorbed
*  particles last.
*
*   particles:  the particle array, required to checking _absorbed
*   pbuffer:    same-size particle buffer array required for parallel sort
*   len:        sorting area-of-interest size (e.g. from previous calls)
*   buffer_len: total array size
*   flag_split: if set, multiply live particles into absorbed slots, up to buffer_len
*   multiplier: output arg, becomes the  SPLIT multiplier if flag_split is set
*/
#ifdef FUNNEL
long sort_absorb_last(_class_particle* particles, _class_particle* pbuffer, long len, long buffer_len, long flag_split, long* multiplier) {
  #define SAL_THREADS 1024 // num parallel sections
  if (len<SAL_THREADS) return sort_absorb_last_serial(particles, len);

  if (multiplier != NULL) *multiplier = -1; // set default out value for multiplier
  long newlen = 0;
  long los[SAL_THREADS]; // target array startidxs
  long lens[SAL_THREADS]; // target array sublens
  long l = floor(len/(SAL_THREADS-1)); // subproblem_len
  long ll = len - l*(SAL_THREADS-1); // last_subproblem_len

  // TODO: The l vs ll is too simplistic, since ll can become much larger
  // than l, resulting in idling. We should distribute lengths more evenly.

  // step 1: sort sub-arrays
  #pragma acc parallel loop present(particles, pbuffer)
  for (unsigned long tidx=0; tidx<SAL_THREADS; tidx++) {
    long lo = l*tidx;
    long loclen = l;
    if (tidx==(SAL_THREADS-1)) loclen = ll; // last sub-problem special case
    long i = lo;
    long j = lo + loclen - 1;

    // write into pbuffer at i and j
    #pragma acc loop seq
    while (i < j) {
      #pragma acc loop seq
      while (!particles[i]._absorbed && i<j) {
        pbuffer[i] = particles[i];
        i++;
      }
      #pragma acc loop seq
      while (particles[j]._absorbed && i<j) {
        pbuffer[j] = particles[j];
        j--;
      }
      if (i < j) {
        pbuffer[j] = particles[i];
        pbuffer[i] = particles[j];
        i++;
        j--;
      }
    }
    // transfer edge case
    if (i==j)
      pbuffer[i] = particles[i];

    lens[tidx] = i - lo;
    if (i==j && !particles[i]._absorbed) lens[tidx]++;
  }

  // determine lo's
  long accumlen = 0;
  #pragma acc loop seq
  for (long idx=0; idx<SAL_THREADS; idx++) {
    los[idx] = accumlen;
    accumlen = accumlen + lens[idx];
  }

  // step 2: write non-absorbed sub-arrays to psorted/output from the left
  #pragma acc parallel loop present(pbuffer)
  for (unsigned long tidx=0; tidx<SAL_THREADS; tidx++) {
    long j, k;
    #pragma acc loop seq
    for (long i=0; i<lens[tidx]; i++) {
      j = i + l*tidx;
      k = i + los[tidx];
      particles[k] = pbuffer[j];
    }
  }
  //for (int ii=0;ii<accumlen;ii++) printf("%ld ", (psorted[ii]->_absorbed));

  // return (no SPLIT)
  if (flag_split != 1)
    return accumlen;

  // SPLIT - repeat the non-absorbed block N-1 times, where len % accumlen = N + R
  int mult = buffer_len / accumlen; // TODO: possibly use a new arg, bufferlen, rather than len

  // not enough space for full-block split, return
  if (mult <= 1)
    return accumlen;

  // copy non-absorbed block
  #pragma acc parallel loop present(particles)
  for (long tidx = 0; tidx < accumlen; tidx++) { // tidx: thread index
    unsigned long randstate[7];
    _class_particle sourcebuffer;
    _class_particle targetbuffer;
    // assign reduced weight to all particles
    particles[tidx].p=particles[tidx].p/mult;
    #pragma acc loop seq
    for (long bidx = 1; bidx < mult; bidx++) { // bidx: block index
      // preserve absorbed particle (for randstate)
      sourcebuffer = particles[bidx*accumlen + tidx];
      // buffer full particle struct
      targetbuffer = particles[tidx];
      // reassign previous randstate
      targetbuffer.randstate[0] = sourcebuffer.randstate[0];
      targetbuffer.randstate[1] = sourcebuffer.randstate[1];
      targetbuffer.randstate[2] = sourcebuffer.randstate[2];
      targetbuffer.randstate[3] = sourcebuffer.randstate[3];
      targetbuffer.randstate[4] = sourcebuffer.randstate[4];
      targetbuffer.randstate[5] = sourcebuffer.randstate[5];
      targetbuffer.randstate[6] = sourcebuffer.randstate[6];
      // apply
      particles[bidx*accumlen + tidx] = targetbuffer;
    }
  }

  // set out split multiplier value
  *multiplier = mult;

  // return expanded array size
  return accumlen * mult;
}

#endif

/*
*  Fallback serial version of the one above.
*/
long sort_absorb_last_serial(_class_particle* particles, long len) {
  long i = 0;
  long j = len - 1;
  _class_particle pbuffer;

  // bubble
  while (i < j) {
    while (!particles[i]._absorbed && i<j) i++;
    while (particles[j]._absorbed && i<j) j--;
    if (i < j) {
      pbuffer = particles[j];
      particles[j] = particles[i];
      particles[i] = pbuffer;
      i++;
      j--;
    }
  }

  // return new length
  if (i==j && !particles[i]._absorbed)
    return i + 1;
  else
    return i;
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void mccoordschange(Coords a, Rotation t, _class_particle *particle)
{
  Coords b, c;

  b.x = particle->x;
  b.y = particle->y;
  b.z = particle->z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  particle->x = b.x;
  particle->y = b.y;
  particle->z = b.z;

#if MCCODE_PARTICLE_CODE == 2112
    if (particle->vz != 0.0 || particle->vx != 0.0 || particle->vy != 0.0)
      mccoordschange_polarisation(t, &(particle->vx), &(particle->vy), &(particle->vz));

    if (particle->sz != 0.0 || particle->sx != 0.0 || particle->sy != 0.0)
      mccoordschange_polarisation(t, &(particle->sx), &(particle->sy), &(particle->sz));
#elif MCCODE_PARTICLE_CODE == 22
    if (particle->kz != 0.0 || particle->kx != 0.0 || particle->ky != 0.0)
      mccoordschange_polarisation(t, &(particle->kx), &(particle->ky), &(particle->kz));

    if (particle->Ez != 0.0 || particle->Ex != 0.0 || particle->Ey != 0.0)
      mccoordschange_polarisation(t, &(particle->Ex), &(particle->Ey), &(particle->Ez));
#endif
}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
void normal_vec(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order_old(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

int solve_2nd_order(double *t0, double *t1, double A, double B, double C){
  int retval=0;
  double sign=copysign(1.0,B);
  double dt0,dt1;

  dt0=0;
  dt1=0;
  if(t1){ *t1=0;}

  /*protect against rounding errors by locally equating DBL_EPSILON with 0*/
  if (fabs(A)<DBL_EPSILON){
    A=0;
  }
  if (fabs(B)<DBL_EPSILON){
    B=0;
  }
  if (fabs(C)<DBL_EPSILON){
    C=0;
  }

  /*check if coefficient are sane*/
  if( A==0  && B==0){
    retval=0;
  }else{
    if(A==0){
      /*equation is linear*/
      dt0=-C/B;
      retval=1;
    }else if (C==0){
      /*one root is 0*/
      if(sign<0){
        dt0=0;dt1=-B/A;
      }else{
        dt0=-B/A;dt1=0;
      }
      retval=2;
    }else{
      /*a regular 2nd order eq. Also works out fine for B==0.*/
      double D;
      D=B*B-4*A*C;
      if (D>=0){
        dt0=(-B - sign*sqrt(B*B-4*A*C))/(2*A);
        dt1=C/(A*dt0);
        retval=2;
      }else{
        /*no real roots*/
        retval=0;
      }
    }
    /*sort the solutions*/
    if (retval==1){
      /*put both solutions in t0 and t1*/
      *t0=dt0;
      if(t1) *t1=dt1;
    }else{
      /*we have two solutions*/
      /*swap if both are positive and t1 smaller than t0 or t1 the only positive*/
      int swap=0;
      if(dt1>0 && ( dt1<dt0 || dt0<=0) ){
        swap=1;
      }
      if (swap){
        *t0=dt1;
        if(t1) *t1=dt0;
      }else{
        *t0=dt0;
        if(t1) *t1=dt0;
      }
    }

  }
  return retval;

} /*solve_2nd_order_improved*/


/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void _randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double radius,
        _class_particle* _particle)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos(1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
}
/* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void _randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi, double width, double height, Rotation A,
        _class_particle* _particle)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle, xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);
}
/* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/
void _randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
        double xi, double yi, double zi,
        double width, double height, Rotation A,
        double lx, double ly, double lz, int order,
        _class_particle* _particle)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
        *solid_angle = *solid_angle * cos_theta;
      }
    }
  }
}
/* randvec_target_rect_real */


/* SECTION: random numbers ==================================================

  How to add a new RNG:

  - Use an rng with a manegable state vector, e.g. of lengt 4 or 7. The state
  will sit on the particle struct as a "randstate_t state[RANDSTATE_LEN]"
  - If the rng has a long state (as MT), set an empty "srandom" and initialize
  it explicitly using the appropriate define (RNG_ALG)
  - Add a seed and a random function (the transforms will be reused)
  - Write the proper defines in mccode-r.h, e.g. randstate_t and RANDSTATE_LEN,
  srandom and random.
  - Compile using -DRNG_ALG=<selector int value>

============================================================================= */


/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */
/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/
#include <stdio.h>
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

unsigned long mt[N]; /* the array for the state vector  */
int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

// required for common rng alg interface (see RNG_ALG usage in mccode-r.h)
void mt_srandom_empty() {}

// initializes mt[N] with a seed
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
/* Initialize by an array with array-length.
   Init_key is the array for initializing keys.
   key_length is its length. */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}
/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
/* End of "Mersenne Twister". */


/*
KISS

 From: http://www.helsbreth.org/random/rng_kiss.html
 Scott Nelson 1999

 Based on Marsaglia's KISS or (KISS+SWB) <http://www.cs.yorku.ca/~oz/marsaglia-
rng.html>

 KISS - Keep it Simple Stupid PRNG

 the idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
        x(n)=a*x(n-1)+1 mod 2^32
        y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
        z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
      2^32*(2^32-1)*(2^63+2^32-1) > 2^127
*/

/* the KISS state is stored as a vector of 7 unsigned long        */
/*   0  1  2  3  4      5  6   */
/* [ x, y, z, w, carry, k, m ] */

unsigned long *kiss_srandom(unsigned long state[7], unsigned long seed) {
  if (seed == 0) seed = 1;
  state[0] = seed | 1; // x
  state[1] = seed | 2; // y
  state[2] = seed | 4; // z
  state[3] = seed | 8; // w
  state[4] = 0;        // carry
  return 0;
}

unsigned long kiss_random(unsigned long state[7]) {
    state[0] = state[0] * 69069 + 1;
    state[1] ^= state[1] << 13;
    state[1] ^= state[1] >> 17;
    state[1] ^= state[1] << 5;
    state[5] = (state[2] >> 2) + (state[3] >> 3) + (state[4] >> 2);
    state[6] = state[3] + state[3] + state[2] + state[4];
    state[2] = state[3];
    state[3] = state[6];
    state[4] = state[5] >> 30;
    return state[0] + state[1] + state[3];
}
/* end of "KISS" rng */


/* FAST KISS in another implementation (Hundt) */

//////////////////////////////////////////////////////////////////////////////
// fast keep it simple stupid generator
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Thomas Mueller hash for initialization of rngs
// http://stackoverflow.com/questions/664014/
//        what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
//////////////////////////////////////////////////////////////////////////////
randstate_t _hash(randstate_t x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x);
  return x;
}


// SECTION: random number transforms ==========================================



// generate a random number from normal law
double _randnorm(randstate_t* state)
{
  static double v1, v2, s; /* removing static breaks comparison with McStas <= 2.5 */
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = _rand01(state);
      u2 = _rand01(state);
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}
// another one
double _randnorm2(randstate_t* state) {
  double x, y, r;
  do {
      x = 2.0 * _rand01(state) - 1.0;
      y = 2.0 * _rand01(state) - 1.0;
      r = x*x + y*y;
  } while (r == 0.0 || r >= 1.0);
  return x * sqrt((-2.0 * log(r)) / r);
}

// Generate a random number from -1 to 1 with triangle distribution
double _randtriangle(randstate_t* state) {
	double randnum = _rand01(state);
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}
double _rand01(randstate_t* state) {
	double randnum;
	randnum = (double) _random();
  // TODO: can we mult instead of div?
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}
// Return a random number between 1 and -1
double _randpm1(randstate_t* state) {
	double randnum;
	randnum = (double) _random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}
// Return a random number between 0 and max.
double _rand0max(double max, randstate_t* state) {
	double randnum;
	randnum = (double) _random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}
// Return a random number between min and max.
double _randminmax(double min, double max, randstate_t* state) {
	return _rand0max(max - min, state) + max;
}


/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", instrument_name, instrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of particles to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --source                   Show the instrument code which was compiled.\n"
#ifdef OPENACC
"\n"
"  --vecsize                  OpenACC vector-size (default: 128)\n"
"  --numgangs                 Number of OpenACC gangs (default: 7813)\n"
"  --gpu_innerloop            Maximum rays to process pr. OpenACC \n"
"                             kernel run (default: 2147483647)\n"
"\n"
#endif
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
#ifdef OPENACC
  fprintf(stderr,
  "This instrument has been compiled with NVIDIA GPU support through OpenACC.\n  Running on systems without such devices will lead to segfaults.\nFurter, fprintf, sprintf and printf have been removed from any component TRACE.\n");
#endif

  if(numipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < numipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(traceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    instrument_name, instrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to numipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((numipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < numipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]) || !strcmp("--version", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]) || !strcmp("--verbose", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strncmp("--vecsize=", argv[i], 10)) {
      vecsize=atoi(&argv[i][10]);
    }    
    else if(!strcmp("--vecsize", argv[i]) && (i + 1) < argc) {
      vecsize=atoi(argv[++i]);
    }
    else if(!strncmp("--numgangs=", argv[i], 11)) {
      numgangs=atoi(&argv[i][11]);
    }
    else if(!strcmp("--numgangs", argv[i]) && (i + 1) < argc) {
      numgangs=atoi(argv[++i]);
    }
    else if(!strncmp("--gpu_innerloop=", argv[i], 16)) {
      gpu_innerloop=(long)strtod(&argv[i][16], NULL);
    }
    else if(!strcmp("--gpu_innerloop", argv[i]) && (i + 1) < argc) {
      gpu_innerloop=(long)strtod(argv[++i], NULL);
    }

    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(!strcmp("--source", argv[i])) {
      printf("/* Source code %s from %s: */\n"
        "/******************************************************************************/\n"
        "%s\n"
        "/******************************************************************************/\n"
        "/* End of source code %s from %s */\n",
        instrument_name, instrument_source, instrument_code,
        instrument_name, instrument_source);
      exit(1);
    }
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < numipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == numipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < numipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", instrument_name, instrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    save(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    finally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  init();

  /* *** parse options *** */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

/* embedding file "mcstas-r.c" */

/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables
*******************************************************************************/
_class_particle mcsetstate(double x, double y, double z, double vx, double vy, double vz,
			   double t, double sx, double sy, double sz, double p, int mcgravitation, void *mcMagnet, int mcallowbackprop)
{
  _class_particle mcneutron;

  mcneutron.x  = x;
  mcneutron.y  = y;
  mcneutron.z  = z;
  mcneutron.vx = vx;
  mcneutron.vy = vy;
  mcneutron.vz = vz;
  mcneutron.t  = t;
  mcneutron.sx = sx;
  mcneutron.sy = sy;
  mcneutron.sz = sz;
  mcneutron.p  = p;
  mcneutron.mcgravitation = mcgravitation;
  mcneutron.mcMagnet = mcMagnet;
  mcneutron.allow_backprop = mcallowbackprop;
  mcneutron._uid       = 0;
  mcneutron._index     = 1;
  mcneutron._absorbed  = 0;
  mcneutron._restore   = 0;
  mcneutron._scattered = 0;

  return(mcneutron);
} /* mcsetstate */

/*******************************************************************************
* mcgetstate: get neutron parameters from particle structure
*******************************************************************************/
_class_particle mcgetstate(_class_particle mcneutron, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
  *x  =  mcneutron.x;
  *y  =  mcneutron.y;
  *z  =  mcneutron.z;
  *vx =  mcneutron.vx;
  *vy =  mcneutron.vy;
  *vz =  mcneutron.vz;
  *t  =  mcneutron.t;
  *sx =  mcneutron.sx;
  *sy =  mcneutron.sy;
  *sz =  mcneutron.sz;
  *p  =  mcneutron.p;

  return(mcneutron);
} /* mcgetstate */


/*******************************************************************************
* mcgenstate: set default neutron parameters
*******************************************************************************/
// Moved to generated code
/* #pragma acc routine seq */
/* _class_particle mcgenstate(void) */
/* { */
/*   return(mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, mcgravitation, mcMagnet, mcallowbackprop)); */
/* } */

/*******************************************************************************
* mccoordschanges: old style rotation routine rot -> (x y z) ,(vx vy vz),(sx,sy,sz)
*******************************************************************************/
void
mccoordschanges(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) )
    mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) )
    mccoordschange_polarisation(t, sx, sy, sz);

}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight)
* return 0 if outside and 1 if inside
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999.
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98
  *******************************************************************************/
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1
 *******************************************************************************/
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */


/* *****************************************************************************
* Start of instrument 'template_body_simple' generated code
***************************************************************************** */

#ifdef MC_TRACE_ENABLED
int traceenabled = 1;
#else
int traceenabled = 0;
#endif
#define MCSTAS "/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//"
int   defaultmain         = 1;
char  instrument_name[]   = "template_body_simple";
char  instrument_source[] = "phonon_eigenvector.instr";
char *instrument_exe      = NULL; /* will be set to argv[0] in main */
char  instrument_code[]   = "Instrument template_body_simple source code phonon_eigenvector.instr is not embedded in this executable.\n  Use --source option when running McStas.\n";

int main(int argc, char *argv[]){return mccode_main(argc, argv);}

/* *****************************************************************************
* instrument 'template_body_simple' and components DECLARE
***************************************************************************** */

/* Instrument parameters: structure and a table for the initialisation
   (Used in e.g. inputparse and I/O function (e.g. detector_out) */

struct _struct_instrument_parameters {
  MCNUM E;
  MCNUM Ef;
  MCNUM Dlambda;
  MCNUM h;
  MCNUM l;
  MCNUM dA3;
  MCNUM Temp;
  MCNUM width;
  MCNUM coll;
  long phononmode;
  long E_steps_high;
  long E_steps_low;
  long Verbose;
  long DISP;
};
typedef struct _struct_instrument_parameters _class_instrument_parameters;

/* instrument SPLIT and GROUP control logic */
struct instrument_logic_struct {
  long Split_phonon_bvk_pg; /* this is the SPLIT counter decremented down to 0 */
  _class_particle Split_phonon_bvk_pg_particle; /* this is the particle to duplicate */
};

struct _instrument_struct {
  char   _name[256]; /* the name of this instrument e.g. 'template_body_simple' */
/* Counters per component instance */
  double counter_AbsorbProp[23]; /* absorbed events in PROP routines */
  double counter_N[23], counter_P[23], counter_P2[23]; /* event counters after each component instance */
  _class_particle _trajectory[23]; /* current trajectory for STORE/RESTORE */
/* Components position table (absolute and relative coords) */
  Coords _position_relative[23]; /* positions of all components */
  Coords _position_absolute[23];
  _class_instrument_parameters _parameters; /* instrument parameters */
  struct instrument_logic_struct logic; /* instrument logic */
} _instrument_var;
struct _instrument_struct *instrument = & _instrument_var;
#pragma acc declare create ( _instrument_var )
#pragma acc declare create ( instrument )

int numipar = 14;
struct mcinputtable_struct mcinputtable[] = {
  "E", &(_instrument_var._parameters.E), instr_type_double, "5", 
  "Ef", &(_instrument_var._parameters.Ef), instr_type_double, "14.7", 
  "Dlambda", &(_instrument_var._parameters.Dlambda), instr_type_double, "0.05", 
  "h", &(_instrument_var._parameters.h), instr_type_double, "0", 
  "l", &(_instrument_var._parameters.l), instr_type_double, "4.2", 
  "dA3", &(_instrument_var._parameters.dA3), instr_type_double, "-90", 
  "Temp", &(_instrument_var._parameters.Temp), instr_type_double, "2", 
  "width", &(_instrument_var._parameters.width), instr_type_double, "0.005", 
  "coll", &(_instrument_var._parameters.coll), instr_type_double, "40", 
  "phononmode", &(_instrument_var._parameters.phononmode), instr_type_int, "0", 
  "E_steps_high", &(_instrument_var._parameters.E_steps_high), instr_type_int, "50", 
  "E_steps_low", &(_instrument_var._parameters.E_steps_low), instr_type_int, "50", 
  "Verbose", &(_instrument_var._parameters.Verbose), instr_type_int, "0", 
  "DISP", &(_instrument_var._parameters.DISP), instr_type_int, "0", 
  NULL, NULL, instr_type_double, ""
};


/* ************************************************************************** */
/*             SHARE user declarations for all components                     */
/* ************************************************************************** */

/* Shared user declarations for all components types 'Source_Maxwell_3'. */
/* A normalised Maxwellian distribution : Integral over all l = 1 */
#pragma acc routine seq
double SM3_Maxwell(double l, double temp)
  {
    double a=949.0/temp;
    return 2*a*a*exp(-a/(l*l))/(l*l*l*l*l);
  }

/* Shared user declarations for all components types 'Monochromator_flat'. */
#ifndef GAUSS
/* Define these arrays only once for all instances. */
/* Values for Gauss quadrature. Taken from Brice Carnahan, H. A. Luther and
James O Wilkes, "Applied numerical methods", Wiley, 1969, page 103.
This reference is available from the Copenhagen UB2 library */
double Gauss_X[] = {-0.987992518020485, -0.937273392400706, -0.848206583410427,
-0.724417731360170, -0.570972172608539, -0.394151347077563,
-0.201194093997435, 0, 0.201194093997435,
0.394151347077563, 0.570972172608539, 0.724417731360170,
0.848206583410427, 0.937273392400706, 0.987992518020485};
double Gauss_W[] = {0.030753241996117, 0.070366047488108, 0.107159220467172,
0.139570677926154, 0.166269205816994, 0.186161000115562,
0.198431485327111, 0.202578241925561, 0.198431485327111,
0.186161000115562, 0.166269205816994, 0.139570677926154,
0.107159220467172, 0.070366047488108, 0.030753241996117};
#pragma acc declare create ( Gauss_X )
#pragma acc declare create ( Gauss_W )

#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2*(rms)*(rms)))/(sqrt(2*PI)*(rms)))
#endif

/* Shared user declarations for all components types 'Phonon_BvK_PG_eigenvector'. */
#ifndef PHONON_SIMPLE
#define PHONON_SIMPLE $Revision$
#define T2E (1/11.605)   /* Kelvin to meV */
#define SMALL_NUMBER 1E-6
#define V_HIGH 8000   // Highest velocity to search` (TO BE CHANGED)
#define MATRIX_TEST 0 // Change to 1 to write matrices to disk
#define TIME_EIGENSYSTEMS 1
#pragma acc routine 


// CONSTANTS
#define pi    3.14159265358979323846264338327950288419716939937510L
#define sqrt3 1.73205080756887729352744634150587236694280525381038L  
#define DIM 12  
    
#define DV 0.001   // Velocity change used for numerical derivative 

    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "/usr/include/linux/types.h"// Change from Mads: On OS X this path is not valid
#include"cnr.h"
//#include <malloc.h>
#include <stdlib.h> 

#include "../src/eigen.h"
   
    extern FILE *matrix_test_file, *parms_test_file;
    extern int matrix_test_count;
    int mode,verbose;


double nbose(double omega, double T)  /* Other name ?? */
  {
    double nb;

    nb= (omega>0) ? 1+1/(exp(omega/(T*T2E))-1) : 1/(exp(-omega/(T*T2E))-1);
    return nb;
  }
#undef T2E
/* Routine types from Numerical Recipies book */
#define UNUSED (-1.11e30)
#define MAXRIDD 60

void fatalerror_cpu(char *s)
  {
    fprintf(stderr,"%s \n",s);
    exit(1);
  }
 
#pragma acc routine 
void fatalerror(char *s)
  {
  #ifndef OPENACC	
    fatalerror_cpu(s);
  #endif
  }

  #pragma acc routine

void bubblesort (int size , double* inputarray , int* index)
//this function modifies the inputarray by bubble sorting, and also modifies the index array as the
//index after the sorting. The index array should be initialized as index[] = {0,1,2,3,4,5,6,7...}
{
    int tempindex = 0;
    int i;
    int j;
    for (i = 0 ; i < size-1 ;  i++)
    {
        for (j = 0 ; j < size-1 ; j++)
        {
            if (inputarray[index[j]] > inputarray[index[j+1]])
            {
                tempindex = index[j];
                index[j] = index[j+1];
                index[j+1] = tempindex;
            }
        }
     }
    return;
}



/* TODO: Det er meget langsomt at bruge pointere. Erstat med flad memory. */
//3*3 matrices multiplication
void matrixmultiplication( double matrix1[3][3], double matrix2[3][3], double result[3][3])
    {
        // multiplies matrix1 by matrix2 and places the result in result
//      double matrix1[3][3]={{*(array1), *(array1+1), *(array1+2)}, {*(array1+3), *(array1+4), *(array1+5)}, {*(array1+6), *(array1+7), *(array1+8)}};
//      double matrix2[3][3]={{*(array2), *(array2+1), *(array2+2)}, {*(array2+3), *(array2+4), *(array2+5)}, {*(array2+6), *(array2+7), *(array2+8)}};

      result[0][0]=matrix2[0][0]*matrix1[0][0]+matrix2[1][0]*matrix1[0][1]+matrix2[2][0]*matrix1[0][2];
      result[0][1]=matrix2[0][1]*matrix1[0][0]+matrix2[1][1]*matrix1[0][1]+matrix2[2][1]*matrix1[0][2];
      result[0][2]=matrix2[0][2]*matrix1[0][0]+matrix2[1][2]*matrix1[0][1]+matrix2[2][2]*matrix1[0][2];
      result[1][0]=matrix2[0][0]*matrix1[1][0]+matrix2[1][0]*matrix1[1][1]+matrix2[2][0]*matrix1[1][2];
      result[1][1]=matrix2[0][1]*matrix1[1][0]+matrix2[1][1]*matrix1[1][1]+matrix2[2][1]*matrix1[1][2];
      result[1][2]=matrix2[0][2]*matrix1[1][0]+matrix2[1][2]*matrix1[1][1]+matrix2[2][2]*matrix1[1][2];
      result[2][0]=matrix2[0][0]*matrix1[2][0]+matrix2[1][0]*matrix1[2][1]+matrix2[2][0]*matrix1[2][2];
      result[2][1]=matrix2[0][1]*matrix1[2][0]+matrix2[1][1]*matrix1[2][1]+matrix2[2][1]*matrix1[2][2];
      result[2][2]=matrix2[0][2]*matrix1[2][0]+matrix2[1][2]*matrix1[2][1]+matrix2[2][2]*matrix1[2][2];

      return;
    }

// TODO: Benchmark mod pointer-snask versioner: matrix3mupltiplication of matrixmultiplication
void matrix3x3_multiply3(const double A[9], const double B[9], const double C[9], double result[9])
{
  // R_{i,j} =  \sum_{k,l=1}^3 A_{i,k}*B_{k,l}*C_{l,j}
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      double ij_sum = 0;
      for(int k=0;k<3;k++)
	for(int l=0;l<3;l++)
	  ij_sum += A[i*3+k]*B[k*3+l]*C[l*3+j];
      result[i*3+j] = ij_sum;
    }
}

/* TODO: Det er meget langsomt at bruge pointere. Erstat med flad memory. */ 
void matrix3multiplication(double matrix1[3][3], double matrix2[3][3], double matrix3[3][3], double result[3][3])
    {
        // Performs the multiplication (matrix1*matrix2)*matrix3 and placed the result in result
//#define TEST_MULTIPLY
        double temp[3][3];
        matrixmultiplication(matrix1,matrix2,temp);
        matrixmultiplication(temp,matrix3,result);
#ifdef TEST_MULTIPLY
        printf("Matrix3multiply called with \n");
        printf("matrix1 = (");
        for (int i=0; i<3; i++)
        {
            printf("( %g %g %g), ",matrix1[i][0], matrix1[i][1], matrix1[i][2]);
        }
        printf(")\n Matrix2 = (");
        for (int i=0; i<3; i++)
        {
            printf("( %g %g %g), ",matrix2[i][0], matrix2[i][1], matrix2[i][2]);
        }
        printf(")\n Matrix3 = (");
                for (int i=0; i<3; i++)
        {
            printf("( %g %g %g), ",matrix3[i][0], matrix3[i][1], matrix3[i][2]);
        }
        printf(")\n Temp = (");
                for (int i=0; i<3; i++)
        {
            printf("( %g %g %g), ",temp[i][0], temp[i][1], temp[i][2]);
        }
        printf(")\n Result = (");
                 for (int i=0; i<3; i++)
        {
            printf("( %g %g %g), ",result[i][0], result[i][1], result[i][2]);
        }
        printf(")\n");
#endif // TEST_MULTIPLY

        return;
    }
    

/* TODO: Erstat med standard complex.h operation */
komplex complexexp (double* q, double* rj)
{
//    double qrjtest
    double qrj = *(q)**(rj)+*(q+1)**(rj+1)+*(q+2)**(rj+2);
//    printf(" q= (%g %g %g), r = (%g %g %g) qrj = %g \n",q[0],q[1],q[2],rj[0],rj[1],rj[2],qrj);
    return cos(qrj)+I*sin(qrj);
}

    // _________________ INITIALIZE VARIABLES omega_q ___________________________
#define a_latt 2.461 //PG lattice constant in AA  paper gives 2.45
#define c_latt 6.708 // PG lattice constant in AA  paper gives 6.70
#define b_length 6.646 // PG scattering legth in fm
#define THz2meV 4.13566
#define Da2kg 1.6605e-27 //Dalton to kg converter
#define dyn2N 1e-3 //Convert dyn/cm to N/m
 
    double M = 12.011 * Da2kg; //atomic mass of C

#define K_l1  3.62e+5      // Force constant longitudinal nn1 - in (a,b) plane; here units of dyn
#define K_t1  1.99e+5      // Force constant transverse nn1 - in (a,b) plane
#define K_l2  1.33e+5      // Force constant longitudinal nn2 - in (a,b) plane
#define K_t2 -0.520e+5     // Force constant transverse nn2 - in (a,b) plane
#define K_l3 -0.037e+5     // Force constant longitudinal nn3 - in (a,b) plane
#define K_t3  0.288e+5     // Force constant transverse nn3 - in (a,b) plane
#define K_l4  0.058e+5     // Force constant longitudinal nn4 - along c
#define K_t4  0.0077e+5    // Force constant transverse nn4 - along c

// #define TEST_MATRICES

    
    //Lattice basis
    
    // PG is hexagonal close packed with 4 atoms per cell
    // Coordinates from the paper Fig. 7: Atoms A, B, C, D
    double Delta[4*3] = {0 , 0 , 0 , a_latt/(2*sqrt3) , a_latt/2 , 0 , 0 , 0 , c_latt/2 ,  -a_latt/(2*sqrt3) , a_latt/2 , c_latt/2};

    // Lattice vectors from paper fig. 7
    double avec[3] = {a_latt*sqrt3/2 , a_latt*0.5 , 0}; // THIS IS NEVER USED ??!!
    double bvec[3] = {a_latt*(-sqrt3/2) , a_latt*0.5 , 0};
    double cvec[3] = {0 , 0 , c_latt};

    // reciprocal lattice vectors
    double astar = 4*PI/(sqrt3*a_latt);  // length of reciprocal lattice vector a*
    double cstar = 2*PI/c_latt;  // length of reciprocal lattice vector c*

    // Rotation matrices
    double Rot120[3][3] = {{-0.5 , sqrt3/2 , 0} , {-sqrt3/2 , -0.5 , 0} , {0 , 0 , 1}};
    double Rot60[3][3] = {{0.5 , sqrt3/2 , 0} , {-sqrt3/2 , 0.5 , 0} , {0 , 0 , 1}};
    double Rot120rev[3][3] = {{-0.5 , -sqrt3/2 , 0} , {sqrt3/2 , -0.5 , 0} , {0 , 0 , 1}}; //rotate -120
    double Rot60rev[3][3] = {{0.5 , -sqrt3/2 , 0} , {sqrt3/2 , 0.5 , 0} , {0 , 0 , 1}}; //rotate -60
    double Rot180[3][3] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}; // rotate 180 or -180
//    double Rot5[3][3] = {{0.5 , sqrt3/2 , 0} , {-sqrt3/2 , 0.5 , 0} , {0 , 0 , 1}}; //rotate 60
//    double Rot5rev[3][3] = {{0.5 , -sqrt3/2 , 0} , {sqrt3/2 , 0.5 , 0} , {0 , 0 , 1}}; //rotate -60
    
    //1st neighbour
    
//    double r_j1[3][3] = {{-a_latt/sqrt3 , 0 , 0} , {a_latt*0.5/sqrt3 , -a_latt*0.5 , 0} , {a_latt*0.5/sqrt3 , a_latt*0.5 , 0}}; // This holds from atoms A and D
//    double r_j2[3][3] = {{a_latt/sqrt3 , 0 , 0} , {-a_latt*0.5/sqrt3 , a_latt*0.5 , 0} , {-a_latt*0.5/sqrt3 , -a_latt*0.5 , 0}}; // This holds from atoms B and C

    double Phi_nn1[3][3] = {{K_l1*dyn2N , 0 , 0} , {0 , K_t1*dyn2N , 0} , {0 , 0 , K_t1*dyn2N}};

    double Phi1[3][3] = {{K_l1*dyn2N , 0 , 0} , {0 , K_t1*dyn2N , 0} , {0 , 0 , K_t1*dyn2N}}; //Phi1 = Phi_nn1

    double Phi2[3][3];

    double Phi3[3][3];
    
     
    //2nd neighbour
//    double r_j11[9][3] = {{a_latt*(-1/sqrt3) , 0 , 0} , {a_latt*0.5/sqrt3 , a_latt*(-0.5) , 0} , {a_latt*0.5/sqrt3 , a_latt*0.5 , 0} , {0 , a_latt , 0} , {-a_latt*sqrt3/2 , a_latt/2 , 0} , {-a_latt*sqrt3/2 , -a_latt/2 , 0} , {0 , -a_latt , 0} , {a_latt*sqrt3/2 , -a_latt/2 , 0} , {a_latt*sqrt3/2 , a_latt/2 , 0}};//r_j1 + r_add
//    double r_j21[9][3] = {{-a_latt*(-1/sqrt3) , -a_latt*0 , 0} , {-a_latt*0.5/sqrt3 , -a_latt*(-0.5) , 0} , {-a_latt*0.5/sqrt3 , -a_latt*0.5 , 0} , {0 , a_latt , 0} , {-a_latt*sqrt3/2 , a_latt/2 , 0} , {-a_latt*sqrt3/2 , -a_latt/2 , 0} , {0 , -a_latt , 0} , {a_latt*sqrt3/2 , -a_latt/2 , 0} , {a_latt*sqrt3/2 , a_latt/2 , 0}};//r_j2 + r_add
    
    double Phi_nn2[3][3] = {{K_t2*dyn2N , 0 , 0} , {0 , K_l2*dyn2N , 0} , {0 , 0 , K_t2*dyn2N}};
    double Phi4[3][3] = {{K_t2*dyn2N , 0 , 0} , {0 , K_l2*dyn2N , 0} , {0 , 0 , K_t2*dyn2N}};
    double Phi5[3][3];
    
    double Phi6[3][3];
        
    double Phi7[3][3];

    double Phi8[3][3];
     
    double Phi9[3][3];

    //3rd neighbour
//    double r_j12[12][3] = {{a_latt*(-1/sqrt3) , 0 , 0} , {a_latt*0.5/sqrt3 , a_latt*(-0.5) , 0} , {a_latt*0.5/sqrt3 , a_latt*0.5 , 0} , {0 , a_latt , 0} , {-a_latt*sqrt3/2 , a_latt/2 , 0} , {-a_latt*sqrt3/2 , -a_latt/2 , 0} , {0 , -a_latt , 0} , {a_latt*sqrt3/2 , -a_latt/2 , 0} , {a_latt*sqrt3/2 , a_latt/2 , 0} , {2*a_latt/sqrt3 , 0 , 0} , {-a_latt/sqrt3 , a_latt , 0} , {-a_latt/sqrt3 , -a_latt , 0}};//r_j11 + r_add2
    double r_j2[12][3] = {{a_latt/sqrt3 , 0 , 0} , {-a_latt*0.5/sqrt3 , a_latt*0.5 , 0} , {-a_latt*0.5/sqrt3 , -a_latt*0.5 , 0} , {0 , a_latt , 0} , {-a_latt*sqrt3/2 , a_latt/2 , 0} , {-a_latt*sqrt3/2 , -a_latt/2 , 0} , {0 , -a_latt , 0} , {a_latt*sqrt3/2 , -a_latt/2 , 0} , {a_latt*sqrt3/2 , a_latt/2 , 0} , {-2*a_latt/sqrt3 , 0 , 0} , {a_latt/sqrt3 , -a_latt , 0} , {a_latt/sqrt3 , a_latt , 0}};//r_j21 - r_add2

     
    double Phi_nn3[3][3] = {{K_l3*dyn2N , 0 , 0} , {0 , K_t3*dyn2N , 0} , {0 , 0 , K_t3*dyn2N}};

    // c2 = 2/sqrt(6);
    // c1 = 1/sqrt(6);
    // c0 = 1/sqrt(2);
    // c3 = 1/sqrt(3);

    double Phi10[3][3] = {{K_l3*dyn2N , 0 , 0} , {0 , K_t3*dyn2N , 0} , {0 , 0 , K_t3*dyn2N}}; //Phi_nn3;

    double Phi11[3][3];
    
    double Phi12[3][3];
    
    //4th neighbour
    double r_j13[14][3] = {{-a_latt/sqrt3 , 0 , 0} , {a_latt*0.5/sqrt3 , -a_latt*0.5 , 0} , {a_latt*0.5/sqrt3 , a_latt*0.5 , 0} , {0 , a_latt , 0} , {-a_latt*sqrt3/2 , a_latt/2 , 0} , {-a_latt*sqrt3/2 , -a_latt/2 , 0} , {0 , -a_latt , 0} , {a_latt*sqrt3/2 , -a_latt/2 , 0} , {a_latt*sqrt3/2 , a_latt/2 , 0} , {2*a_latt/sqrt3 , 0 , 0} , {-a_latt/sqrt3 , a_latt , 0} , {-a_latt/sqrt3 , -a_latt , 0} , {0 , 0 , c_latt/2} , {0 , 0 , -c_latt/2}};//r_j12 + r_add3 Sublattice B and D do not have this coupling

    double Phi_nn4[3][3] = {{K_t4*dyn2N , 0 , 0} , {0 , K_t4*dyn2N , 0} , {0 , 0 , K_l4*dyn2N}};
     
    double Phi13[3][3] = {{K_t4*dyn2N , 0 , 0} , {0 , K_t4*dyn2N , 0} , {0 , 0 , K_l4*dyn2N}}; //Phi_nn4;
    double Phi14[3][3] = {{K_t4*dyn2N , 0 , 0} , {0 , K_t4*dyn2N , 0} , {0 , 0 , K_l4*dyn2N}}; //Phi_nn4;
    
    
    // _______________ END INITIALIZE omega_q ___________________________
 
void initialize_omega_q()
{
    // Make the matrix algebra that finilizes the force constant set-up
    
    matrix3multiplication(Rot120,Phi_nn1,Rot120rev,Phi2); //Rot2*Phi_nn1*Rot2^(-1) 
    matrix3multiplication(Rot120rev,Phi_nn1,Rot120,Phi3);
    matrix3multiplication(Rot60,Phi_nn2,Rot60rev,Phi5);
    matrix3multiplication(Rot120,Phi_nn2,Rot120rev,Phi6); //Rot6*Phi_nn2*Rot6^(-1);
    matrix3multiplication(Rot180,Phi_nn2,Rot180,Phi7);    //Rot7*Phi_nn2*Rot7^(-1);
    matrix3multiplication(Rot120rev,Phi_nn2,Rot120,Phi8); //Rot8*Phi_nn2*Rot8^(-1);
    matrix3multiplication(Rot60rev,Phi_nn2,Rot60,Phi9);    //Rot9*Phi_nn2*Rot9^(-1);
    matrix3multiplication(Rot120,Phi_nn3,Rot120rev,Phi11);  //Rot11*Phi_nn3*Rot11^(-1);
    matrix3multiplication(Rot120rev,Phi_nn3,Rot120,Phi12); //Rot12*Phi_nn3*Rot12^(-1);
    
    return;

}
    
    double omega_q(komplex* parms, komplex* eigenvectorgood) {
    
    // _______________ DETERMINE Q ____________________________
    
      double vi, vf, vv_x, vv_y, vv_z, vi_x, vi_y, vi_z;
      double q, qx, qy, qz, Jq, hw_phonon, hw_neutron;
      double h,k,l;
      int disp;

// #define ASSERT_NONANS_COMPLEX(A,n) for(int ii=0;ii<n;ii++) if(isnan(creal(A[ii]) || isnan(cimag(A[ii])))){ fprintf(stderr,"Array " #A " contains a nan at position %d.\n",ii); abort(); }
      
//       ASSERT_NONANS_COMPLEX(parms,8);

      for(int ii=0;ii<8;ii++) if(isnan(creal(parms[ii])) || isnan(cimag(parms[ii]))){
	  fprintf(stderr,"omega_q: Array parms contains a nan at position %d.\n",ii);
	  fprintf(stderr,"parms = ["); for(int i=0;i<8;i++) fprintf(stderr,"%g + i%g%s",creal(parms[i]), cimag(parms[i]), i+1<8?", ":"]\n");
	  abort();
	}

      //      fprintf(stderr,"omega_q: no NaN's\n");      
      vf=parms[0];
      vi=parms[1];
      vv_x=parms[2];
      vv_y=parms[3];
      vv_z=parms[4];
      vi_x=parms[5];
      vi_y=parms[6];
      vi_z=parms[7];
      h = parms[9];
      k = parms[10];
      l = parms[11];
      disp = (int)parms[12];

 
      qx=V2K*(vi_x-vf*vv_x);
      qy=V2K*(vi_y-vf*vv_y);
      qz=V2K*(vi_z-vf*vv_z);
    
      if (disp==1) /* calculate dispersion directly */
      {
          printf("Dispersion calculation: ");
          qx = h*astar;
          qy = k*astar;
          qz = l*cstar;
      }
      
    double qvec[3]={qx,qy,qz};

    
    if (verbose>=3)
        printf("Enter omega_q; q = (%g, %g, %g) \n",qx,qy,qz);
    
    // ________________ END DETERMINE Q ________________-
    

    /* TODO: Det er meget langsommere med pointere. Erstat med flad memory. [3][3] -> [9] over alt.*/
    /* TODO: Ryd op og streamline! Pænere OG hurtigere. */
    komplex Phi_diag[3][3];
    int i;
    int j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_diag[i][j] = Phi1[i][j]+Phi2[i][j]+Phi3[i][j]+Phi4[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[3][0]))+Phi5[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[4][0]))+Phi6[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[5][0]))+Phi7[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[6][0]))+Phi8[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[7][0]))+Phi9[i][j]*(1-(*complexexp)(&qvec[0], &r_j13[8][0]))+Phi10[i][j]+Phi11[i][j]+Phi12[i][j];
        }
    }
    
    komplex Phi_diag1[3][3];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_diag1[i][j] = Phi13[i][j]+Phi14[i][j];
        }
    }
    
    komplex Phi_offdiag1[3][3];
    komplex Phi_offdiag2[3][3];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_offdiag1[i][j] =-Phi1[i][j]*((*complexexp)(&qvec[0], &r_j13[0][0]))-Phi2[i][j]*((*complexexp)(&qvec[0], &r_j13[1][0]))-Phi3[i][j]*((*complexexp)(&qvec[0], &r_j13[2][0]))-Phi10[i][j]*((*complexexp)(&qvec[0], &r_j13[9][0]))-Phi11[i][j]*((*complexexp)(&qvec[0], &r_j13[10][0]))-Phi12[i][j]*((*complexexp)(&qvec[0], &r_j13[11][0]));
            Phi_offdiag2[i][j] = conj(Phi_offdiag1[i][j]);
        }
    }
    
    komplex Phi_6Dim[6][6];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_6Dim[i][j] = Phi_diag[i][j]+Phi_diag1[i][j];
            Phi_6Dim[i][j+3] = Phi_offdiag1[i][j];
            Phi_6Dim[i+3][j] = Phi_offdiag2[i][j];
            Phi_6Dim[i+3][j+3] = Phi_diag[i][j];  // Because Phi_diag1 does not appear for the B,C, sublattices
        }
    }
    
    komplex Phi_AC[3][3];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_AC[i][j] = -Phi13[i][j]*((*complexexp)(&qvec[0], &r_j13[12][0]))-Phi14[i][j]*((*complexexp)(&qvec[0], &r_j13[13][0]));
        }
    }
    
    komplex Zero_3Dim[3][3] = {{0 , 0 , 0},{0 , 0 , 0},{0 , 0 , 0}};

    komplex Phi_6Dim_offdiag[6][6];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Phi_6Dim_offdiag[i][j] = Phi_AC[i][j];
            Phi_6Dim_offdiag[i][j+3] = Zero_3Dim[i][j];
            Phi_6Dim_offdiag[i+3][j] = Zero_3Dim[i][j];
            Phi_6Dim_offdiag[i+3][j+3] = Zero_3Dim[i][j];
        }
    }

    komplex Phi_12Dim[12][12];
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            Phi_12Dim[i][j] = Phi_6Dim[i][j];
            Phi_12Dim[i][j+6] = Phi_6Dim_offdiag[i][j];
            Phi_12Dim[i+6][j] = conj(Phi_6Dim_offdiag[i][j]);
            Phi_12Dim[i+6][j+6] = Phi_6Dim[i][j];
        }
    }
    
#ifdef TEST_MATRICES
    printf("Phi_12Dim=\n");
    for (i = 0; i < 12; i++) {
      for (j = 0; j< 12; j++) {
	if (j == 0) {
	  printf("{");
	}
	printf("%g+(%g)i  ", real(Phi_12Dim[i][j]) , imag((Phi_12Dim[i][j])));
	if (j == 11) {
	  printf("}\n");
	}
      }
    }
    
    printf("Phi_6Dim=\n");
    for (i = 0; i < 6; i++) {
      for (j = 0; j< 6; j++) {
	if (j == 0) {
	  printf("{");
	}
	printf("%g+(%g)i  ", real(Phi_6Dim[i][j]) , imag((Phi_6Dim[i][j])));
	if (j == 5) {
	  printf("}\n");
	}
      }
    }

    printf("Phi_diag=\n");
    for (i = 0; i < 3; i++) {
      for (j = 0; j< 3; j++) {
	if (j == 0) {
	  printf("{");
	}
	printf("%g+(%g)i  ", real(Phi_diag[i][j]) , imag((Phi_diag[i][j])));
	if (j == 2) {
	  printf("}\n");
	}
      }
    }
    printf("Phi2={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}}\n",Phi2[0][0],Phi2[0][1],Phi2[0][2],Phi2[1][0],Phi2[1][1],Phi2[1][2],Phi2[2][0],Phi2[2][1],Phi2[2][2]);
    printf("Phi3={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi3[0][0],Phi3[0][1],Phi3[0][2],Phi3[1][0],Phi3[1][1],Phi3[1][2],Phi3[2][0],Phi3[2][1],Phi3[2][2]);
    printf("Phi5={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi5[0][0],Phi5[0][1],Phi5[0][2],Phi5[1][0],Phi5[1][1],Phi5[1][2],Phi5[2][0],Phi5[2][1],Phi5[2][2]);
    printf("Phi6={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi6[0][0],Phi6[0][1],Phi6[0][2],Phi6[1][0],Phi6[1][1],Phi6[1][2],Phi6[2][0],Phi6[2][1],Phi6[2][2]);
    printf("Phi7={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi7[0][0],Phi7[0][1],Phi7[0][2],Phi7[1][0],Phi7[1][1],Phi7[1][2],Phi7[2][0],Phi7[2][1],Phi7[2][2]);
    printf("Phi8={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi8[0][0],Phi8[0][1],Phi8[0][2],Phi8[1][0],Phi8[1][1],Phi8[1][2],Phi8[2][0],Phi8[2][1],Phi8[2][2]);
    printf("Phi9={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi9[0][0],Phi9[0][1],Phi9[0][2],Phi9[1][0],Phi9[1][1],Phi9[1][2],Phi9[2][0],Phi9[2][1],Phi9[2][2]);
    printf("Phi11={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi11[0][0],Phi11[0][1],Phi11[0][2],Phi11[1][0],Phi11[1][1],Phi11[1][2],Phi11[2][0],Phi11[2][1],Phi11[2][2]);
    printf("Phi12={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi12[0][0],Phi12[0][1],Phi12[0][2],Phi12[1][0],Phi12[1][1],Phi12[1][2],Phi12[2][0],Phi12[2][1],Phi12[2][2]);
    printf("Phi14={{%lf,%lf,%lf},{%lf,%lf,%lf},{%lf,%lf,%lf}\n}",Phi14[0][0],Phi14[0][1],Phi14[0][2],Phi14[1][0],Phi14[1][1],Phi14[1][2],Phi14[2][0],Phi14[2][1],Phi14[2][2]);

    printf(" 1-exp(q.r_j13[3]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[3][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[3][0])));
    printf(" 1-exp(q.r_j13[4]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[4][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[4][0])));
    printf(" 1-exp(q.r_j13[5]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[5][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[5][0])));
    printf(" 1-exp(q.r_j13[6]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[6][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[6][0])));
    printf(" 1-exp(q.r_j13[7]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[7][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[7][0])));
    printf(" 1-exp(q.r_j13[8]) = %lf + i %lf \n",real(1-(*complexexp)(&qvec[0], &r_j13[8][0])),imag(1-(*complexexp)(&qvec[0], &r_j13[8][0])));
    
#endif // TEST_MATRICES
    
  double eigenvalue[DIM]; 
  komplex Q[DIM*DIM];
  komplex matrix[DIM*DIM];
  int index[DIM];
  
  
  for (i=0; i<DIM; i++)
    for (j=i; j<DIM; j++){
      matrix[i*DIM+j] = Phi_12Dim[i][j];
      matrix[j*DIM+i] = conj(matrix[i*DIM+j]);
    }

#if MATRIX_TEST
    fwrite(matrix,sizeof(komplex),DIM*DIM,matrix_test_file);
    fwrite(parms, sizeof(komplex),8,    parms_test_file);
#endif  
  

 eigensystem_hermitian(matrix, DIM, eigenvalue, Q); /* Produces unitary transformation matrix Q: eigenvectors are columns */

   for (i=0; i<DIM; i++)  index[i]=i;

  if (verbose >= 6)
  {
      printf("Before bubblesort\n");
      for (i=0; i<DIM; i++)
      {
          printf("i = %i eigenvalue[i] = %g index[i] %i \n",i,eigenvalue[i],index[i]);
      }
  }
  // Sort eigenvalues from lowest to highest
  //mergesort(eigenvalue,DIM,index);  

  bubblesort(DIM, eigenvalue, index);
  if (verbose >= 5)
  {
      printf("After bubblesort\n");
      for (i=0; i<DIM; i++)
      {
          printf("i = %i eigenvalue[i] = %g index[i] %i \n",i,eigenvalue[i],index[i]);
      }
  }   

  double eigenvaluegood[DIM];
  
    for (j = 0; j < DIM ; j++){
      eigenvectorgood[j] = Q[index[mode]*DIM + j]; /* Eigenvectors are columns of Q */ // No, they are rows!
    }

      if (verbose>=7)
      {
          printf("Q = [");
          for (i = 0; i < DIM ; i++)
          {
              printf("\n (");
              for (j=0; j<DIM; j++)
                  printf(" %g +i %g, ",real(Q[i*DIM+j]),imag(Q[i*DIM+j]));
              printf(")");
          }
          printf("]\n");
      }

      for (i=0; i<DIM; i++)
        eigenvaluegood[i] = sqrt(real(eigenvalue[index[i]]))/sqrt(M)/(2*PI*1E12)*THz2meV; // convert to units of THz    
	
  

    if (verbose>=5)
    {
        printf("Diagonalization done, eigenvalues (energies) are: \n (");
        for (i=0; i<DIM; i++)
            printf(" %g, ",eigenvaluegood[i]);
        printf(" )\n");
        printf("mode = %i , eigenvectorgood[mode] = (",mode);
        for (j=0; j<DIM; j++)
            printf(" %g +i %g, ",real(eigenvectorgood[j]),imag(eigenvectorgood[j]));
        printf(" )\n");
    }

        
  // ___________________ RETURN RESULTS TO CALLING FUNCTION __________________
  
  hw_neutron = fabs(VS2E*(vi*vi-vf*vf)); // neutron energy transfer 
  // fabs used to find both positive and negative solutions, controlled by findroots()  
      
  hw_phonon = eigenvaluegood[mode];  // phonon energy

  if(verbose >= 4){
    printf("Return from omega_q \n");
    printf("hw_neutron = %g (vi = %g, vf = %g)\n",hw_neutron, vi, vf);
    printf("hw_phonon  = %g = eigenvaluegood[%d] = %g = eigenvalues[%d] = %g\n",
	    hw_phonon,mode,eigenvaluegood[mode],index[mode], eigenvalue[index[mode]]);
//    printf("old parms  = ["); for(int i=0;i<8;i++) fprintf(stderr,"%g+i%g%s",creal(parms[i]), cimag(parms[i]), i+1<8?", ":"]\n");
  }
  
  parms[8] = hw_phonon;
  
  if (disp==1) {
      printf("mode = %d, (h,k,l) = (%g,%g,%g), hw_q = %g",mode,h,k,l,hw_phonon);
      printf(" polarization: (%g + i %g , %g + i %g , %g + i %g) \n",real(eigenvectorgood[0]),imag(eigenvectorgood[0]),real(eigenvectorgood[1]),imag(eigenvectorgood[1]), real(eigenvectorgood[2]),imag(eigenvectorgood[2]));
  }
  
// if (abs(hw_phonon - hw_neutron) < SMALL_NUMBER)
//      printf("q = ( %g , %g , %g ), q_abs = ( %g %g %g),  hw_q = %g \n",qx/astar,qy/astar,qz/cstar,qx,qy,qz,hw_phonon);

  return (hw_phonon - hw_neutron);
}

    // ------------------START OF THE DIAGONALIZATION CODE--------------------

#define SIGN_CONTROLLED(a,b) ((b)<0 ? -fabs(a) : fabs(a))
//#define TEST_DIAGONALIZE
//#define TEST_DIAGONALIZE_LV2
//#define FIND_EIGENSTATE



void nrerror(const char *error_text){
  printf("Numerical Recipes run-time error ...\n");
  printf("%s\n",error_text);
  printf("... now exiting to system.\n");
  exit(1);
}

double last_fh, last_x2;

double zridd(double (*func)(komplex*,komplex*), double x1, double x2, komplex *parms, komplex* vector, double xacc)
    {
      int j;
      double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

      //      fprintf(stderr,"i zridd(%g,%g, parms,%g)\n",x1,x2,xacc);
      parms[0]=x1;
      //      fprintf(stderr,"zridd fl = omega_q\n");      
      if (xacc>0)
      {
          fl=(*func)(parms,vector);
      }
      else
      {
          fl = last_fh;
          xacc = -xacc;
          if (abs(x1-last_x2)>xacc) 
              printf("*** error in zridd *** x1: %g last_x2: %g \n",x1,last_x2);
      }
      //      fprintf(stderr,"fl = %g\n",fl);                  
      parms[0]=x2;
      //      fprintf(stderr,"zridd fh = omega_q; parms[0] = x2 = %g+i%g = %g\n",creal(parms[0]), cimag(parms[0]),x2);
      last_fh=fh=(*func)(parms,vector);
      last_x2 = x2;
      //      fprintf(stderr,"fh = %g\n",fh);            
      if (fl*fh >= 0)
      {
        if (fl==0) return x1;
        if (fh==0) return x2;
        return UNUSED;
      }
      else
      {
        xl=x1;
        xh=x2;
//        printf("zridd sign change: v_low : %g v_high: %g \n",xl,xh);
        ans=UNUSED;
        for (j=1; j<MAXRIDD; j++)
        {
          xm=0.5*(xl+xh);
          parms[0]=xm;
	  //	  fprintf(stderr,"zridd fm = omega_q\n");	  
          fm=(*func)(parms,vector);
          s=sqrt(fm*fm-fl*fh);
          if (s == 0.0)
            return ans;
          xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
          if (fabs(xnew-ans) <= xacc)
            return ans;
          ans=xnew;
          parms[0]=ans;
	  //	  fprintf(stderr,"s = %g;  fl, fm, fh = %g,%g,%g; fnew = %g; fm*fm-fl*fh = %g\n",
	  //		  s, fl, fm,fh, fnew, fm*fm-fl*fh);
	  //fprintf(stderr,"zridd fnew = omega_q\n");
          fnew=(*func)(parms,vector);
          if (fnew == 0.0) return ans;
          if (fabs(fm)*SIGN(fnew) != fm)
          {
            xl=xm;
            fl=fm;
            xh=ans;
            fh=fnew;
          }
          else
            if (fabs(fl)*SIGN(fnew) != fl)
            {
              xh=ans;
              fh=fnew;
            }
            else
              if(fabs(fh)*SIGN(fnew) != fh)
              {
                xl=ans;
                fl=fnew;
              }
              else
                fatalerror("never get here in zridd");
          if (fabs(xh-xl) <= xacc)
            return ans;
        }
        fatalerror("zridd exceeded maximum iterations");
      }
      return 0.0;  /* Never get here */
    }

#pragma acc routine 
double zridd_gpu(double x1, double x2, komplex *parms, komplex* vector, double xacc)
    {
      int j;
      double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

      parms[0]=x1;
      //      fprintf(stderr,"fl=omega_q(%g+i%g)\n",creal(parms[0]),cimag(parms[0]));                  
      fl=omega_q(parms,vector);
      parms[0]=x2;
      //      fprintf(stderr,"fh=omega_q(%g+i%g)\n",creal(parms[0]),cimag(parms[0]));                  
      fh=omega_q(parms,vector);
      if (fl*fh >= 0)
      {
        if (fl==0) return x1;
        if (fh==0) return x2;
        return UNUSED;
      }
      else
      {
        xl=x1;
        xh=x2;
        ans=UNUSED;
        for (j=1; j<MAXRIDD; j++)
        {
          xm=0.5*(xl+xh);
          parms[0]=xm;
	  //	  fprintf(stderr,"fm=omega_q(%g+i%g)\n",creal(parms[0]),cimag(parms[0]));            	  
          fm=omega_q(parms,vector);
          s=sqrt(fm*fm-fl*fh);
          if (s == 0.0)
            return ans;
          xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
          if (fabs(xnew-ans) <= xacc)
            return ans;
          ans=xnew;
          parms[0]=ans;
	  //	  fprintf(stderr,"fnew=omega_q(%g+i%g)\n",creal(parms[0]),cimag(parms[0]));            	  
          fnew=omega_q(parms,vector);
          if (fnew == 0.0) return ans;
          if (fabs(fm)*SIGN(fnew) != fm)
          {
            xl=xm;
            fl=fm;
            xh=ans;
            fh=fnew;
          }
          else
            if (fabs(fl)*SIGN(fnew) != fl)
            {
              xh=ans;
              fh=fnew;
            }
            else
              if(fabs(fh)*SIGN(fnew) != fh)
              {
                xl=ans;
                fl=fnew;
              }
              else
                fatalerror("never get here in zridd");
          if (fabs(xh-xl) <= xacc)
            return ans;
        }
        fatalerror("zridd exceeded maximum iterations");
      }
      return 0.0;  /* Never get here */
    }

 
#define ROOTACC 1e-8

  int findroots(double brack_low, double brack_mid, double brack_high, double *list, int* index, double (*f)(komplex*,komplex*), komplex *parms, komplex* vector, int low_steps, int high_steps)
    {
      double root,range;
      int i;

     if (verbose==2)
          printf("Enter findroots() \n");
     
 // First, find roots in energy loss side, if any
     range = brack_mid-brack_low;
     for (i=0; i<(low_steps-1); i++) 
     {
        if (i==0)
           root = zridd(f, brack_low+range*i/low_steps,
                   brack_low+range*(i+1)/low_steps,
                   (komplex *)parms, (komplex *)vector, ROOTACC);
        else
           root = zridd(f, brack_low+range*i/low_steps,
                   brack_low+range*(i+1)/low_steps,
                   (komplex *)parms, (komplex *)vector, -ROOTACC);  // Re-use the previous fh value
      if (root != UNUSED)
      {
        list[(*index)++]=root;
        if (verbose >=3)
            printf("--- findroots returned a root on the energy loss side: vf = %g \n",root);
      }
     }
   
     range = (brack_mid-brack_low)/(double)low_steps;
     for (i=0; i<low_steps; i++)  // Small steps close to zero energy transfer
     {
         if (low_steps==1)  // then the previous loop did not execute
            root = zridd(f, brack_low+range*(low_steps-1+i/(double)low_steps),
                   brack_low+range*(low_steps-1+(i+1)/(double)low_steps),
                   (komplex *)parms, (komplex *)vector, ROOTACC);
        else
            root = zridd(f, brack_low+range*(low_steps-1+i/(double)low_steps),
                   brack_low+range*(low_steps-1+(i+1)/(double)low_steps),
                   (komplex *)parms, (komplex *)vector, -ROOTACC);  // Re-use the previous fh value  
      if (root != UNUSED)
      {
        list[(*index)++]=root;
        if (verbose >=3)
            printf("--- findroots returned a root on the energy loss side: vf = %g \n",root);
      }
     }

 // Second, find roots in energy gain side, there is always some
     range = (brack_high-brack_mid)/(double)high_steps;
     for (i=0; i<high_steps; i++) // Close to zero: small steps
     {
         if (i==0)
            root = zridd(f, brack_mid+range*i/(double)high_steps, 
                   brack_mid+range*(i+1)/(double)high_steps, 
                   (komplex *)parms, (komplex *)vector, ROOTACC);
        else
            root = zridd(f, brack_mid+range*i/(double)high_steps, 
                   brack_mid+range*(i+1)/(double)high_steps, 
                   (komplex *)parms, (komplex *)vector, -ROOTACC);  // Re-use the previous fh value
      if (root != UNUSED)
      {
        list[(*index)++]=root;
        if (verbose >= 3)
            printf("*** findroots returned a root on the energy gain side: vf = %g \n",root);
      }
     }
     
     range = brack_high-brack_mid;
     for (i=1; i<high_steps; i++)  // Larger steps away from zero
     {
      root = zridd(f, brack_mid+range*i/(double)high_steps, 
                   brack_mid+range*(i+1)/(double)high_steps, 
                   (komplex *)parms, (komplex *)vector, -ROOTACC);  // Re-use the previous fh value - there was always one
      if (root != UNUSED)
      {
        list[(*index)++]=root;
        if (verbose >= 3)
            printf("*** findroots returned a root on the energy gain side: vf = %g \n",root);
      }
     }
     
    }
  
#pragma acc routine 
  int findroots_gpu(double brack_low, double brack_mid, double brack_high, double *list, int* index, double *parms, komplex* vector, int low_steps, int high_steps)
    {
      double root,range=brack_mid-brack_low;
      int i;

     for (i=0; i<low_steps; i++)
     {
       root = zridd_gpu(brack_low+range*i/(int)low_steps,
                   brack_low+range*(i+1)/(int)low_steps,
                   (komplex *)parms, (komplex *)vector, ROOTACC);
      if (root != UNUSED)
      {
        list[(*index)++]=root;
      }
     } 
    for (i=0; i<high_steps; i++)
     {
      root = zridd_gpu(brack_mid+range*i/(int)high_steps,
                   brack_high+range*(i+1)/(int)high_steps,
                   (komplex *)parms, (komplex *)vector, ROOTACC);
       if (root != UNUSED)
      {
        list[(*index)++]=root;
      }
     }
    }
  
#undef UNUSED
#undef MAXRIDD
#endif

/* Shared user declarations for all components types 'Slit'. */
void slit_print_if(int condition, char* level, char* message, char* component){
  if (condition) fprintf(stderr, "Slit: %s: %s: %s\n", component, level, message);
} 
void slit_error_if(int condition, char* message, char* component){
  slit_print_if(condition, "Error", message, component);
  if (condition) exit(-1);
}
void slit_warning_if(int condition, char* message, char* component){
  slit_print_if(condition, "Warning", message, component);
}



/* ************************************************************************** */
/*             End of SHARE user declarations for all components              */
/* ************************************************************************** */


/* ********************** component definition declarations. **************** */

/* component origin=Progress_bar() [1] DECLARE */
/* Parameter definition for component type 'Progress_bar' */
struct _struct_Progress_bar_parameters {
  /* Component type 'Progress_bar' setting parameters */
  char profile[16384];
  MCNUM percent;
  MCNUM flag_save;
  MCNUM minutes;
  /* Component type 'Progress_bar' private parameters */
  double  IntermediateCnts;
  time_t  StartTime;
  time_t  EndTime;
  time_t  CurrentTime;
}; /* _struct_Progress_bar_parameters */
typedef struct _struct_Progress_bar_parameters _class_Progress_bar_parameters;

/* Parameters for component type 'Progress_bar' */
struct _struct_Progress_bar {
  char     _name[256]; /* e.g. origin */
  char     _type[256]; /* Progress_bar */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Progress_bar_parameters _parameters;
};
typedef struct _struct_Progress_bar _class_Progress_bar;
_class_Progress_bar _origin_var;
#pragma acc declare create ( _origin_var )

/* component source=Source_Maxwell_3() [2] DECLARE */
/* Parameter definition for component type 'Source_Maxwell_3' */
struct _struct_Source_Maxwell_3_parameters {
  /* Component type 'Source_Maxwell_3' setting parameters */
  MCNUM size;
  MCNUM yheight;
  MCNUM xwidth;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM dist;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM T1;
  MCNUM T2;
  MCNUM T3;
  MCNUM I1;
  MCNUM I2;
  MCNUM I3;
  long target_index;
  MCNUM lambda0;
  MCNUM dlambda;
  /* Component type 'Source_Maxwell_3' private parameters */
  double  l_range;
  double  w_mult;
  double  w_source;
  double  h_source;
}; /* _struct_Source_Maxwell_3_parameters */
typedef struct _struct_Source_Maxwell_3_parameters _class_Source_Maxwell_3_parameters;

/* Parameters for component type 'Source_Maxwell_3' */
struct _struct_Source_Maxwell_3 {
  char     _name[256]; /* e.g. source */
  char     _type[256]; /* Source_Maxwell_3 */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Source_Maxwell_3_parameters _parameters;
};
typedef struct _struct_Source_Maxwell_3 _class_Source_Maxwell_3;
_class_Source_Maxwell_3 _source_var;
#pragma acc declare create ( _source_var )

/* component l_monitor_before_mono=L_monitor() [3] DECLARE */
/* Parameter definition for component type 'L_monitor' */
struct _struct_L_monitor_parameters {
  /* Component type 'L_monitor' setting parameters */
  MCNUM nL;
  char filename[16384];
  long nowritefile;
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM Lmin;
  MCNUM Lmax;
  MCNUM restore_neutron;
  /* Component type 'L_monitor' private parameters */
  DArray1d  L_N;
  DArray1d  L_p;
  DArray1d  L_p2;
}; /* _struct_L_monitor_parameters */
typedef struct _struct_L_monitor_parameters _class_L_monitor_parameters;

/* Parameters for component type 'L_monitor' */
struct _struct_L_monitor {
  char     _name[256]; /* e.g. l_monitor_before_mono */
  char     _type[256]; /* L_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_L_monitor_parameters _parameters;
};
typedef struct _struct_L_monitor _class_L_monitor;
_class_L_monitor _l_monitor_before_mono_var;
#pragma acc declare create ( _l_monitor_before_mono_var )

/* component monochromator_flat=Monochromator_flat() [4] DECLARE */
/* Parameter definition for component type 'Monochromator_flat' */
struct _struct_Monochromator_flat_parameters {
  /* Component type 'Monochromator_flat' setting parameters */
  MCNUM zmin;
  MCNUM zmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM zwidth;
  MCNUM yheight;
  MCNUM mosaich;
  MCNUM mosaicv;
  MCNUM r0;
  MCNUM Q;
  MCNUM DM;
  /* Component type 'Monochromator_flat' private parameters */
  double  mos_rms_y;
  double  mos_rms_z;
  double  mos_rms_max;
  double  mono_Q;
}; /* _struct_Monochromator_flat_parameters */
typedef struct _struct_Monochromator_flat_parameters _class_Monochromator_flat_parameters;

/* Parameters for component type 'Monochromator_flat' */
struct _struct_Monochromator_flat {
  char     _name[256]; /* e.g. monochromator_flat */
  char     _type[256]; /* Monochromator_flat */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Monochromator_flat_parameters _parameters;
};
typedef struct _struct_Monochromator_flat _class_Monochromator_flat;
_class_Monochromator_flat _monochromator_flat_var;
#pragma acc declare create ( _monochromator_flat_var )

/* component arm1=Arm() [5] DECLARE */
/* Parameter definition for component type 'Arm' */
struct _struct_Arm_parameters {
  char Arm_has_no_parameters;
}; /* _struct_Arm_parameters */
typedef struct _struct_Arm_parameters _class_Arm_parameters;

/* Parameters for component type 'Arm' */
struct _struct_Arm {
  char     _name[256]; /* e.g. arm1 */
  char     _type[256]; /* Arm */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Arm_parameters _parameters;
};
typedef struct _struct_Arm _class_Arm;
_class_Arm _arm1_var;
#pragma acc declare create ( _arm1_var )

/* component collimator_linear1=Collimator_linear() [6] DECLARE */
/* Parameter definition for component type 'Collimator_linear' */
struct _struct_Collimator_linear_parameters {
  /* Component type 'Collimator_linear' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM length;
  MCNUM divergence;
  MCNUM transmission;
  MCNUM divergenceV;
  /* Component type 'Collimator_linear' private parameters */
  double  slope;
  double  slopeV;
}; /* _struct_Collimator_linear_parameters */
typedef struct _struct_Collimator_linear_parameters _class_Collimator_linear_parameters;

/* Parameters for component type 'Collimator_linear' */
struct _struct_Collimator_linear {
  char     _name[256]; /* e.g. collimator_linear1 */
  char     _type[256]; /* Collimator_linear */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Collimator_linear_parameters _parameters;
};
typedef struct _struct_Collimator_linear _class_Collimator_linear;
_class_Collimator_linear _collimator_linear1_var;
#pragma acc declare create ( _collimator_linear1_var )

_class_L_monitor _l_monitor_before_sample_var;
#pragma acc declare create ( _l_monitor_before_sample_var )

/* component e_monitor_before_sample=E_monitor() [8] DECLARE */
/* Parameter definition for component type 'E_monitor' */
struct _struct_E_monitor_parameters {
  /* Component type 'E_monitor' setting parameters */
  MCNUM nE;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  long nowritefile;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM Emin;
  MCNUM Emax;
  MCNUM restore_neutron;
  /* Component type 'E_monitor' private parameters */
  DArray1d  E_N;
  DArray1d  E_p;
  DArray1d  E_p2;
  double  S_p;
  double  S_pE;
  double  S_pE2;
}; /* _struct_E_monitor_parameters */
typedef struct _struct_E_monitor_parameters _class_E_monitor_parameters;

/* Parameters for component type 'E_monitor' */
struct _struct_E_monitor {
  char     _name[256]; /* e.g. e_monitor_before_sample */
  char     _type[256]; /* E_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_E_monitor_parameters _parameters;
};
typedef struct _struct_E_monitor _class_E_monitor;
_class_E_monitor _e_monitor_before_sample_var;
#pragma acc declare create ( _e_monitor_before_sample_var )

/* component PSD_monitor_before_sample=PSD_monitor() [9] DECLARE */
/* Parameter definition for component type 'PSD_monitor' */
struct _struct_PSD_monitor_parameters {
  /* Component type 'PSD_monitor' setting parameters */
  MCNUM nx;
  MCNUM ny;
  char filename[16384];
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM xwidth;
  MCNUM yheight;
  MCNUM restore_neutron;
  long nowritefile;
  /* Component type 'PSD_monitor' private parameters */
  DArray2d  PSD_N;
  DArray2d  PSD_p;
  DArray2d  PSD_p2;
}; /* _struct_PSD_monitor_parameters */
typedef struct _struct_PSD_monitor_parameters _class_PSD_monitor_parameters;

/* Parameters for component type 'PSD_monitor' */
struct _struct_PSD_monitor {
  char     _name[256]; /* e.g. PSD_monitor_before_sample */
  char     _type[256]; /* PSD_monitor */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_PSD_monitor_parameters _parameters;
};
typedef struct _struct_PSD_monitor _class_PSD_monitor;
_class_PSD_monitor _PSD_monitor_before_sample_var;
#pragma acc declare create ( _PSD_monitor_before_sample_var )

_class_Arm _arm2_var;
#pragma acc declare create ( _arm2_var )

/* component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() [11] DECLARE */
/* Parameter definition for component type 'Phonon_BvK_PG_eigenvector' */
struct _struct_Phonon_BvK_PG_eigenvector_parameters {
  /* Component type 'Phonon_BvK_PG_eigenvector' setting parameters */
  MCNUM hh;
  MCNUM kk;
  MCNUM ll;
  MCNUM radius;
  MCNUM yheight;
  MCNUM sigma_abs;
  MCNUM sigma_inc;
  MCNUM DW;
  MCNUM T;
  MCNUM focus_r;
  MCNUM focus_xw;
  MCNUM focus_yh;
  MCNUM focus_aw;
  MCNUM focus_ah;
  MCNUM target_x;
  MCNUM target_y;
  MCNUM target_z;
  long target_index;
  long mode_input;
  long e_steps_low;
  long e_steps_high;
  long verbose_input;
  long dispersion;
  /* Component type 'Phonon_BvK_PG_eigenvector' private parameters */
  double  V_rho;
  double  V_my_s;
  double  V_my_a_v;
  komplex ** Matrix;
  double  q[3];
  double  qx;
  double  qy;
  double  qz;
  double  q_x;
  double  q_y;
  double  q_z;
  komplex  eigenvectormode[DIM];
  komplex  p_call[15];
}; /* _struct_Phonon_BvK_PG_eigenvector_parameters */
typedef struct _struct_Phonon_BvK_PG_eigenvector_parameters _class_Phonon_BvK_PG_eigenvector_parameters;

/* Parameters for component type 'Phonon_BvK_PG_eigenvector' */
struct _struct_Phonon_BvK_PG_eigenvector {
  char     _name[256]; /* e.g. phonon_bvk_pg */
  char     _type[256]; /* Phonon_BvK_PG_eigenvector */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Phonon_BvK_PG_eigenvector_parameters _parameters;
};
typedef struct _struct_Phonon_BvK_PG_eigenvector _class_Phonon_BvK_PG_eigenvector;
_class_Phonon_BvK_PG_eigenvector _phonon_bvk_pg_var;
#pragma acc declare create ( _phonon_bvk_pg_var )

_class_Arm _arm3_var;
#pragma acc declare create ( _arm3_var )

_class_Collimator_linear _collimator_linear2_var;
#pragma acc declare create ( _collimator_linear2_var )

_class_PSD_monitor _psd_monitoraftersample_var;
#pragma acc declare create ( _psd_monitoraftersample_var )

/* component slit1=Slit() [15] DECLARE */
/* Parameter definition for component type 'Slit' */
struct _struct_Slit_parameters {
  /* Component type 'Slit' setting parameters */
  MCNUM xmin;
  MCNUM xmax;
  MCNUM ymin;
  MCNUM ymax;
  MCNUM radius;
  MCNUM xwidth;
  MCNUM yheight;
  /* Component type 'Slit' private parameters */
  char  isradial;
}; /* _struct_Slit_parameters */
typedef struct _struct_Slit_parameters _class_Slit_parameters;

/* Parameters for component type 'Slit' */
struct _struct_Slit {
  char     _name[256]; /* e.g. slit1 */
  char     _type[256]; /* Slit */
  long     _index; /* e.g. 2 index in TRACE list */
  Coords   _position_absolute;
  Coords   _position_relative; /* wrt PREVIOUS */
  Rotation _rotation_absolute;
  Rotation _rotation_relative; /* wrt PREVIOUS */
  int      _rotation_is_identity;
  int      _position_relative_is_zero;
  _class_Slit_parameters _parameters;
};
typedef struct _struct_Slit _class_Slit;
_class_Slit _slit1_var;
#pragma acc declare create ( _slit1_var )

_class_E_monitor _e_monitor_beforeana_var;
#pragma acc declare create ( _e_monitor_beforeana_var )

_class_Monochromator_flat _analyzer_var;
#pragma acc declare create ( _analyzer_var )

_class_Arm _arm4_var;
#pragma acc declare create ( _arm4_var )

_class_Slit _slit2_var;
#pragma acc declare create ( _slit2_var )

_class_E_monitor _Emon_after_analyzer_var;
#pragma acc declare create ( _Emon_after_analyzer_var )

_class_PSD_monitor _psd_detector_var;
#pragma acc declare create ( _psd_detector_var )

int mcNUMCOMP = 21;

/* User declarations from instrument definition. Can define functions. */
  double Gqx,Gqy,Gqz;
//double Ef=24.8
double Ei;
//meV
double qx,qy,qz;
double thetaM;
double twothetaS;
double thetaA;
double A3;
double QM;
double alpha;
double lambda_i;
double SMALL__NUMBER;

#undef compcurname
#undef compcurtype
#undef compcurindex
/* end of instrument 'template_body_simple' and components DECLARE */

/* *****************************************************************************
* instrument 'template_body_simple' and components INITIALISE
***************************************************************************** */

/* component origin=Progress_bar() SETTING, POSITION/ROTATION */
int _origin_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_origin_setpos] component origin=Progress_bar() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:57]");
  stracpy(_origin_var._name, "origin", 16384);
  stracpy(_origin_var._type, "Progress_bar", 16384);
  _origin_var._index=1;
  if("NULL" && strlen("NULL"))
    stracpy(_origin_var._parameters.profile, "NULL" ? "NULL" : "", 16384);
  else 
  _origin_var._parameters.profile[0]='\0';
  _origin_var._parameters.percent = 10;
  _origin_var._parameters.flag_save = 0;
  _origin_var._parameters.minutes = 0;


  /* component origin=Progress_bar() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(_origin_var._rotation_absolute,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_copy(_origin_var._rotation_relative, _origin_var._rotation_absolute);
    _origin_var._rotation_is_identity =  rot_test_identity(_origin_var._rotation_relative);
    _origin_var._position_absolute = coords_set(
      0, 0, 0);
    tc1 = coords_neg(_origin_var._position_absolute);
    _origin_var._position_relative = rot_apply(_origin_var._rotation_absolute, tc1);
  } /* origin=Progress_bar() AT ROTATED */
  DEBUG_COMPONENT("origin", _origin_var._position_absolute, _origin_var._rotation_absolute);
  instrument->_position_absolute[1] = _origin_var._position_absolute;
  instrument->_position_relative[1] = _origin_var._position_relative;
    _origin_var._position_relative_is_zero =  coords_test_zero(_origin_var._position_relative);
  instrument->counter_N[1]  = instrument->counter_P[1] = instrument->counter_P2[1] = 0;
  instrument->counter_AbsorbProp[1]= 0;
  return(0);
} /* _origin_setpos */

/* component source=Source_Maxwell_3() SETTING, POSITION/ROTATION */
int _source_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_source_setpos] component source=Source_Maxwell_3() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//sources/Source_Maxwell_3.comp:88]");
  stracpy(_source_var._name, "source", 16384);
  stracpy(_source_var._type, "Source_Maxwell_3", 16384);
  _source_var._index=2;
  _source_var._parameters.size = 0;
  _source_var._parameters.yheight = _instrument_var._parameters.width;
  _source_var._parameters.xwidth = _instrument_var._parameters.width;
  _source_var._parameters.Lmin = lambda_i - _instrument_var._parameters.Dlambda / 2;
  _source_var._parameters.Lmax = lambda_i + _instrument_var._parameters.Dlambda / 2;
  _source_var._parameters.dist = 7.5;
  _source_var._parameters.focus_xw = _instrument_var._parameters.width;
  _source_var._parameters.focus_yh = _instrument_var._parameters.width;
  _source_var._parameters.T1 = 300;
  _source_var._parameters.T2 = 300;
  _source_var._parameters.T3 = 300;
  _source_var._parameters.I1 = 1E15;
  _source_var._parameters.I2 = 1E15;
  _source_var._parameters.I3 = 1E15;
  _source_var._parameters.target_index = + 1;
  _source_var._parameters.lambda0 = 0;
  _source_var._parameters.dlambda = 0;


  /* component source=Source_Maxwell_3() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _origin_var._rotation_absolute, _source_var._rotation_absolute);
    rot_transpose(_origin_var._rotation_absolute, tr1);
    rot_mul(_source_var._rotation_absolute, tr1, _source_var._rotation_relative);
    _source_var._rotation_is_identity =  rot_test_identity(_source_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_origin_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _source_var._position_absolute = coords_add(_origin_var._position_absolute, tc2);
    tc1 = coords_sub(_origin_var._position_absolute, _source_var._position_absolute);
    _source_var._position_relative = rot_apply(_source_var._rotation_absolute, tc1);
  } /* source=Source_Maxwell_3() AT ROTATED */
  DEBUG_COMPONENT("source", _source_var._position_absolute, _source_var._rotation_absolute);
  instrument->_position_absolute[2] = _source_var._position_absolute;
  instrument->_position_relative[2] = _source_var._position_relative;
    _source_var._position_relative_is_zero =  coords_test_zero(_source_var._position_relative);
  instrument->counter_N[2]  = instrument->counter_P[2] = instrument->counter_P2[2] = 0;
  instrument->counter_AbsorbProp[2]= 0;
  return(0);
} /* _source_setpos */

/* component l_monitor_before_mono=L_monitor() SETTING, POSITION/ROTATION */
int _l_monitor_before_mono_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_l_monitor_before_mono_setpos] component l_monitor_before_mono=L_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:66]");
  stracpy(_l_monitor_before_mono_var._name, "l_monitor_before_mono", 16384);
  stracpy(_l_monitor_before_mono_var._type, "L_monitor", 16384);
  _l_monitor_before_mono_var._index=3;
  _l_monitor_before_mono_var._parameters.nL = 200;
  if("lbeforemono" && strlen("lbeforemono"))
    stracpy(_l_monitor_before_mono_var._parameters.filename, "lbeforemono" ? "lbeforemono" : "", 16384);
  else 
  _l_monitor_before_mono_var._parameters.filename[0]='\0';
  _l_monitor_before_mono_var._parameters.nowritefile = 0;
  _l_monitor_before_mono_var._parameters.xmin = -0.05;
  _l_monitor_before_mono_var._parameters.xmax = 0.05;
  _l_monitor_before_mono_var._parameters.ymin = -0.05;
  _l_monitor_before_mono_var._parameters.ymax = 0.05;
  _l_monitor_before_mono_var._parameters.xwidth = 0.16;
  _l_monitor_before_mono_var._parameters.yheight = 0.25;
  _l_monitor_before_mono_var._parameters.Lmin = 0;
  _l_monitor_before_mono_var._parameters.Lmax = 4;
  _l_monitor_before_mono_var._parameters.restore_neutron = 1;


  /* component l_monitor_before_mono=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _l_monitor_before_mono_var._rotation_absolute);
    rot_transpose(_source_var._rotation_absolute, tr1);
    rot_mul(_l_monitor_before_mono_var._rotation_absolute, tr1, _l_monitor_before_mono_var._rotation_relative);
    _l_monitor_before_mono_var._rotation_is_identity =  rot_test_identity(_l_monitor_before_mono_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 7);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _l_monitor_before_mono_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_source_var._position_absolute, _l_monitor_before_mono_var._position_absolute);
    _l_monitor_before_mono_var._position_relative = rot_apply(_l_monitor_before_mono_var._rotation_absolute, tc1);
  } /* l_monitor_before_mono=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("l_monitor_before_mono", _l_monitor_before_mono_var._position_absolute, _l_monitor_before_mono_var._rotation_absolute);
  instrument->_position_absolute[3] = _l_monitor_before_mono_var._position_absolute;
  instrument->_position_relative[3] = _l_monitor_before_mono_var._position_relative;
    _l_monitor_before_mono_var._position_relative_is_zero =  coords_test_zero(_l_monitor_before_mono_var._position_relative);
  instrument->counter_N[3]  = instrument->counter_P[3] = instrument->counter_P2[3] = 0;
  instrument->counter_AbsorbProp[3]= 0;
  return(0);
} /* _l_monitor_before_mono_setpos */

/* component monochromator_flat=Monochromator_flat() SETTING, POSITION/ROTATION */
int _monochromator_flat_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_monochromator_flat_setpos] component monochromator_flat=Monochromator_flat() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Monochromator_flat.comp:103]");
  stracpy(_monochromator_flat_var._name, "monochromator_flat", 16384);
  stracpy(_monochromator_flat_var._type, "Monochromator_flat", 16384);
  _monochromator_flat_var._index=4;
  _monochromator_flat_var._parameters.zmin = -0.05;
  _monochromator_flat_var._parameters.zmax = 0.05;
  _monochromator_flat_var._parameters.ymin = -0.05;
  _monochromator_flat_var._parameters.ymax = 0.05;
  _monochromator_flat_var._parameters.zwidth = _instrument_var._parameters.width;
  _monochromator_flat_var._parameters.yheight = _instrument_var._parameters.width;
  _monochromator_flat_var._parameters.mosaich = 30;
  _monochromator_flat_var._parameters.mosaicv = 30;
  _monochromator_flat_var._parameters.r0 = 0.7;
  _monochromator_flat_var._parameters.Q = QM;
  _monochromator_flat_var._parameters.DM = 0;


  /* component monochromator_flat=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (thetaM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _source_var._rotation_absolute, _monochromator_flat_var._rotation_absolute);
    rot_transpose(_l_monitor_before_mono_var._rotation_absolute, tr1);
    rot_mul(_monochromator_flat_var._rotation_absolute, tr1, _monochromator_flat_var._rotation_relative);
    _monochromator_flat_var._rotation_is_identity =  rot_test_identity(_monochromator_flat_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 7.5);
    rot_transpose(_source_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _monochromator_flat_var._position_absolute = coords_add(_source_var._position_absolute, tc2);
    tc1 = coords_sub(_l_monitor_before_mono_var._position_absolute, _monochromator_flat_var._position_absolute);
    _monochromator_flat_var._position_relative = rot_apply(_monochromator_flat_var._rotation_absolute, tc1);
  } /* monochromator_flat=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("monochromator_flat", _monochromator_flat_var._position_absolute, _monochromator_flat_var._rotation_absolute);
  instrument->_position_absolute[4] = _monochromator_flat_var._position_absolute;
  instrument->_position_relative[4] = _monochromator_flat_var._position_relative;
    _monochromator_flat_var._position_relative_is_zero =  coords_test_zero(_monochromator_flat_var._position_relative);
  instrument->counter_N[4]  = instrument->counter_P[4] = instrument->counter_P2[4] = 0;
  instrument->counter_AbsorbProp[4]= 0;
  return(0);
} /* _monochromator_flat_setpos */

/* component arm1=Arm() SETTING, POSITION/ROTATION */
int _arm1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_arm1_setpos] component arm1=Arm() SETTING [Arm:0]");
  stracpy(_arm1_var._name, "arm1", 16384);
  stracpy(_arm1_var._type, "Arm", 16384);
  _arm1_var._index=5;
  /* component arm1=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (thetaM)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _monochromator_flat_var._rotation_absolute, _arm1_var._rotation_absolute);
    rot_transpose(_monochromator_flat_var._rotation_absolute, tr1);
    rot_mul(_arm1_var._rotation_absolute, tr1, _arm1_var._rotation_relative);
    _arm1_var._rotation_is_identity =  rot_test_identity(_arm1_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_monochromator_flat_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _arm1_var._position_absolute = coords_add(_monochromator_flat_var._position_absolute, tc2);
    tc1 = coords_sub(_monochromator_flat_var._position_absolute, _arm1_var._position_absolute);
    _arm1_var._position_relative = rot_apply(_arm1_var._rotation_absolute, tc1);
  } /* arm1=Arm() AT ROTATED */
  DEBUG_COMPONENT("arm1", _arm1_var._position_absolute, _arm1_var._rotation_absolute);
  instrument->_position_absolute[5] = _arm1_var._position_absolute;
  instrument->_position_relative[5] = _arm1_var._position_relative;
    _arm1_var._position_relative_is_zero =  coords_test_zero(_arm1_var._position_relative);
  instrument->counter_N[5]  = instrument->counter_P[5] = instrument->counter_P2[5] = 0;
  instrument->counter_AbsorbProp[5]= 0;
  return(0);
} /* _arm1_setpos */

/* component collimator_linear1=Collimator_linear() SETTING, POSITION/ROTATION */
int _collimator_linear1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_collimator_linear1_setpos] component collimator_linear1=Collimator_linear() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Collimator_linear.comp:55]");
  stracpy(_collimator_linear1_var._name, "collimator_linear1", 16384);
  stracpy(_collimator_linear1_var._type, "Collimator_linear", 16384);
  _collimator_linear1_var._index=6;
  _collimator_linear1_var._parameters.xmin = -0.02;
  _collimator_linear1_var._parameters.xmax = 0.02;
  _collimator_linear1_var._parameters.ymin = -0.05;
  _collimator_linear1_var._parameters.ymax = 0.05;
  _collimator_linear1_var._parameters.xwidth = 0.1;
  _collimator_linear1_var._parameters.yheight = 0.2;
  _collimator_linear1_var._parameters.length = 0.2;
  _collimator_linear1_var._parameters.divergence = _instrument_var._parameters.coll;
  _collimator_linear1_var._parameters.transmission = 1;
  _collimator_linear1_var._parameters.divergenceV = 0;


  /* component collimator_linear1=Collimator_linear() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _arm1_var._rotation_absolute, _collimator_linear1_var._rotation_absolute);
    rot_transpose(_monochromator_flat_var._rotation_absolute, tr1);
    rot_mul(_collimator_linear1_var._rotation_absolute, tr1, _collimator_linear1_var._rotation_relative);
    _collimator_linear1_var._rotation_is_identity =  rot_test_identity(_collimator_linear1_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.7);
    rot_transpose(_arm1_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _collimator_linear1_var._position_absolute = coords_add(_arm1_var._position_absolute, tc2);
    tc1 = coords_sub(_monochromator_flat_var._position_absolute, _collimator_linear1_var._position_absolute);
    _collimator_linear1_var._position_relative = rot_apply(_collimator_linear1_var._rotation_absolute, tc1);
  } /* collimator_linear1=Collimator_linear() AT ROTATED */
  DEBUG_COMPONENT("collimator_linear1", _collimator_linear1_var._position_absolute, _collimator_linear1_var._rotation_absolute);
  instrument->_position_absolute[6] = _collimator_linear1_var._position_absolute;
  instrument->_position_relative[6] = _collimator_linear1_var._position_relative;
    _collimator_linear1_var._position_relative_is_zero =  coords_test_zero(_collimator_linear1_var._position_relative);
  instrument->counter_N[6]  = instrument->counter_P[6] = instrument->counter_P2[6] = 0;
  instrument->counter_AbsorbProp[6]= 0;
  return(0);
} /* _collimator_linear1_setpos */

/* component l_monitor_before_sample=L_monitor() SETTING, POSITION/ROTATION */
int _l_monitor_before_sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_l_monitor_before_sample_setpos] component l_monitor_before_sample=L_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:66]");
  stracpy(_l_monitor_before_sample_var._name, "l_monitor_before_sample", 16384);
  stracpy(_l_monitor_before_sample_var._type, "L_monitor", 16384);
  _l_monitor_before_sample_var._index=7;
  _l_monitor_before_sample_var._parameters.nL = 300;
  if("L_beforesampla" && strlen("L_beforesampla"))
    stracpy(_l_monitor_before_sample_var._parameters.filename, "L_beforesampla" ? "L_beforesampla" : "", 16384);
  else 
  _l_monitor_before_sample_var._parameters.filename[0]='\0';
  _l_monitor_before_sample_var._parameters.nowritefile = 0;
  _l_monitor_before_sample_var._parameters.xmin = -0.05;
  _l_monitor_before_sample_var._parameters.xmax = 0.05;
  _l_monitor_before_sample_var._parameters.ymin = -0.05;
  _l_monitor_before_sample_var._parameters.ymax = 0.05;
  _l_monitor_before_sample_var._parameters.xwidth = 0;
  _l_monitor_before_sample_var._parameters.yheight = 0;
  _l_monitor_before_sample_var._parameters.Lmin = lambda_i - _instrument_var._parameters.Dlambda;
  _l_monitor_before_sample_var._parameters.Lmax = lambda_i + _instrument_var._parameters.Dlambda;
  _l_monitor_before_sample_var._parameters.restore_neutron = 1;


  /* component l_monitor_before_sample=L_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _arm1_var._rotation_absolute, _l_monitor_before_sample_var._rotation_absolute);
    rot_transpose(_collimator_linear1_var._rotation_absolute, tr1);
    rot_mul(_l_monitor_before_sample_var._rotation_absolute, tr1, _l_monitor_before_sample_var._rotation_relative);
    _l_monitor_before_sample_var._rotation_is_identity =  rot_test_identity(_l_monitor_before_sample_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.9);
    rot_transpose(_arm1_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _l_monitor_before_sample_var._position_absolute = coords_add(_arm1_var._position_absolute, tc2);
    tc1 = coords_sub(_collimator_linear1_var._position_absolute, _l_monitor_before_sample_var._position_absolute);
    _l_monitor_before_sample_var._position_relative = rot_apply(_l_monitor_before_sample_var._rotation_absolute, tc1);
  } /* l_monitor_before_sample=L_monitor() AT ROTATED */
  DEBUG_COMPONENT("l_monitor_before_sample", _l_monitor_before_sample_var._position_absolute, _l_monitor_before_sample_var._rotation_absolute);
  instrument->_position_absolute[7] = _l_monitor_before_sample_var._position_absolute;
  instrument->_position_relative[7] = _l_monitor_before_sample_var._position_relative;
    _l_monitor_before_sample_var._position_relative_is_zero =  coords_test_zero(_l_monitor_before_sample_var._position_relative);
  instrument->counter_N[7]  = instrument->counter_P[7] = instrument->counter_P2[7] = 0;
  instrument->counter_AbsorbProp[7]= 0;
  return(0);
} /* _l_monitor_before_sample_setpos */

/* component e_monitor_before_sample=E_monitor() SETTING, POSITION/ROTATION */
int _e_monitor_before_sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_e_monitor_before_sample_setpos] component e_monitor_before_sample=E_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:69]");
  stracpy(_e_monitor_before_sample_var._name, "e_monitor_before_sample", 16384);
  stracpy(_e_monitor_before_sample_var._type, "E_monitor", 16384);
  _e_monitor_before_sample_var._index=8;
  _e_monitor_before_sample_var._parameters.nE = 100;
  if("E_before_sample" && strlen("E_before_sample"))
    stracpy(_e_monitor_before_sample_var._parameters.filename, "E_before_sample" ? "E_before_sample" : "", 16384);
  else 
  _e_monitor_before_sample_var._parameters.filename[0]='\0';
  _e_monitor_before_sample_var._parameters.xmin = -0.05;
  _e_monitor_before_sample_var._parameters.xmax = 0.05;
  _e_monitor_before_sample_var._parameters.ymin = -0.05;
  _e_monitor_before_sample_var._parameters.ymax = 0.05;
  _e_monitor_before_sample_var._parameters.nowritefile = 0;
  _e_monitor_before_sample_var._parameters.xwidth = 0;
  _e_monitor_before_sample_var._parameters.yheight = 0;
  _e_monitor_before_sample_var._parameters.Emin = Ei -2;
  _e_monitor_before_sample_var._parameters.Emax = Ei + 2;
  _e_monitor_before_sample_var._parameters.restore_neutron = 1;


  /* component e_monitor_before_sample=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _l_monitor_before_sample_var._rotation_absolute, _e_monitor_before_sample_var._rotation_absolute);
    rot_transpose(_l_monitor_before_sample_var._rotation_absolute, tr1);
    rot_mul(_e_monitor_before_sample_var._rotation_absolute, tr1, _e_monitor_before_sample_var._rotation_relative);
    _e_monitor_before_sample_var._rotation_is_identity =  rot_test_identity(_e_monitor_before_sample_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.01);
    rot_transpose(_l_monitor_before_sample_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _e_monitor_before_sample_var._position_absolute = coords_add(_l_monitor_before_sample_var._position_absolute, tc2);
    tc1 = coords_sub(_l_monitor_before_sample_var._position_absolute, _e_monitor_before_sample_var._position_absolute);
    _e_monitor_before_sample_var._position_relative = rot_apply(_e_monitor_before_sample_var._rotation_absolute, tc1);
  } /* e_monitor_before_sample=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("e_monitor_before_sample", _e_monitor_before_sample_var._position_absolute, _e_monitor_before_sample_var._rotation_absolute);
  instrument->_position_absolute[8] = _e_monitor_before_sample_var._position_absolute;
  instrument->_position_relative[8] = _e_monitor_before_sample_var._position_relative;
    _e_monitor_before_sample_var._position_relative_is_zero =  coords_test_zero(_e_monitor_before_sample_var._position_relative);
  instrument->counter_N[8]  = instrument->counter_P[8] = instrument->counter_P2[8] = 0;
  instrument->counter_AbsorbProp[8]= 0;
  return(0);
} /* _e_monitor_before_sample_setpos */

/* component PSD_monitor_before_sample=PSD_monitor() SETTING, POSITION/ROTATION */
int _PSD_monitor_before_sample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_PSD_monitor_before_sample_setpos] component PSD_monitor_before_sample=PSD_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:62]");
  stracpy(_PSD_monitor_before_sample_var._name, "PSD_monitor_before_sample", 16384);
  stracpy(_PSD_monitor_before_sample_var._type, "PSD_monitor", 16384);
  _PSD_monitor_before_sample_var._index=9;
  _PSD_monitor_before_sample_var._parameters.nx = 200;
  _PSD_monitor_before_sample_var._parameters.ny = 200;
  if("PSD_before_sample" && strlen("PSD_before_sample"))
    stracpy(_PSD_monitor_before_sample_var._parameters.filename, "PSD_before_sample" ? "PSD_before_sample" : "", 16384);
  else 
  _PSD_monitor_before_sample_var._parameters.filename[0]='\0';
  _PSD_monitor_before_sample_var._parameters.xmin = -0.05;
  _PSD_monitor_before_sample_var._parameters.xmax = 0.05;
  _PSD_monitor_before_sample_var._parameters.ymin = -0.05;
  _PSD_monitor_before_sample_var._parameters.ymax = 0.05;
  _PSD_monitor_before_sample_var._parameters.xwidth = 0.1;
  _PSD_monitor_before_sample_var._parameters.yheight = 0.1;
  _PSD_monitor_before_sample_var._parameters.restore_neutron = 1;
  _PSD_monitor_before_sample_var._parameters.nowritefile = 0;


  /* component PSD_monitor_before_sample=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _e_monitor_before_sample_var._rotation_absolute, _PSD_monitor_before_sample_var._rotation_absolute);
    rot_transpose(_e_monitor_before_sample_var._rotation_absolute, tr1);
    rot_mul(_PSD_monitor_before_sample_var._rotation_absolute, tr1, _PSD_monitor_before_sample_var._rotation_relative);
    _PSD_monitor_before_sample_var._rotation_is_identity =  rot_test_identity(_PSD_monitor_before_sample_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.01);
    rot_transpose(_e_monitor_before_sample_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _PSD_monitor_before_sample_var._position_absolute = coords_add(_e_monitor_before_sample_var._position_absolute, tc2);
    tc1 = coords_sub(_e_monitor_before_sample_var._position_absolute, _PSD_monitor_before_sample_var._position_absolute);
    _PSD_monitor_before_sample_var._position_relative = rot_apply(_PSD_monitor_before_sample_var._rotation_absolute, tc1);
  } /* PSD_monitor_before_sample=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("PSD_monitor_before_sample", _PSD_monitor_before_sample_var._position_absolute, _PSD_monitor_before_sample_var._rotation_absolute);
  instrument->_position_absolute[9] = _PSD_monitor_before_sample_var._position_absolute;
  instrument->_position_relative[9] = _PSD_monitor_before_sample_var._position_relative;
    _PSD_monitor_before_sample_var._position_relative_is_zero =  coords_test_zero(_PSD_monitor_before_sample_var._position_relative);
  instrument->counter_N[9]  = instrument->counter_P[9] = instrument->counter_P2[9] = 0;
  instrument->counter_AbsorbProp[9]= 0;
  return(0);
} /* _PSD_monitor_before_sample_setpos */

/* component arm2=Arm() SETTING, POSITION/ROTATION */
int _arm2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_arm2_setpos] component arm2=Arm() SETTING [Arm:0]");
  stracpy(_arm2_var._name, "arm2", 16384);
  stracpy(_arm2_var._type, "Arm", 16384);
  _arm2_var._index=10;
  /* component arm2=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _arm1_var._rotation_absolute, _arm2_var._rotation_absolute);
    rot_transpose(_PSD_monitor_before_sample_var._rotation_absolute, tr1);
    rot_mul(_arm2_var._rotation_absolute, tr1, _arm2_var._rotation_relative);
    _arm2_var._rotation_is_identity =  rot_test_identity(_arm2_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1);
    rot_transpose(_arm1_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _arm2_var._position_absolute = coords_add(_arm1_var._position_absolute, tc2);
    tc1 = coords_sub(_PSD_monitor_before_sample_var._position_absolute, _arm2_var._position_absolute);
    _arm2_var._position_relative = rot_apply(_arm2_var._rotation_absolute, tc1);
  } /* arm2=Arm() AT ROTATED */
  DEBUG_COMPONENT("arm2", _arm2_var._position_absolute, _arm2_var._rotation_absolute);
  instrument->_position_absolute[10] = _arm2_var._position_absolute;
  instrument->_position_relative[10] = _arm2_var._position_relative;
    _arm2_var._position_relative_is_zero =  coords_test_zero(_arm2_var._position_relative);
  instrument->counter_N[10]  = instrument->counter_P[10] = instrument->counter_P2[10] = 0;
  instrument->counter_AbsorbProp[10]= 0;
  return(0);
} /* _arm2_setpos */

/* component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() SETTING, POSITION/ROTATION */
int _phonon_bvk_pg_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_phonon_bvk_pg_setpos] component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() SETTING [Phonon_BvK_PG_eigenvector.comp:964]");
  stracpy(_phonon_bvk_pg_var._name, "phonon_bvk_pg", 16384);
  stracpy(_phonon_bvk_pg_var._type, "Phonon_BvK_PG_eigenvector", 16384);
  _phonon_bvk_pg_var._index=11;
  _phonon_bvk_pg_var._parameters.hh = _instrument_var._parameters.h;
  _phonon_bvk_pg_var._parameters.kk = 0;
  _phonon_bvk_pg_var._parameters.ll = _instrument_var._parameters.l;
  _phonon_bvk_pg_var._parameters.radius = _instrument_var._parameters.width / 2;
  _phonon_bvk_pg_var._parameters.yheight = 2 * _instrument_var._parameters.width;
  _phonon_bvk_pg_var._parameters.sigma_abs = 0;
  _phonon_bvk_pg_var._parameters.sigma_inc = 0;
  _phonon_bvk_pg_var._parameters.DW = 1;
  _phonon_bvk_pg_var._parameters.T = _instrument_var._parameters.Temp;
  _phonon_bvk_pg_var._parameters.focus_r = 0;
  _phonon_bvk_pg_var._parameters.focus_xw = _instrument_var._parameters.width;
  _phonon_bvk_pg_var._parameters.focus_yh = _instrument_var._parameters.width;
  _phonon_bvk_pg_var._parameters.focus_aw = 0;
  _phonon_bvk_pg_var._parameters.focus_ah = 0;
  _phonon_bvk_pg_var._parameters.target_x = 0;
  _phonon_bvk_pg_var._parameters.target_y = 0;
  _phonon_bvk_pg_var._parameters.target_z = 0;
  _phonon_bvk_pg_var._parameters.target_index = 4;
  _phonon_bvk_pg_var._parameters.mode_input = _instrument_var._parameters.phononmode;
  _phonon_bvk_pg_var._parameters.e_steps_low = _instrument_var._parameters.E_steps_low;
  _phonon_bvk_pg_var._parameters.e_steps_high = _instrument_var._parameters.E_steps_high;
  _phonon_bvk_pg_var._parameters.verbose_input = _instrument_var._parameters.Verbose;
  _phonon_bvk_pg_var._parameters.dispersion = _instrument_var._parameters.DISP;


  /* component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (- A3 + _instrument_var._parameters.dA3)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _arm2_var._rotation_absolute, _phonon_bvk_pg_var._rotation_absolute);
    rot_transpose(_PSD_monitor_before_sample_var._rotation_absolute, tr1);
    rot_mul(_phonon_bvk_pg_var._rotation_absolute, tr1, _phonon_bvk_pg_var._rotation_relative);
    _phonon_bvk_pg_var._rotation_is_identity =  rot_test_identity(_phonon_bvk_pg_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_arm2_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _phonon_bvk_pg_var._position_absolute = coords_add(_arm2_var._position_absolute, tc2);
    tc1 = coords_sub(_PSD_monitor_before_sample_var._position_absolute, _phonon_bvk_pg_var._position_absolute);
    _phonon_bvk_pg_var._position_relative = rot_apply(_phonon_bvk_pg_var._rotation_absolute, tc1);
  } /* phonon_bvk_pg=Phonon_BvK_PG_eigenvector() AT ROTATED */
  DEBUG_COMPONENT("phonon_bvk_pg", _phonon_bvk_pg_var._position_absolute, _phonon_bvk_pg_var._rotation_absolute);
  instrument->_position_absolute[11] = _phonon_bvk_pg_var._position_absolute;
  instrument->_position_relative[11] = _phonon_bvk_pg_var._position_relative;
    _phonon_bvk_pg_var._position_relative_is_zero =  coords_test_zero(_phonon_bvk_pg_var._position_relative);
  instrument->counter_N[11]  = instrument->counter_P[11] = instrument->counter_P2[11] = 0;
  instrument->counter_AbsorbProp[11]= 0;
  return(0);
} /* _phonon_bvk_pg_setpos */

/* component arm3=Arm() SETTING, POSITION/ROTATION */
int _arm3_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_arm3_setpos] component arm3=Arm() SETTING [Arm:0]");
  stracpy(_arm3_var._name, "arm3", 16384);
  stracpy(_arm3_var._type, "Arm", 16384);
  _arm3_var._index=12;
  /* component arm3=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (- twothetaS)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _arm2_var._rotation_absolute, _arm3_var._rotation_absolute);
    rot_transpose(_phonon_bvk_pg_var._rotation_absolute, tr1);
    rot_mul(_arm3_var._rotation_absolute, tr1, _arm3_var._rotation_relative);
    _arm3_var._rotation_is_identity =  rot_test_identity(_arm3_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_phonon_bvk_pg_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _arm3_var._position_absolute = coords_add(_phonon_bvk_pg_var._position_absolute, tc2);
    tc1 = coords_sub(_phonon_bvk_pg_var._position_absolute, _arm3_var._position_absolute);
    _arm3_var._position_relative = rot_apply(_arm3_var._rotation_absolute, tc1);
  } /* arm3=Arm() AT ROTATED */
  DEBUG_COMPONENT("arm3", _arm3_var._position_absolute, _arm3_var._rotation_absolute);
  instrument->_position_absolute[12] = _arm3_var._position_absolute;
  instrument->_position_relative[12] = _arm3_var._position_relative;
    _arm3_var._position_relative_is_zero =  coords_test_zero(_arm3_var._position_relative);
  instrument->counter_N[12]  = instrument->counter_P[12] = instrument->counter_P2[12] = 0;
  instrument->counter_AbsorbProp[12]= 0;
  return(0);
} /* _arm3_setpos */

/* component collimator_linear2=Collimator_linear() SETTING, POSITION/ROTATION */
int _collimator_linear2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_collimator_linear2_setpos] component collimator_linear2=Collimator_linear() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Collimator_linear.comp:55]");
  stracpy(_collimator_linear2_var._name, "collimator_linear2", 16384);
  stracpy(_collimator_linear2_var._type, "Collimator_linear", 16384);
  _collimator_linear2_var._index=13;
  _collimator_linear2_var._parameters.xmin = -0.02;
  _collimator_linear2_var._parameters.xmax = 0.02;
  _collimator_linear2_var._parameters.ymin = -0.05;
  _collimator_linear2_var._parameters.ymax = 0.05;
  _collimator_linear2_var._parameters.xwidth = 0.1;
  _collimator_linear2_var._parameters.yheight = 0.2;
  _collimator_linear2_var._parameters.length = 0.2;
  _collimator_linear2_var._parameters.divergence = _instrument_var._parameters.coll;
  _collimator_linear2_var._parameters.transmission = 1;
  _collimator_linear2_var._parameters.divergenceV = 0;


  /* component collimator_linear2=Collimator_linear() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _arm3_var._rotation_absolute, _collimator_linear2_var._rotation_absolute);
    rot_transpose(_phonon_bvk_pg_var._rotation_absolute, tr1);
    rot_mul(_collimator_linear2_var._rotation_absolute, tr1, _collimator_linear2_var._rotation_relative);
    _collimator_linear2_var._rotation_is_identity =  rot_test_identity(_collimator_linear2_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.5);
    rot_transpose(_arm3_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _collimator_linear2_var._position_absolute = coords_add(_arm3_var._position_absolute, tc2);
    tc1 = coords_sub(_phonon_bvk_pg_var._position_absolute, _collimator_linear2_var._position_absolute);
    _collimator_linear2_var._position_relative = rot_apply(_collimator_linear2_var._rotation_absolute, tc1);
  } /* collimator_linear2=Collimator_linear() AT ROTATED */
  DEBUG_COMPONENT("collimator_linear2", _collimator_linear2_var._position_absolute, _collimator_linear2_var._rotation_absolute);
  instrument->_position_absolute[13] = _collimator_linear2_var._position_absolute;
  instrument->_position_relative[13] = _collimator_linear2_var._position_relative;
    _collimator_linear2_var._position_relative_is_zero =  coords_test_zero(_collimator_linear2_var._position_relative);
  instrument->counter_N[13]  = instrument->counter_P[13] = instrument->counter_P2[13] = 0;
  instrument->counter_AbsorbProp[13]= 0;
  return(0);
} /* _collimator_linear2_setpos */

/* component psd_monitoraftersample=PSD_monitor() SETTING, POSITION/ROTATION */
int _psd_monitoraftersample_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psd_monitoraftersample_setpos] component psd_monitoraftersample=PSD_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:62]");
  stracpy(_psd_monitoraftersample_var._name, "psd_monitoraftersample", 16384);
  stracpy(_psd_monitoraftersample_var._type, "PSD_monitor", 16384);
  _psd_monitoraftersample_var._index=14;
  _psd_monitoraftersample_var._parameters.nx = 200;
  _psd_monitoraftersample_var._parameters.ny = 200;
  if("psdaftersample" && strlen("psdaftersample"))
    stracpy(_psd_monitoraftersample_var._parameters.filename, "psdaftersample" ? "psdaftersample" : "", 16384);
  else 
  _psd_monitoraftersample_var._parameters.filename[0]='\0';
  _psd_monitoraftersample_var._parameters.xmin = -0.05;
  _psd_monitoraftersample_var._parameters.xmax = 0.05;
  _psd_monitoraftersample_var._parameters.ymin = -0.05;
  _psd_monitoraftersample_var._parameters.ymax = 0.05;
  _psd_monitoraftersample_var._parameters.xwidth = 0;
  _psd_monitoraftersample_var._parameters.yheight = 0;
  _psd_monitoraftersample_var._parameters.restore_neutron = 1;
  _psd_monitoraftersample_var._parameters.nowritefile = 0;


  /* component psd_monitoraftersample=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _arm3_var._rotation_absolute, _psd_monitoraftersample_var._rotation_absolute);
    rot_transpose(_collimator_linear2_var._rotation_absolute, tr1);
    rot_mul(_psd_monitoraftersample_var._rotation_absolute, tr1, _psd_monitoraftersample_var._rotation_relative);
    _psd_monitoraftersample_var._rotation_is_identity =  rot_test_identity(_psd_monitoraftersample_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1);
    rot_transpose(_arm3_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psd_monitoraftersample_var._position_absolute = coords_add(_arm3_var._position_absolute, tc2);
    tc1 = coords_sub(_collimator_linear2_var._position_absolute, _psd_monitoraftersample_var._position_absolute);
    _psd_monitoraftersample_var._position_relative = rot_apply(_psd_monitoraftersample_var._rotation_absolute, tc1);
  } /* psd_monitoraftersample=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psd_monitoraftersample", _psd_monitoraftersample_var._position_absolute, _psd_monitoraftersample_var._rotation_absolute);
  instrument->_position_absolute[14] = _psd_monitoraftersample_var._position_absolute;
  instrument->_position_relative[14] = _psd_monitoraftersample_var._position_relative;
    _psd_monitoraftersample_var._position_relative_is_zero =  coords_test_zero(_psd_monitoraftersample_var._position_relative);
  instrument->counter_N[14]  = instrument->counter_P[14] = instrument->counter_P2[14] = 0;
  instrument->counter_AbsorbProp[14]= 0;
  return(0);
} /* _psd_monitoraftersample_setpos */

/* component slit1=Slit() SETTING, POSITION/ROTATION */
int _slit1_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit1_setpos] component slit1=Slit() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Slit.comp:66]");
  stracpy(_slit1_var._name, "slit1", 16384);
  stracpy(_slit1_var._type, "Slit", 16384);
  _slit1_var._index=15;
  _slit1_var._parameters.xmin = UNSET;
  _slit1_var._parameters.xmax = UNSET;
  _slit1_var._parameters.ymin = UNSET;
  _slit1_var._parameters.ymax = UNSET;
  _slit1_var._parameters.radius = UNSET;
  _slit1_var._parameters.xwidth = _instrument_var._parameters.width;
  _slit1_var._parameters.yheight = _instrument_var._parameters.width;


  /* component slit1=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _psd_monitoraftersample_var._rotation_absolute, _slit1_var._rotation_absolute);
    rot_transpose(_psd_monitoraftersample_var._rotation_absolute, tr1);
    rot_mul(_slit1_var._rotation_absolute, tr1, _slit1_var._rotation_relative);
    _slit1_var._rotation_is_identity =  rot_test_identity(_slit1_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.001);
    rot_transpose(_psd_monitoraftersample_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit1_var._position_absolute = coords_add(_psd_monitoraftersample_var._position_absolute, tc2);
    tc1 = coords_sub(_psd_monitoraftersample_var._position_absolute, _slit1_var._position_absolute);
    _slit1_var._position_relative = rot_apply(_slit1_var._rotation_absolute, tc1);
  } /* slit1=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit1", _slit1_var._position_absolute, _slit1_var._rotation_absolute);
  instrument->_position_absolute[15] = _slit1_var._position_absolute;
  instrument->_position_relative[15] = _slit1_var._position_relative;
    _slit1_var._position_relative_is_zero =  coords_test_zero(_slit1_var._position_relative);
  instrument->counter_N[15]  = instrument->counter_P[15] = instrument->counter_P2[15] = 0;
  instrument->counter_AbsorbProp[15]= 0;
  return(0);
} /* _slit1_setpos */

/* component e_monitor_beforeana=E_monitor() SETTING, POSITION/ROTATION */
int _e_monitor_beforeana_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_e_monitor_beforeana_setpos] component e_monitor_beforeana=E_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:69]");
  stracpy(_e_monitor_beforeana_var._name, "e_monitor_beforeana", 16384);
  stracpy(_e_monitor_beforeana_var._type, "E_monitor", 16384);
  _e_monitor_beforeana_var._index=16;
  _e_monitor_beforeana_var._parameters.nE = 250;
  if("Ebeforeana" && strlen("Ebeforeana"))
    stracpy(_e_monitor_beforeana_var._parameters.filename, "Ebeforeana" ? "Ebeforeana" : "", 16384);
  else 
  _e_monitor_beforeana_var._parameters.filename[0]='\0';
  _e_monitor_beforeana_var._parameters.xmin = -0.05;
  _e_monitor_beforeana_var._parameters.xmax = 0.05;
  _e_monitor_beforeana_var._parameters.ymin = -0.05;
  _e_monitor_beforeana_var._parameters.ymax = 0.05;
  _e_monitor_beforeana_var._parameters.nowritefile = 0;
  _e_monitor_beforeana_var._parameters.xwidth = _instrument_var._parameters.width;
  _e_monitor_beforeana_var._parameters.yheight = _instrument_var._parameters.width;
  _e_monitor_beforeana_var._parameters.Emin = 0;
  _e_monitor_beforeana_var._parameters.Emax = 25;
  _e_monitor_beforeana_var._parameters.restore_neutron = 0;


  /* component e_monitor_beforeana=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _slit1_var._rotation_absolute, _e_monitor_beforeana_var._rotation_absolute);
    rot_transpose(_slit1_var._rotation_absolute, tr1);
    rot_mul(_e_monitor_beforeana_var._rotation_absolute, tr1, _e_monitor_beforeana_var._rotation_relative);
    _e_monitor_beforeana_var._rotation_is_identity =  rot_test_identity(_e_monitor_beforeana_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.001);
    rot_transpose(_slit1_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _e_monitor_beforeana_var._position_absolute = coords_add(_slit1_var._position_absolute, tc2);
    tc1 = coords_sub(_slit1_var._position_absolute, _e_monitor_beforeana_var._position_absolute);
    _e_monitor_beforeana_var._position_relative = rot_apply(_e_monitor_beforeana_var._rotation_absolute, tc1);
  } /* e_monitor_beforeana=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("e_monitor_beforeana", _e_monitor_beforeana_var._position_absolute, _e_monitor_beforeana_var._rotation_absolute);
  instrument->_position_absolute[16] = _e_monitor_beforeana_var._position_absolute;
  instrument->_position_relative[16] = _e_monitor_beforeana_var._position_relative;
    _e_monitor_beforeana_var._position_relative_is_zero =  coords_test_zero(_e_monitor_beforeana_var._position_relative);
  instrument->counter_N[16]  = instrument->counter_P[16] = instrument->counter_P2[16] = 0;
  instrument->counter_AbsorbProp[16]= 0;
  return(0);
} /* _e_monitor_beforeana_setpos */

/* component analyzer=Monochromator_flat() SETTING, POSITION/ROTATION */
int _analyzer_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_analyzer_setpos] component analyzer=Monochromator_flat() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Monochromator_flat.comp:103]");
  stracpy(_analyzer_var._name, "analyzer", 16384);
  stracpy(_analyzer_var._type, "Monochromator_flat", 16384);
  _analyzer_var._index=17;
  _analyzer_var._parameters.zmin = -0.05;
  _analyzer_var._parameters.zmax = 0.05;
  _analyzer_var._parameters.ymin = -0.05;
  _analyzer_var._parameters.ymax = 0.05;
  _analyzer_var._parameters.zwidth = 4 * _instrument_var._parameters.width;
  _analyzer_var._parameters.yheight = 4 * _instrument_var._parameters.width;
  _analyzer_var._parameters.mosaich = 30;
  _analyzer_var._parameters.mosaicv = 30;
  _analyzer_var._parameters.r0 = 0.7;
  _analyzer_var._parameters.Q = QM;
  _analyzer_var._parameters.DM = 0;


  /* component analyzer=Monochromator_flat() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (thetaA)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _e_monitor_beforeana_var._rotation_absolute, _analyzer_var._rotation_absolute);
    rot_transpose(_e_monitor_beforeana_var._rotation_absolute, tr1);
    rot_mul(_analyzer_var._rotation_absolute, tr1, _analyzer_var._rotation_relative);
    _analyzer_var._rotation_is_identity =  rot_test_identity(_analyzer_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.1);
    rot_transpose(_e_monitor_beforeana_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _analyzer_var._position_absolute = coords_add(_e_monitor_beforeana_var._position_absolute, tc2);
    tc1 = coords_sub(_e_monitor_beforeana_var._position_absolute, _analyzer_var._position_absolute);
    _analyzer_var._position_relative = rot_apply(_analyzer_var._rotation_absolute, tc1);
  } /* analyzer=Monochromator_flat() AT ROTATED */
  DEBUG_COMPONENT("analyzer", _analyzer_var._position_absolute, _analyzer_var._rotation_absolute);
  instrument->_position_absolute[17] = _analyzer_var._position_absolute;
  instrument->_position_relative[17] = _analyzer_var._position_relative;
    _analyzer_var._position_relative_is_zero =  coords_test_zero(_analyzer_var._position_relative);
  instrument->counter_N[17]  = instrument->counter_P[17] = instrument->counter_P2[17] = 0;
  instrument->counter_AbsorbProp[17]= 0;
  return(0);
} /* _analyzer_setpos */

/* component arm4=Arm() SETTING, POSITION/ROTATION */
int _arm4_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_arm4_setpos] component arm4=Arm() SETTING [Arm:0]");
  stracpy(_arm4_var._name, "arm4", 16384);
  stracpy(_arm4_var._type, "Arm", 16384);
  _arm4_var._index=18;
  /* component arm4=Arm() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0)*DEG2RAD, (thetaA)*DEG2RAD, (0)*DEG2RAD);
    rot_mul(tr1, _analyzer_var._rotation_absolute, _arm4_var._rotation_absolute);
    rot_transpose(_analyzer_var._rotation_absolute, tr1);
    rot_mul(_arm4_var._rotation_absolute, tr1, _arm4_var._rotation_relative);
    _arm4_var._rotation_is_identity =  rot_test_identity(_arm4_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0);
    rot_transpose(_analyzer_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _arm4_var._position_absolute = coords_add(_analyzer_var._position_absolute, tc2);
    tc1 = coords_sub(_analyzer_var._position_absolute, _arm4_var._position_absolute);
    _arm4_var._position_relative = rot_apply(_arm4_var._rotation_absolute, tc1);
  } /* arm4=Arm() AT ROTATED */
  DEBUG_COMPONENT("arm4", _arm4_var._position_absolute, _arm4_var._rotation_absolute);
  instrument->_position_absolute[18] = _arm4_var._position_absolute;
  instrument->_position_relative[18] = _arm4_var._position_relative;
    _arm4_var._position_relative_is_zero =  coords_test_zero(_arm4_var._position_relative);
  instrument->counter_N[18]  = instrument->counter_P[18] = instrument->counter_P2[18] = 0;
  instrument->counter_AbsorbProp[18]= 0;
  return(0);
} /* _arm4_setpos */

/* component slit2=Slit() SETTING, POSITION/ROTATION */
int _slit2_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_slit2_setpos] component slit2=Slit() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Slit.comp:66]");
  stracpy(_slit2_var._name, "slit2", 16384);
  stracpy(_slit2_var._type, "Slit", 16384);
  _slit2_var._index=19;
  _slit2_var._parameters.xmin = UNSET;
  _slit2_var._parameters.xmax = UNSET;
  _slit2_var._parameters.ymin = UNSET;
  _slit2_var._parameters.ymax = UNSET;
  _slit2_var._parameters.radius = UNSET;
  _slit2_var._parameters.xwidth = _instrument_var._parameters.width;
  _slit2_var._parameters.yheight = _instrument_var._parameters.width;


  /* component slit2=Slit() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _arm4_var._rotation_absolute, _slit2_var._rotation_absolute);
    rot_transpose(_analyzer_var._rotation_absolute, tr1);
    rot_mul(_slit2_var._rotation_absolute, tr1, _slit2_var._rotation_relative);
    _slit2_var._rotation_is_identity =  rot_test_identity(_slit2_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 1);
    rot_transpose(_arm4_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _slit2_var._position_absolute = coords_add(_arm4_var._position_absolute, tc2);
    tc1 = coords_sub(_analyzer_var._position_absolute, _slit2_var._position_absolute);
    _slit2_var._position_relative = rot_apply(_slit2_var._rotation_absolute, tc1);
  } /* slit2=Slit() AT ROTATED */
  DEBUG_COMPONENT("slit2", _slit2_var._position_absolute, _slit2_var._rotation_absolute);
  instrument->_position_absolute[19] = _slit2_var._position_absolute;
  instrument->_position_relative[19] = _slit2_var._position_relative;
    _slit2_var._position_relative_is_zero =  coords_test_zero(_slit2_var._position_relative);
  instrument->counter_N[19]  = instrument->counter_P[19] = instrument->counter_P2[19] = 0;
  instrument->counter_AbsorbProp[19]= 0;
  return(0);
} /* _slit2_setpos */

/* component Emon_after_analyzer=E_monitor() SETTING, POSITION/ROTATION */
int _Emon_after_analyzer_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_Emon_after_analyzer_setpos] component Emon_after_analyzer=E_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:69]");
  stracpy(_Emon_after_analyzer_var._name, "Emon_after_analyzer", 16384);
  stracpy(_Emon_after_analyzer_var._type, "E_monitor", 16384);
  _Emon_after_analyzer_var._index=20;
  _Emon_after_analyzer_var._parameters.nE = 250;
  if("Eafteranalyzer" && strlen("Eafteranalyzer"))
    stracpy(_Emon_after_analyzer_var._parameters.filename, "Eafteranalyzer" ? "Eafteranalyzer" : "", 16384);
  else 
  _Emon_after_analyzer_var._parameters.filename[0]='\0';
  _Emon_after_analyzer_var._parameters.xmin = -0.05;
  _Emon_after_analyzer_var._parameters.xmax = 0.05;
  _Emon_after_analyzer_var._parameters.ymin = -0.05;
  _Emon_after_analyzer_var._parameters.ymax = 0.05;
  _Emon_after_analyzer_var._parameters.nowritefile = 0;
  _Emon_after_analyzer_var._parameters.xwidth = 0.05;
  _Emon_after_analyzer_var._parameters.yheight = 0.10;
  _Emon_after_analyzer_var._parameters.Emin = 0;
  _Emon_after_analyzer_var._parameters.Emax = 75;
  _Emon_after_analyzer_var._parameters.restore_neutron = 0;


  /* component Emon_after_analyzer=E_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _slit2_var._rotation_absolute, _Emon_after_analyzer_var._rotation_absolute);
    rot_transpose(_slit2_var._rotation_absolute, tr1);
    rot_mul(_Emon_after_analyzer_var._rotation_absolute, tr1, _Emon_after_analyzer_var._rotation_relative);
    _Emon_after_analyzer_var._rotation_is_identity =  rot_test_identity(_Emon_after_analyzer_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.001);
    rot_transpose(_slit2_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _Emon_after_analyzer_var._position_absolute = coords_add(_slit2_var._position_absolute, tc2);
    tc1 = coords_sub(_slit2_var._position_absolute, _Emon_after_analyzer_var._position_absolute);
    _Emon_after_analyzer_var._position_relative = rot_apply(_Emon_after_analyzer_var._rotation_absolute, tc1);
  } /* Emon_after_analyzer=E_monitor() AT ROTATED */
  DEBUG_COMPONENT("Emon_after_analyzer", _Emon_after_analyzer_var._position_absolute, _Emon_after_analyzer_var._rotation_absolute);
  instrument->_position_absolute[20] = _Emon_after_analyzer_var._position_absolute;
  instrument->_position_relative[20] = _Emon_after_analyzer_var._position_relative;
    _Emon_after_analyzer_var._position_relative_is_zero =  coords_test_zero(_Emon_after_analyzer_var._position_relative);
  instrument->counter_N[20]  = instrument->counter_P[20] = instrument->counter_P2[20] = 0;
  instrument->counter_AbsorbProp[20]= 0;
  return(0);
} /* _Emon_after_analyzer_setpos */

/* component psd_detector=PSD_monitor() SETTING, POSITION/ROTATION */
int _psd_detector_setpos(void)
{ /* sets initial component parameters, position and rotation */
  SIG_MESSAGE("[_psd_detector_setpos] component psd_detector=PSD_monitor() SETTING [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:62]");
  stracpy(_psd_detector_var._name, "psd_detector", 16384);
  stracpy(_psd_detector_var._type, "PSD_monitor", 16384);
  _psd_detector_var._index=21;
  _psd_detector_var._parameters.nx = 128;
  _psd_detector_var._parameters.ny = 128;
  if("PSDdetector" && strlen("PSDdetector"))
    stracpy(_psd_detector_var._parameters.filename, "PSDdetector" ? "PSDdetector" : "", 16384);
  else 
  _psd_detector_var._parameters.filename[0]='\0';
  _psd_detector_var._parameters.xmin = -0.05;
  _psd_detector_var._parameters.xmax = 0.05;
  _psd_detector_var._parameters.ymin = -0.05;
  _psd_detector_var._parameters.ymax = 0.05;
  _psd_detector_var._parameters.xwidth = 0.05;
  _psd_detector_var._parameters.yheight = 0.10;
  _psd_detector_var._parameters.restore_neutron = 0;
  _psd_detector_var._parameters.nowritefile = 0;


  /* component psd_detector=PSD_monitor() AT ROTATED */
  {
    Coords tc1, tc2;
    tc1 = coords_set(0,0,0);
    tc2 = coords_set(0,0,0);
    Rotation tr1;
    rot_set_rotation(tr1,0,0,0);
    rot_set_rotation(tr1,
      (0.0)*DEG2RAD, (0.0)*DEG2RAD, (0.0)*DEG2RAD);
    rot_mul(tr1, _Emon_after_analyzer_var._rotation_absolute, _psd_detector_var._rotation_absolute);
    rot_transpose(_Emon_after_analyzer_var._rotation_absolute, tr1);
    rot_mul(_psd_detector_var._rotation_absolute, tr1, _psd_detector_var._rotation_relative);
    _psd_detector_var._rotation_is_identity =  rot_test_identity(_psd_detector_var._rotation_relative);
    tc1 = coords_set(
      0, 0, 0.001);
    rot_transpose(_Emon_after_analyzer_var._rotation_absolute, tr1);
    tc2 = rot_apply(tr1, tc1);
    _psd_detector_var._position_absolute = coords_add(_Emon_after_analyzer_var._position_absolute, tc2);
    tc1 = coords_sub(_Emon_after_analyzer_var._position_absolute, _psd_detector_var._position_absolute);
    _psd_detector_var._position_relative = rot_apply(_psd_detector_var._rotation_absolute, tc1);
  } /* psd_detector=PSD_monitor() AT ROTATED */
  DEBUG_COMPONENT("psd_detector", _psd_detector_var._position_absolute, _psd_detector_var._rotation_absolute);
  instrument->_position_absolute[21] = _psd_detector_var._position_absolute;
  instrument->_position_relative[21] = _psd_detector_var._position_relative;
    _psd_detector_var._position_relative_is_zero =  coords_test_zero(_psd_detector_var._position_relative);
  instrument->counter_N[21]  = instrument->counter_P[21] = instrument->counter_P2[21] = 0;
  instrument->counter_AbsorbProp[21]= 0;
  return(0);
} /* _psd_detector_setpos */

_class_Progress_bar *class_Progress_bar_init(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_init] component origin=Progress_bar() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:57]");

IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", instrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
  #ifdef OPENACC
  time(&StartTime);
  #endif
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_init */

_class_Source_Maxwell_3 *class_Source_Maxwell_3_init(_class_Source_Maxwell_3 *_comp
) {
  #define size (_comp->_parameters.size)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define T1 (_comp->_parameters.T1)
  #define T2 (_comp->_parameters.T2)
  #define T3 (_comp->_parameters.T3)
  #define I1 (_comp->_parameters.I1)
  #define I2 (_comp->_parameters.I2)
  #define I3 (_comp->_parameters.I3)
  #define target_index (_comp->_parameters.target_index)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_source (_comp->_parameters.w_source)
  #define h_source (_comp->_parameters.h_source)
  SIG_MESSAGE("[_source_init] component source=Source_Maxwell_3() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//sources/Source_Maxwell_3.comp:88]");

  if (target_index && !dist)
  {
    Coords ToTarget;
    double tx,ty,tz;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index),POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &tx, &ty, &tz);
    dist=sqrt(tx*tx+ty*ty+tz*tz);
  }

  if (size>0) {
    w_source = h_source = size;
  } else {
    w_source = xwidth;
    h_source = yheight;
  }
  if (lambda0) {
    Lmin=lambda0-dlambda;
    Lmax=lambda0+dlambda;
  }
  l_range = Lmax-Lmin;
  w_mult = w_source*h_source*1.0e4;     /* source area correction */
  w_mult *= l_range;            /* wavelength range correction */
  w_mult *= 1.0/mcget_ncount();   /* correct for # neutron rays */

  if (w_source <0 || h_source < 0 || Lmin <= 0 || Lmax <= 0 || dist <= 0 || T1 <= 0 || T2 <= 0|| T3 <= 0 || Lmax<=Lmin) {
      printf("Source_Maxwell_3: %s: Error in input parameter values!\n"
             "ERROR          Exiting\n",
           NAME_CURRENT_COMP);
      exit(0);
  }

  #undef size
  #undef yheight
  #undef xwidth
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef T1
  #undef T2
  #undef T3
  #undef I1
  #undef I2
  #undef I3
  #undef target_index
  #undef lambda0
  #undef dlambda
  #undef l_range
  #undef w_mult
  #undef w_source
  #undef h_source
  return(_comp);
} /* class_Source_Maxwell_3_init */

_class_L_monitor *class_L_monitor_init(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  SIG_MESSAGE("[_l_monitor_before_mono_init] component l_monitor_before_mono=L_monitor() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:66]");

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("L_monitor: %s: Null detection area !\n"
      "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
      NAME_CURRENT_COMP);
    exit(0);
  }

  L_N = create_darr1d(nL);
  L_p = create_darr1d(nL);
  L_p2 = create_darr1d(nL);

  int i;
  for (i=0; i<nL; i++)
  {
    L_N[i] = 0;
    L_p[i] = 0;
    L_p2[i] = 0;
  }

  // Use instance name for monitor output if no input was given
  if (!strcmp(filename,"\0")) sprintf(filename,NAME_CURRENT_COMP);
  #undef nL
  #undef filename
  #undef nowritefile
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_init */

_class_Monochromator_flat *class_Monochromator_flat_init(_class_Monochromator_flat *_comp
) {
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  SIG_MESSAGE("[_monochromator_flat_init] component monochromator_flat=Monochromator_flat() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Monochromator_flat.comp:103]");

  mos_rms_y = MIN2RAD*mosaicv/sqrt(8*log(2));
  mos_rms_z = MIN2RAD*mosaich/sqrt(8*log(2));
  mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

  mono_Q = Q;
  if (DM != 0) mono_Q = 2*PI/DM;

  if (zwidth>0)  { zmax = zwidth/2;  zmin=-zmax; }
  if (yheight>0) { ymax = yheight/2; ymin=-ymax; }

  if (zmin==zmax || ymin==ymax)
    exit(fprintf(stderr, "Monochromator_flat: %s : Surface is null (zmin,zmax,ymin,ymax)\n", NAME_CURRENT_COMP));
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_init */

_class_Collimator_linear *class_Collimator_linear_init(_class_Collimator_linear *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
  SIG_MESSAGE("[_collimator_linear1_init] component collimator_linear1=Collimator_linear() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Collimator_linear.comp:55]");

slope = tan(MIN2RAD*divergence);
  slopeV= tan(MIN2RAD*divergenceV);
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("Collimator_linear: %s: Null slit opening area !\n"
	         "ERROR              (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_init */

_class_E_monitor *class_E_monitor_init(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_before_sample_init] component e_monitor_before_sample=E_monitor() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:69]");

  int i;

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)) {
    printf("E_monitor: %s: Null detection area !\n"
           "ERROR      (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
           NAME_CURRENT_COMP);
    exit(0);
  }

  E_N = create_darr1d(nE);
  E_p = create_darr1d(nE);
  E_p2 = create_darr1d(nE);

  for (i=0; i<nE; i++)
  {
    E_N[i] = 0;
    E_p[i] = 0;
    E_p2[i] = 0;
  }
  S_p = S_pE = S_pE2 = 0;

  // Use instance name for monitor output if no input was given
  if (!strcmp(filename,"\0")) sprintf(filename,NAME_CURRENT_COMP);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_init */

_class_PSD_monitor *class_PSD_monitor_init(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_PSD_monitor_before_sample_init] component PSD_monitor_before_sample=PSD_monitor() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:62]");

  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(0);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  // Use instance name for monitor output if no input was given
  if (!strcmp(filename,"\0")) sprintf(filename,NAME_CURRENT_COMP);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_init */

_class_Phonon_BvK_PG_eigenvector *class_Phonon_BvK_PG_eigenvector_init(_class_Phonon_BvK_PG_eigenvector *_comp
) {
  #define hh (_comp->_parameters.hh)
  #define kk (_comp->_parameters.kk)
  #define ll (_comp->_parameters.ll)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define DW (_comp->_parameters.DW)
  #define T (_comp->_parameters.T)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define target_index (_comp->_parameters.target_index)
  #define mode_input (_comp->_parameters.mode_input)
  #define e_steps_low (_comp->_parameters.e_steps_low)
  #define e_steps_high (_comp->_parameters.e_steps_high)
  #define verbose_input (_comp->_parameters.verbose_input)
  #define dispersion (_comp->_parameters.dispersion)
  #define V_rho (_comp->_parameters.V_rho)
  #define V_my_s (_comp->_parameters.V_my_s)
  #define V_my_a_v (_comp->_parameters.V_my_a_v)
  #define Matrix (_comp->_parameters.Matrix)
  #define q (_comp->_parameters.q)
  #define qx (_comp->_parameters.qx)
  #define qy (_comp->_parameters.qy)
  #define qz (_comp->_parameters.qz)
  #define q_x (_comp->_parameters.q_x)
  #define q_y (_comp->_parameters.q_y)
  #define q_z (_comp->_parameters.q_z)
  #define eigenvectormode (_comp->_parameters.eigenvectormode)
  #define p_call (_comp->_parameters.p_call)
  SIG_MESSAGE("[_phonon_bvk_pg_init] component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() INITIALISE [Phonon_BvK_PG_eigenvector.comp:964]");

  V_rho = 4/(a_latt*a_latt*c_latt*sqrt3/2);  // (unit: Å^-3)
  V_my_s = (V_rho * 100 * sigma_inc);
  V_my_a_v = (V_rho * 100 * sigma_abs * 2200);

  /* now compute target coords if a component index is supplied */
  if (!target_index && !target_x && !target_y && !target_z) target_index=1;
  if (target_index){
    Coords ToTarget;
    ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+target_index), POS_A_CURRENT_COMP);
    ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
    coords_get(ToTarget, &target_x, &target_y, &target_z);
  }
  if (!(target_x || target_y || target_z)) {
    printf("Phonon_simple: %s: The target is not defined. Using direct beam (Z-axis).\n",
      NAME_CURRENT_COMP);
    target_z=1;
  }
  initialize_omega_q();
  
  /* Now fill global variables for functions */
        verbose = verbose_input;

  #undef hh
  #undef kk
  #undef ll
  #undef radius
  #undef yheight
  #undef sigma_abs
  #undef sigma_inc
  #undef DW
  #undef T
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_x
  #undef target_y
  #undef target_z
  #undef target_index
  #undef mode_input
  #undef e_steps_low
  #undef e_steps_high
  #undef verbose_input
  #undef dispersion
  #undef V_rho
  #undef V_my_s
  #undef V_my_a_v
  #undef Matrix
  #undef q
  #undef qx
  #undef qy
  #undef qz
  #undef q_x
  #undef q_y
  #undef q_z
  #undef eigenvectormode
  #undef p_call
  return(_comp);
} /* class_Phonon_BvK_PG_eigenvector_init */

_class_Slit *class_Slit_init(_class_Slit *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define isradial (_comp->_parameters.isradial)
  SIG_MESSAGE("[_slit1_init] component slit1=Slit() INITIALISE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Slit.comp:66]");

if (is_unset(radius)){
  isradial=0;
  if (all_set(3, xwidth, xmin, xmax)){
    slit_error_if(xwidth != xmax - xmin, "specifying xwidth, xmin and xmax requires consistent parameters", NAME_CURRENT_COMP);
  } else {
    slit_error_if(is_unset(xwidth) && any_unset(2, xmin, xmax), "specify either xwidth or xmin & xmax", NAME_CURRENT_COMP);
  }
  if (all_set(3, yheight, ymin, ymax)){
    slit_error_if(yheight != ymax - ymin, "specifying yheight, ymin and ymax requires consistent parameters", NAME_CURRENT_COMP);
  } else {
    slit_error_if(is_unset(yheight) && any_unset(2, ymin, ymax), "specify either yheight or ymin & ymax", NAME_CURRENT_COMP);
  }
  if (is_unset(xmin)) { // xmax also unset but xwidth *is* set
    xmax = xwidth/2;
    xmin = -xmax;
  }
  if (is_unset(ymin)) { // ymax also unset but yheight *is* set
    ymax = yheight/2;
    ymin = -ymax;
  }
  slit_warning_if(xmin == xmax || ymin == ymax, "Running with CLOSED rectangular slit - is this intentional?", NAME_CURRENT_COMP);
} else {
  isradial=1;
  slit_error_if(any_set(6, xwidth, xmin, xmax, yheight, ymin, ymax), 
                "specify radius OR width and height parameters", NAME_CURRENT_COMP);
  slit_warning_if(radius == 0., "Running with CLOSED radial slit - is this intentional?", NAME_CURRENT_COMP);
}

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  #undef isradial
  return(_comp);
} /* class_Slit_init */



int init(void) { /* called by mccode_main for template_body_simple:INITIALISE */
  DEBUG_INSTR();

  /* code_main/parseoptions/readparams sets instrument parameters value */
  stracpy(instrument->_name, "template_body_simple", 256);

  /* Instrument 'template_body_simple' INITIALISE */
  SIG_MESSAGE("[template_body_simple] INITIALISE [phonon_eigenvector.instr:52]");
  #define E (instrument->_parameters.E)
  #define Ef (instrument->_parameters.Ef)
  #define Dlambda (instrument->_parameters.Dlambda)
  #define h (instrument->_parameters.h)
  #define l (instrument->_parameters.l)
  #define dA3 (instrument->_parameters.dA3)
  #define Temp (instrument->_parameters.Temp)
  #define width (instrument->_parameters.width)
  #define coll (instrument->_parameters.coll)
  #define phononmode (instrument->_parameters.phononmode)
  #define E_steps_high (instrument->_parameters.E_steps_high)
  #define E_steps_low (instrument->_parameters.E_steps_low)
  #define Verbose (instrument->_parameters.Verbose)
  #define DISP (instrument->_parameters.DISP)
{
// Set monochromator/analyzer Q-value for PG
QM = 1.8734;
SMALL__NUMBER = 1e-6;

double a = 2.461;
double c =  6.708;
//PG lattice constants, in AA, same as in Phonon_BvK_PG.comp

double astar = 4*pi/(sqrt(3)*a);
double cstar = 2*pi/c;
//PG reciprocal lattice vectors, in 1/AA

//calculate Ei
Ei=Ef+E;

//calculate ki, kf, lambda_i, q
double ki = V2K*SE2V*sqrt(Ei);
double kf = V2K*SE2V*sqrt(Ef);
lambda_i=2*pi/ki;
double q = sqrt(h*h*astar*astar+l*l*cstar*cstar);

//calculate 2thetaM and 2thetaA
thetaM = asin(QM/(2*ki))*RAD2DEG;
thetaA = asin(QM/(2*kf))*RAD2DEG;

//calculate scattering angle and sample rotation
twothetaS = acos((q*q-ki*ki-kf*kf)/(-2*ki*kf))*RAD2DEG;
alpha = acos((kf*kf-ki*ki-q*q)/(-2*ki*q));
A3=(asin(l*cstar/(q+SMALL__NUMBER))-alpha)*RAD2DEG;

//Printout calculations
printf("a = %g c = %g astar = %g cstar = %g \n",a,c,astar,cstar);
printf("ki = %g kf = %g q = %g lambda_i = %g\n",ki,kf,q,lambda_i);
printf("thetaA = %g thetaM = %g\n",thetaA,thetaM);
printf("twothetaS = %g \n",twothetaS);
printf("A3 = %g \n",A3);
printf("alpha = %g \n",alpha);

}
  #undef E
  #undef Ef
  #undef Dlambda
  #undef h
  #undef l
  #undef dA3
  #undef Temp
  #undef width
  #undef coll
  #undef phononmode
  #undef E_steps_high
  #undef E_steps_low
  #undef Verbose
  #undef DISP
  _origin_setpos(); /* type Progress_bar */
  _source_setpos(); /* type Source_Maxwell_3 */
  _l_monitor_before_mono_setpos(); /* type L_monitor */
  _monochromator_flat_setpos(); /* type Monochromator_flat */
  _arm1_setpos(); /* type Arm */
  _collimator_linear1_setpos(); /* type Collimator_linear */
  _l_monitor_before_sample_setpos(); /* type L_monitor */
  _e_monitor_before_sample_setpos(); /* type E_monitor */
  _PSD_monitor_before_sample_setpos(); /* type PSD_monitor */
  _arm2_setpos(); /* type Arm */
  _phonon_bvk_pg_setpos(); /* type Phonon_BvK_PG_eigenvector */
  _arm3_setpos(); /* type Arm */
  _collimator_linear2_setpos(); /* type Collimator_linear */
  _psd_monitoraftersample_setpos(); /* type PSD_monitor */
  _slit1_setpos(); /* type Slit */
  _e_monitor_beforeana_setpos(); /* type E_monitor */
  _analyzer_setpos(); /* type Monochromator_flat */
  _arm4_setpos(); /* type Arm */
  _slit2_setpos(); /* type Slit */
  _Emon_after_analyzer_setpos(); /* type E_monitor */
  _psd_detector_setpos(); /* type PSD_monitor */

  /* call iteratively all components INITIALISE */
  class_Progress_bar_init(&_origin_var);

  class_Source_Maxwell_3_init(&_source_var);

  class_L_monitor_init(&_l_monitor_before_mono_var);

  class_Monochromator_flat_init(&_monochromator_flat_var);


  class_Collimator_linear_init(&_collimator_linear1_var);

  class_L_monitor_init(&_l_monitor_before_sample_var);

  class_E_monitor_init(&_e_monitor_before_sample_var);

  class_PSD_monitor_init(&_PSD_monitor_before_sample_var);


  class_Phonon_BvK_PG_eigenvector_init(&_phonon_bvk_pg_var);


  class_Collimator_linear_init(&_collimator_linear2_var);

  class_PSD_monitor_init(&_psd_monitoraftersample_var);

  class_Slit_init(&_slit1_var);

  class_E_monitor_init(&_e_monitor_beforeana_var);

  class_Monochromator_flat_init(&_analyzer_var);


  class_Slit_init(&_slit2_var);

  class_E_monitor_init(&_Emon_after_analyzer_var);

  class_PSD_monitor_init(&_psd_detector_var);

  if (mcdotrace) display();
  DEBUG_INSTR_END();

#ifdef OPENACC
#include <openacc.h>
#pragma acc update device(_origin_var)
#pragma acc update device(_source_var)
#pragma acc update device(_l_monitor_before_mono_var)
#pragma acc update device(_monochromator_flat_var)
#pragma acc update device(_arm1_var)
#pragma acc update device(_collimator_linear1_var)
#pragma acc update device(_l_monitor_before_sample_var)
#pragma acc update device(_e_monitor_before_sample_var)
#pragma acc update device(_PSD_monitor_before_sample_var)
#pragma acc update device(_arm2_var)
#pragma acc update device(_phonon_bvk_pg_var)
#pragma acc update device(_arm3_var)
#pragma acc update device(_collimator_linear2_var)
#pragma acc update device(_psd_monitoraftersample_var)
#pragma acc update device(_slit1_var)
#pragma acc update device(_e_monitor_beforeana_var)
#pragma acc update device(_analyzer_var)
#pragma acc update device(_arm4_var)
#pragma acc update device(_slit2_var)
#pragma acc update device(_Emon_after_analyzer_var)
#pragma acc update device(_psd_detector_var)
#pragma acc update device(_instrument_var)
#endif

  return(0);
} /* init */

/*******************************************************************************
* components TRACE
*******************************************************************************/

#define x (_particle->x)
#define y (_particle->y)
#define z (_particle->z)
#define vx (_particle->vx)
#define vy (_particle->vy)
#define vz (_particle->vz)
#define t (_particle->t)
#define sx (_particle->sx)
#define sy (_particle->sy)
#define sz (_particle->sz)
#define p (_particle->p)
#define mcgravitation (_particle->mcgravitation)
#define mcMagnet (_particle->mcMagnet)
#define allow_backprop (_particle->allow_backprop)
/* if on GPU, globally nullify sprintf,fprintf,printfs   */
/* (Similar defines are available in each comp trace but */
/*  those are not enough to handle external libs etc. )  */
#ifdef OPENACC
#ifndef MULTICORE
#define fprintf(stderr,...) printf(__VA_ARGS__)
#define sprintf(string,...) printf(__VA_ARGS__)
#define exit(...) noprintf()
#define strcmp(a,b) str_comp(a,b)
#define strlen(a) str_len(a)
#endif
#endif
#define SCATTERED (_particle->_scattered)
#define RESTORE (_particle->_restore)
#define RESTORE_NEUTRON(_index, ...) _particle->_restore = _index;
#define ABSORBED (_particle->_absorbed)
#define mcget_run_num() _particle->_uid
#define ABSORB0 do { DEBUG_STATE(); DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(_comp); } while(0)
#define ABSORB ABSORB0
#pragma acc routine
_class_Progress_bar *class_Progress_bar_trace(_class_Progress_bar *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_trace] component origin=Progress_bar() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:73]");

#ifndef OPENACC
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) save(NULL);
  }
#endif
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_trace */

#pragma acc routine
_class_Source_Maxwell_3 *class_Source_Maxwell_3_trace(_class_Source_Maxwell_3 *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define size (_comp->_parameters.size)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define T1 (_comp->_parameters.T1)
  #define T2 (_comp->_parameters.T2)
  #define T3 (_comp->_parameters.T3)
  #define I1 (_comp->_parameters.I1)
  #define I2 (_comp->_parameters.I2)
  #define I3 (_comp->_parameters.I3)
  #define target_index (_comp->_parameters.target_index)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_source (_comp->_parameters.w_source)
  #define h_source (_comp->_parameters.h_source)
  SIG_MESSAGE("[_source_trace] component source=Source_Maxwell_3() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//sources/Source_Maxwell_3.comp:124]");

  double v,tau_l,E,lambda,k,r,xf,yf,dx,dy,w_focus;

  t=0;
  z=0;
  x = 0.5*w_source*randpm1();
  y = 0.5*h_source*randpm1();         /* Choose initial position */

  randvec_target_rect_real(&xf, &yf, &r, &w_focus,
		      0, 0, dist, focus_xw, focus_yh, ROT_A_CURRENT_COMP, x, y, z, 2);

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);

  lambda = Lmin+l_range*rand01();    /* Choose from uniform distribution */
  k = 2*PI/lambda;
  v = K2V*k;

  vz = v*dist/r;
  vy = v*dy/r;
  vx = v*dx/r;


/*  printf("pos0 (%g %g %g), pos1 (%g %g %g), r: %g, v (%g %g %g), v %g\n",
  x,y,z,xf,yf,dist,r,vx,vy,vz, v);
  printf("l %g, w_focus %g \n", lambda, w_focus);  */

  p *= w_mult*w_focus;                /* Correct for target focusing etc */
  p *= I1*SM3_Maxwell(lambda,T1)+I2*SM3_Maxwell(lambda,T2)+I3*SM3_Maxwell(lambda,T3);
                                        /* Calculate true intensity */
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef size
  #undef yheight
  #undef xwidth
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef T1
  #undef T2
  #undef T3
  #undef I1
  #undef I2
  #undef I3
  #undef target_index
  #undef lambda0
  #undef dlambda
  #undef l_range
  #undef w_mult
  #undef w_source
  #undef h_source
  return(_comp);
} /* class_Source_Maxwell_3_trace */

#pragma acc routine
_class_L_monitor *class_L_monitor_trace(_class_L_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  SIG_MESSAGE("[_l_monitor_before_mono_trace] component l_monitor_before_mono=L_monitor() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:94]");

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    double L = (2*PI/V2K)/sqrt(vx*vx + vy*vy + vz*vz);
    int i = floor((L-Lmin)*nL/(Lmax-Lmin));
    if(i >= 0 && i < nL)
    {
      double p2 = p*p;
      #pragma acc atomic
      L_N[i] = L_N[i] +1;
      #pragma acc atomic
      L_p[i] = L_p[i] + p;
      #pragma acc atomic
      L_p2[i] = L_p2[i] + p2;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef nL
  #undef filename
  #undef nowritefile
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_trace */

#pragma acc routine
_class_Monochromator_flat *class_Monochromator_flat_trace(_class_Monochromator_flat *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  SIG_MESSAGE("[_monochromator_flat_trace] component monochromator_flat=Monochromator_flat() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Monochromator_flat.comp:119]");

  double y1,z1,t1,dt,kix,kiy,kiz,ratio,order,q0x,k,q0,theta;
  double bx,by,bz,kux,kuy,kuz,ax,ay,az,phi;
  double cos_2theta,k_sin_2theta,cos_phi,sin_phi,q_x,q_y,q_z;
  double delta,p_reflect,total,c1x,c1y,c1z,width,mos_sample;
  int i;

  if(vx != 0.0 && (dt = -x/vx) >= 0.0)
  {                             /* Moving towards crystal? */
    y1 = y + vy*dt;             /* Propagate to crystal plane */
    z1 = z + vz*dt;
    t1 = t + dt;
    if (z1>zmin && z1<zmax && y1>ymin && y1<ymax)
    {                           /* Intersect the crystal? */
      kix = V2K*vx;             /* Initial wave vector */
      kiy = V2K*vy;
      kiz = V2K*vz;
      /* Get reflection order and corresponding nominal scattering vector q0
         of correct length and direction. Only the order with the closest
         scattering vector is considered */
      ratio = -2*kix/mono_Q;
      order = floor(ratio + .5);
      if(order == 0.0)
        order = ratio < 0 ? -1 : 1;
      /* Order will be negative when the neutron enters from the back, in
         which case the direction of Q0 is flipped. */
      if(order < 0)
        order = -order;
      /* Make sure the order is small enough to allow Bragg scattering at the
         given neutron wavelength */
      k = sqrt(kix*kix + kiy*kiy + kiz*kiz);
      kux = kix/k;              /* Unit vector along ki */
      kuy = kiy/k;
      kuz = kiz/k;
      if(order > 2*k/mono_Q)
        order--;
      if(order > 0)             /* Bragg scattering possible? */
      {
        q0 = order*mono_Q;
        q0x = ratio < 0 ? -q0 : q0;
        theta = asin(q0/(2*k)); /* Actual bragg angle */
        /* Make MC choice: reflect or transmit? */
        delta = asin(fabs(kux)) - theta;
        p_reflect = r0*exp(-kiy*kiy/(kiy*kiy + kiz*kiz)*(delta*delta)/
                           (2*mos_rms_y*mos_rms_y))*
                       exp(-kiz*kiz/(kiy*kiy + kiz*kiz)*(delta*delta)/
                           (2*mos_rms_z*mos_rms_z));
        if(rand01() < p_reflect)
        {                       /* Reflect */
          cos_2theta = cos(2*theta);
          k_sin_2theta = k*sin(2*theta);
          /* Get unit normal to plane containing ki and most probable kf */
          vec_prod(bx, by, bz, kix, kiy, kiz, q0x, 0, 0);
          NORM(bx,by,bz);
          bx *= k_sin_2theta;
          by *= k_sin_2theta;
          bz *= k_sin_2theta;
          /* Get unit vector normal to ki and b */
          vec_prod(ax, ay, az, bx, by, bz, kux, kuy, kuz);
          /* Compute the total scattering probability at this ki */
          total = 0;
          /* Choose width of Gaussian distribution to sample the angle
           * phi on the Debye-Scherrer cone for the scattered neutron.
           * The radius of the Debye-Scherrer cone is smaller by a
           * factor 1/cos(theta) than the radius of the (partial) sphere
           * describing the possible orientations of Q due to mosaicity, so we
           * start with a width 1/cos(theta) greater than the largest of
           * the two mosaics. */
          mos_sample = mos_rms_max/cos(theta);
          c1x = kix*(cos_2theta-1);
          c1y = kiy*(cos_2theta-1);
          c1z = kiz*(cos_2theta-1);
          /* Loop, repeatedly reducing the sample width until it is small
           * enough to avoid sampling scattering directions with
           * ridiculously low scattering probability.
           * Use a cut-off at 5 times the gauss width for considering
           * scattering probability as well as for integration limits
           * when integrating the sampled distribution below. */
          for(i=0; i<100; i++) {
            width = 5*mos_sample;
            cos_phi = cos(width);
            sin_phi = sin(width);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = (c1y + cos_phi*ay + sin_phi*by)/mos_rms_y;
            q_z = (c1z + cos_phi*az + sin_phi*bz)/mos_rms_z;
            /* Stop when we get near a factor of 25=5^2. */
            if(q_z*q_z + q_y*q_y < (25/(2.0/3.0))*(q_x*q_x))
              break;
            mos_sample *= (2.0/3.0);
          }
          /* Now integrate the chosen sampling distribution, using a
           * cut-off at five times sigma. */
          for(i = 0; i < (sizeof(Gauss_X)/sizeof(double)); i++)
          {
            phi = width*Gauss_X[i];
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            q_x = c1x + cos_phi*ax + sin_phi*bx;
            q_y = c1y + cos_phi*ay + sin_phi*by;
            q_z = c1z + cos_phi*az + sin_phi*bz;
            p_reflect = GAUSS((q_y/q_x),0,mos_rms_y)*
                        GAUSS((q_z/q_x),0,mos_rms_z);
            total += Gauss_W[i]*p_reflect;
          }
          total *= width;
          /* Choose point on Debye-Scherrer cone. Sample from a Gaussian of
           * width 1/cos(theta) greater than the mosaic and correct for any
           * error by adjusting the neutron weight later. */
          phi = mos_sample*randnorm();
          /* Compute final wave vector kf and scattering vector q = ki - kf */
          cos_phi = cos(phi);
          sin_phi = sin(phi);
          q_x = c1x + cos_phi*ax + sin_phi*bx;
          q_y = c1y + cos_phi*ay + sin_phi*by;
          q_z = c1z + cos_phi*az + sin_phi*bz;
          p_reflect = GAUSS((q_y/q_x),0,mos_rms_y)*
                      GAUSS((q_z/q_x),0,mos_rms_z);
          x = 0;
          y = y1;
          z = z1;
          t = t1;
          vx = K2V*(kix+q_x);
          vy = K2V*(kiy+q_y);
          vz = K2V*(kiz+q_z);
          p_reflect /= total*GAUSS(phi,0,mos_sample);
          if (p_reflect <= 0) ABSORB;
          if (p_reflect > 1)  p_reflect = 1;
          p *= p_reflect;
          SCATTER;
        } /* End MC choice to reflect or transmit neutron */
      } /* End bragg scattering possible */
    } /* End intersect the crystal */
  } /* End neutron moving towards crystal */
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_trace */

#pragma acc routine
_class_Collimator_linear *class_Collimator_linear_trace(_class_Collimator_linear *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
  SIG_MESSAGE("[_collimator_linear1_trace] component collimator_linear1=Collimator_linear() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Collimator_linear.comp:70]");

    double phi, dt;

    PROP_Z0;
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;
    dt = length/vz;
    PROP_DT(dt);
    if (x<xmin || x>xmax || y<ymin || y>ymax)
      ABSORB;

    if(slope > 0.0)
    {
      phi = fabs(vx/vz);
      if (phi > slope)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slope);
      SCATTER;
    }
    if (slopeV > 0) {
      phi = fabs(vy/vz);
      if (phi > slopeV)
        ABSORB;
      else
        p *= transmission*(1.0 - phi/slopeV);
      SCATTER;
    }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_trace */

#pragma acc routine
_class_E_monitor *class_E_monitor_trace(_class_E_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_before_sample_trace] component e_monitor_before_sample=E_monitor() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:99]");

  int i;
  double E;

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax)
  {
    E = VS2E*(vx*vx + vy*vy + vz*vz);

    S_p += p;
    S_pE += p*E;
    S_pE2 += p*E*E;

    i = floor((E-Emin)*nE/(Emax-Emin));
    if(i >= 0 && i < nE)
    {
      double p2 = p*p;
      #pragma acc atomic
      E_N[i] = E_N[i] +1;
      #pragma acc atomic
      E_p[i] = E_p[i] + p;
      #pragma acc atomic
      E_p2[i] = E_p2[i] + p2;
      SCATTER;
    }
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_trace */

#pragma acc routine
_class_PSD_monitor *class_PSD_monitor_trace(_class_PSD_monitor *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_PSD_monitor_before_sample_trace] component PSD_monitor_before_sample=PSD_monitor() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:82]");

  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));

    double p2 = p*p;
    #pragma acc atomic
    PSD_N[i][j] = PSD_N[i][j]+1;

    #pragma acc atomic
    PSD_p[i][j] = PSD_p[i][j]+p;
    
    #pragma acc atomic
    PSD_p2[i][j] = PSD_p2[i][j] + p2;
    
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_trace */

#pragma acc routine
_class_Phonon_BvK_PG_eigenvector *class_Phonon_BvK_PG_eigenvector_trace(_class_Phonon_BvK_PG_eigenvector *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define hh (_comp->_parameters.hh)
  #define kk (_comp->_parameters.kk)
  #define ll (_comp->_parameters.ll)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define DW (_comp->_parameters.DW)
  #define T (_comp->_parameters.T)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define target_index (_comp->_parameters.target_index)
  #define mode_input (_comp->_parameters.mode_input)
  #define e_steps_low (_comp->_parameters.e_steps_low)
  #define e_steps_high (_comp->_parameters.e_steps_high)
  #define verbose_input (_comp->_parameters.verbose_input)
  #define dispersion (_comp->_parameters.dispersion)
  #define V_rho (_comp->_parameters.V_rho)
  #define V_my_s (_comp->_parameters.V_my_s)
  #define V_my_a_v (_comp->_parameters.V_my_a_v)
  #define Matrix (_comp->_parameters.Matrix)
  #define q (_comp->_parameters.q)
  #define qx (_comp->_parameters.qx)
  #define qy (_comp->_parameters.qy)
  #define qz (_comp->_parameters.qz)
  #define q_x (_comp->_parameters.q_x)
  #define q_y (_comp->_parameters.q_y)
  #define q_z (_comp->_parameters.q_z)
  #define eigenvectormode (_comp->_parameters.eigenvectormode)
  #define p_call (_comp->_parameters.p_call)
  SIG_MESSAGE("[_phonon_bvk_pg_trace] component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() TRACE [Phonon_BvK_PG_eigenvector.comp:989]");

  double t0, t1;            /* Entry/exit time for cylinder */
  double v_i, v_f;          /* Neutron velocities: initial, final */
  double vx_i, vy_i, vz_i;  /* Neutron initial velocity vector */
  double dt0, dt;           /* Flight times through sample */
  double l_full;            /* Flight path length for non-scattered neutron */
  double l_i, l_o;          /* Flight path lenght in/out for scattered neutron */
  double my_a_i;            /* Initial attenuation factor */
  double my_a_f;            /* Final attenuation factor */
  double solid_angle;       /* Solid angle of target as seen from scattering point */
  double aim_x=0, aim_y=0, aim_z=1;   /* Position of target relative to scattering point */
  double omega;             /* Energy transfer */
  double qsquare;           /* Square of the scattering vector */
  double bose_factor;       /* Calculated value of the Bose factor */
  int nf, index;            /* Number of allowed final velocities */
  double vf_list[20];       /* List of allowed final velocities */
  double J_factor;          /* Jacobian from delta fnc.s in cross section */
  double f1, f2;            /* Probed values of omega_q minus omega */
  double p1,p2,p3,p4,p5;    /* Temporary multipliers */
  int i,j;
  
  if(cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, radius, yheight))
  {
    //    fprintf(stderr,"intersect\n");      
    if(t0 < 0)
      ABSORB; /* Neutron came from the sample or begins inside */

    if (verbose>=2) 
        printf("neutron entered Phonon_BvK_PG \n");
    /* Neutron enters at t=t0. */
    dt0 = t1-t0;                /* Time in sample */
    v_i = sqrt(vx*vx + vy*vy + vz*vz);
    l_full = v_i * dt0;   /* Length of path through sample if not scattered */
    dt = rand01()*dt0;    /* Time of scattering (relative to t0) */
    l_i = v_i*dt;                 /* Penetration in sample at scattering */
    vx_i=vx;
    vy_i=vy;
    vz_i=vz;
    PROP_DT(dt+t0);             /* Point of scattering */

    aim_x = target_x-x;         /* Vector pointing at target (e.g. analyzer) */
    aim_y = target_y-y;
    aim_z = target_z-z;

    if(focus_aw && focus_ah) {
      randvec_target_rect_angular(&vx, &vy, &vz, &solid_angle,
        aim_x, aim_y, aim_z, focus_aw, focus_ah, ROT_A_CURRENT_COMP);
    } else if(focus_xw && focus_yh) {
      randvec_target_rect(&vx, &vy, &vz, &solid_angle,
        aim_x, aim_y, aim_z, focus_xw, focus_yh, ROT_A_CURRENT_COMP);
    } else {
      randvec_target_sphere(&vx,&vy,&vz,&solid_angle,aim_x,aim_y,aim_z, focus_r);
    }
    NORM(vx, vy, vz);
    nf=0;
    if (mode_input == 12)
    {
        mode = round(12*rand01()-0.5); 
    }
    else
    {
        mode = mode_input;
    }
    if ((mode < 0) || (mode > 11))
    {
        printf("mode = %d ",mode);
        nrerror("illegal value of mode");
    }
    
      p_call[0]=-1;
      p_call[1]=v_i;
      p_call[2]=vx;
      p_call[3]=vy;
      p_call[4]=vz;
      p_call[5]=vx_i;
      p_call[6]=vy_i;
      p_call[7]=vz_i;
      p_call[9]=hh;
      p_call[10]=kk;
      p_call[11]=ll;
      p_call[12] = (double)dispersion;
      
      if (dispersion==1)
      {
          f1=omega_q(p_call,eigenvectormode);
          p_call[12]=0.0;
      }
      
    #ifndef OPENACC
    if (verbose >=2) 
        printf("Call findroots \n");      
    findroots(0, v_i, v_i+V_HIGH, vf_list, &nf, omega_q, p_call, eigenvectormode, e_steps_low, e_steps_high);
    if (verbose >=2)
        printf( "Findroots returned %d roots \n",nf);
    #else
    findroots_gpu(0, v_i, v_i+V_HIGH, vf_list, &nf, p_call, eigenvectormode, e_steps_low, e_steps_high);
    #endif

    if (verbose>=3)
    {
        printf("Findroots done, mode %i , last phonon energy %g \n",mode,real(p_call[8]));  
    }
    
    // ________________ ALL FROM HERE IS CALCULATION OF INTENSITY - SHOULD BE REVISITED ___________________
   
        index=(int)floor(rand01()*nf); // Select random solution
        v_f=vf_list[index];
        p_call[0]=v_f;      // transfer choice of v_f to omega_q

        f1=omega_q(p_call,eigenvectormode);
        if (verbose >= 2)
        {
            printf("Chosen solution is %i; v_f = %g, hw = %g, energy match is %g; eigenvector is: (", index, v_f, real(p_call[8]), f1);
            for (j=0; j<DIM; j++)
                printf(" %g +i %g, ",real(eigenvectormode[j]),imag(eigenvectormode[j]));
            printf(" )\n");
            printf("--------------------------------------------\n");
        }
      p_call[0]=v_f-DV;
      //      fprintf(stderr,"f1=omega_q(%g+i%g)\n",creal(p_call[0]),cimag(p_call[0]));            
      f1=omega_q(p_call,eigenvectormode);
      //      fprintf(stderr,"retur fra f1=omega_q\n");            
      p_call[0]=v_f+DV;
      //      fprintf(stderr,"f2=omega_q(%g+i%g)\n",creal(p_call[0]),cimag(p_call[0]));                  
      f2=omega_q(p_call,eigenvectormode);
      //      fprintf(stderr,"retur fra f2=omega_q\n");                  
      J_factor = fabs(f2-f1)/(2*DV); // Removed K2V, double check!
      omega=VS2E*(v_i*v_i-v_f*v_f);
      vx *= v_f;
      vy *= v_f;
      vz *= v_f;

      q_x=q[0]=V2K*(vx_i-vx);
      q_y=q[1]=V2K*(vy_i-vy);
      q_z=q[2]=V2K*(vz_i-vz);

      if (dispersion==1)
      {
          printf("Dispersion found:  (h,k,l) = (%g, %g, %g), hw = %g \n",q[0]/astar,q[1]/astar,q[2]/cstar,omega);
      }
      
      qsquare=1;
      
      komplex Fn = 0;
      komplex q_dot_e;
      double q_dot_r;
      for (i=0; i<4; i++) // loop over atoms in unit cell
      {
          q_dot_e = (komplex)0;
          q_dot_r = 0;
          for (j=0; j<3; j++) // loop over x,y,z
          {
              q_dot_e += (komplex)q[j]*eigenvectormode[3*i+j]; 
              q_dot_r += q[j]*Delta[3*i+j];
          }
          if (verbose >= 7)
          {
              printf("i = %i , q dot r = %g , q dot e = (%g + i %g) Fn contrib = (%g + i %g)",i,q_dot_r,real(q_dot_e),imag(q_dot_e),real(b_length*q_dot_e*eksp(I*q_dot_r)),imag(b_length*q_dot_e*eksp(I*q_dot_r)));
          }
          Fn += b_length*q_dot_e*eksp(I*q_dot_r);
      }
      if (verbose >= 7) 
          printf("\n");
      
      if(!cylinder_intersect(&t0, &t1, x, y, z, vx, vy, vz, radius, yheight))
	{
	  printf("FATAL ERROR: Did not hit cylinder from inside.\n");
	  exit(1);
	}
	
	if (verbose>=2)
        printf("neutron left Phonon_BvK_PG, p_init = %g,",p);
      dt = t1;
      l_o = v_f*dt;
      
      my_a_i = V_my_a_v/v_i;
      my_a_f = V_my_a_v/v_f;
      bose_factor=nbose(omega,T);
      p1 = exp(-(V_my_s*(l_i+l_o)+my_a_i*l_i+my_a_f*l_o)); /* Absorption factor (units: 1) */
      p2 = nf*solid_angle*l_full*V_rho/(4*PI);     /* Focusing factors; assume random choice of n_f possibilities (units: m/Å^3) */
      p3 = (v_f/v_i)*DW/fabs(omega)*bose_factor;   /* Cross section factor 1 (units: s) */
      p4 = 2*VS2E*v_f/J_factor;  /* Jacobian of delta functions in cross section */
      p5 = Fn*conj(Fn)/M;  /* Cross section factor 2. !!! TEST !!!  (units: fm^2 Å^2) */
      p *= p1*p2*p3*p4*p5;
      if (mode_input == 12)
          p *= 12;  //The factor 12 is due to MC choice between the 12 modes

      if (verbose>=2) 
        printf(" p1 = %g, p2= %g, p3 = %g, p4 = %g, p5 = %g, bose_factor = %g, p_final = %g\n",p1,p2,p3,p4,p5,bose_factor,p);
  } /* else transmit: Neutron did not hit the sample */
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif
  #undef hh
  #undef kk
  #undef ll
  #undef radius
  #undef yheight
  #undef sigma_abs
  #undef sigma_inc
  #undef DW
  #undef T
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_x
  #undef target_y
  #undef target_z
  #undef target_index
  #undef mode_input
  #undef e_steps_low
  #undef e_steps_high
  #undef verbose_input
  #undef dispersion
  #undef V_rho
  #undef V_my_s
  #undef V_my_a_v
  #undef Matrix
  #undef q
  #undef qx
  #undef qy
  #undef qz
  #undef q_x
  #undef q_y
  #undef q_z
  #undef eigenvectormode
  #undef p_call
  return(_comp);
} /* class_Phonon_BvK_PG_eigenvector_trace */

#pragma acc routine
_class_Slit *class_Slit_trace(_class_Slit *_comp
  , _class_particle *_particle) {
  ABSORBED=SCATTERED=RESTORE=0;
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define isradial (_comp->_parameters.isradial)
  SIG_MESSAGE("[_slit1_trace] component slit1=Slit() TRACE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Slit.comp:98]");

    PROP_Z0;
    if (!isradial ? (x < xmin || x > xmax || y < ymin || y > ymax) : (x * x + y * y > radius * radius))
      ABSORB;
    else
      SCATTER;
#ifndef NOABSORB_INF_NAN
  /* Check for nan or inf particle parms */ 
  if(isnan(p)  ||  isinf(p)) ABSORB;
  if(isnan(t)  ||  isinf(t)) ABSORB;
  if(isnan(vx) || isinf(vx)) ABSORB;
  if(isnan(vy) || isinf(vy)) ABSORB;
  if(isnan(vz) || isinf(vz)) ABSORB;
  if(isnan(x)  ||  isinf(x)) ABSORB;
  if(isnan(y)  ||  isinf(y)) ABSORB;
  if(isnan(z)  ||  isinf(z)) ABSORB;
#else
  if(isnan(p)  ||  isinf(p)) printf("NAN or INF found in p,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(t)  ||  isinf(t)) printf("NAN or INF found in t,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vx) || isinf(vx)) printf("NAN or INF found in vx, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vy) || isinf(vy)) printf("NAN or INF found in vy, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(vz) || isinf(vz)) printf("NAN or INF found in vz, %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(x)  ||  isinf(x)) printf("NAN or INF found in x,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(y)  ||  isinf(y)) printf("NAN or INF found in y,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
  if(isnan(z)  ||  isinf(z)) printf("NAN or INF found in z,  %s (particle %lld)\n",_comp->_name,_particle->_uid);
#endif

if (_comp->_index == 19) { // EXTEND 'slit2'
//    qx = MC_GETPAR(phonon_bvk_pg,q_x);
//    qy = MC_GETPAR(phonon_bvk_bg,q_y);
//    qz = MC_GETPAR(phonon_bvk_bg,q_z);
//    printf("hit; q = (%g, %g, %g)\n",0,0,0);
}

  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  #undef isradial
  return(_comp);
} /* class_Slit_trace */

/* *****************************************************************************
* instrument 'template_body_simple' TRACE
***************************************************************************** */

#ifndef FUNNEL
#pragma acc routine
int raytrace(_class_particle* _particle) { /* single event propagation, called by mccode_main for template_body_simple:TRACE */

  /* init variables and counters for TRACE */
  #undef ABSORB0
  #undef ABSORB
  #define ABSORB0 do { DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; return(ABSORBED);} while(0)
  #define ABSORB ABSORB0
  DEBUG_ENTER();
  DEBUG_STATE();
  _particle->flag_nocoordschange=0; /* Init */
  _class_particle _particle_save;
  /* the main iteration loop for one incoming event */
  while (!ABSORBED) { /* iterate event until absorbed */
    /* send particle event to component instance, one after the other */
    /* begin component origin=Progress_bar() [1] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_origin_var._rotation_is_identity) {
        if(!_origin_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _origin_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_origin_var._position_relative, _origin_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 1) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_origin_var._name);
      DEBUG_STATE();
      class_Progress_bar_trace(&_origin_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component origin [1] */
    /* begin component source=Source_Maxwell_3() [2] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_source_var._rotation_is_identity) {
        if(!_source_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _source_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_source_var._position_relative, _source_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 2) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_source_var._name);
      DEBUG_STATE();
      class_Source_Maxwell_3_trace(&_source_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component source [2] */
    /* begin component l_monitor_before_mono=L_monitor() [3] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_l_monitor_before_mono_var._rotation_is_identity) {
        if(!_l_monitor_before_mono_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _l_monitor_before_mono_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_l_monitor_before_mono_var._position_relative, _l_monitor_before_mono_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 3) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_l_monitor_before_mono_var._name);
      DEBUG_STATE();
      class_L_monitor_trace(&_l_monitor_before_mono_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component l_monitor_before_mono [3] */
    /* begin component monochromator_flat=Monochromator_flat() [4] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_monochromator_flat_var._rotation_is_identity) {
        if(!_monochromator_flat_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _monochromator_flat_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_monochromator_flat_var._position_relative, _monochromator_flat_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 4) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_monochromator_flat_var._name);
      DEBUG_STATE();
      class_Monochromator_flat_trace(&_monochromator_flat_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component monochromator_flat [4] */
    /* begin component arm1=Arm() [5] */
    if (!ABSORBED && _particle->_index == 5) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component arm1 [5] */
    /* begin component collimator_linear1=Collimator_linear() [6] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_collimator_linear1_var._rotation_is_identity) {
        if(!_collimator_linear1_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _collimator_linear1_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_collimator_linear1_var._position_relative, _collimator_linear1_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 6) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_collimator_linear1_var._name);
      DEBUG_STATE();
      class_Collimator_linear_trace(&_collimator_linear1_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component collimator_linear1 [6] */
    /* begin component l_monitor_before_sample=L_monitor() [7] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_l_monitor_before_sample_var._rotation_is_identity) {
        if(!_l_monitor_before_sample_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _l_monitor_before_sample_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_l_monitor_before_sample_var._position_relative, _l_monitor_before_sample_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 7) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_l_monitor_before_sample_var._name);
      DEBUG_STATE();
      class_L_monitor_trace(&_l_monitor_before_sample_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component l_monitor_before_sample [7] */
    /* begin component e_monitor_before_sample=E_monitor() [8] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_e_monitor_before_sample_var._rotation_is_identity) {
        if(!_e_monitor_before_sample_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_before_sample_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_e_monitor_before_sample_var._position_relative, _e_monitor_before_sample_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 8) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_e_monitor_before_sample_var._name);
      DEBUG_STATE();
      class_E_monitor_trace(&_e_monitor_before_sample_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component e_monitor_before_sample [8] */
    /* begin component PSD_monitor_before_sample=PSD_monitor() [9] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_PSD_monitor_before_sample_var._rotation_is_identity) {
        if(!_PSD_monitor_before_sample_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _PSD_monitor_before_sample_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_PSD_monitor_before_sample_var._position_relative, _PSD_monitor_before_sample_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 9) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_PSD_monitor_before_sample_var._name);
      DEBUG_STATE();
      class_PSD_monitor_trace(&_PSD_monitor_before_sample_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component PSD_monitor_before_sample [9] */
    /* begin component arm2=Arm() [10] */
    if (!ABSORBED && _particle->_index == 10) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component arm2 [10] */
#ifndef NOSPLIT
    /* start SPLIT at phonon_bvk_pg */
    if (!ABSORBED) {
    _class_particle Split_phonon_bvk_pg_particle=*_particle;
    int Split_phonon_bvk_pg_counter;
    int SplitS_phonon_bvk_pg = 1;
    for (Split_phonon_bvk_pg_counter = 0; Split_phonon_bvk_pg_counter< SplitS_phonon_bvk_pg; Split_phonon_bvk_pg_counter++) {
      randstate_t randbackup = *_particle->randstate;
      *_particle=Split_phonon_bvk_pg_particle;
      *_particle->randstate = randbackup;
      srandom(_hash(_particle->_uid* 11 *(1+Split_phonon_bvk_pg_counter)));
      p /= SplitS_phonon_bvk_pg > 0 ? SplitS_phonon_bvk_pg : 1;
#endif
    /* begin component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() [11] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_phonon_bvk_pg_var._rotation_is_identity) {
        if(!_phonon_bvk_pg_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _phonon_bvk_pg_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_phonon_bvk_pg_var._position_relative, _phonon_bvk_pg_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 11) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_phonon_bvk_pg_var._name);
      DEBUG_STATE();
      class_Phonon_BvK_PG_eigenvector_trace(&_phonon_bvk_pg_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component phonon_bvk_pg [11] */
    /* begin component arm3=Arm() [12] */
    if (!ABSORBED && _particle->_index == 12) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component arm3 [12] */
    /* begin component collimator_linear2=Collimator_linear() [13] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_collimator_linear2_var._rotation_is_identity) {
        if(!_collimator_linear2_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _collimator_linear2_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_collimator_linear2_var._position_relative, _collimator_linear2_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 13) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_collimator_linear2_var._name);
      DEBUG_STATE();
      class_Collimator_linear_trace(&_collimator_linear2_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component collimator_linear2 [13] */
    /* begin component psd_monitoraftersample=PSD_monitor() [14] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_psd_monitoraftersample_var._rotation_is_identity) {
        if(!_psd_monitoraftersample_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _psd_monitoraftersample_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_psd_monitoraftersample_var._position_relative, _psd_monitoraftersample_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 14) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_psd_monitoraftersample_var._name);
      DEBUG_STATE();
      class_PSD_monitor_trace(&_psd_monitoraftersample_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component psd_monitoraftersample [14] */
    /* begin component slit1=Slit() [15] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_slit1_var._rotation_is_identity) {
        if(!_slit1_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _slit1_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_slit1_var._position_relative, _slit1_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 15) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_slit1_var._name);
      DEBUG_STATE();
      class_Slit_trace(&_slit1_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component slit1 [15] */
    /* begin component e_monitor_beforeana=E_monitor() [16] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_e_monitor_beforeana_var._rotation_is_identity) {
        if(!_e_monitor_beforeana_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_beforeana_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_e_monitor_beforeana_var._position_relative, _e_monitor_beforeana_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 16) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_e_monitor_beforeana_var._name);
      DEBUG_STATE();
      class_E_monitor_trace(&_e_monitor_beforeana_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component e_monitor_beforeana [16] */
    /* begin component analyzer=Monochromator_flat() [17] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_analyzer_var._rotation_is_identity) {
        if(!_analyzer_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _analyzer_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_analyzer_var._position_relative, _analyzer_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 17) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_analyzer_var._name);
      DEBUG_STATE();
      class_Monochromator_flat_trace(&_analyzer_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component analyzer [17] */
    /* begin component arm4=Arm() [18] */
    if (!ABSORBED && _particle->_index == 18) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle->_index++;
    } /* end component arm4 [18] */
    /* begin component slit2=Slit() [19] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_slit2_var._rotation_is_identity) {
        if(!_slit2_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _slit2_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_slit2_var._position_relative, _slit2_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 19) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_slit2_var._name);
      DEBUG_STATE();
      class_Slit_trace(&_slit2_var, _particle); /* contains EXTEND code */
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component slit2 [19] */
    /* begin component Emon_after_analyzer=E_monitor() [20] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_Emon_after_analyzer_var._rotation_is_identity) {
        if(!_Emon_after_analyzer_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _Emon_after_analyzer_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_Emon_after_analyzer_var._position_relative, _Emon_after_analyzer_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 20) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_Emon_after_analyzer_var._name);
      DEBUG_STATE();
      class_E_monitor_trace(&_Emon_after_analyzer_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component Emon_after_analyzer [20] */
    /* begin component psd_detector=PSD_monitor() [21] */
    if (!_particle->flag_nocoordschange) { // flag activated by JUMP to pass coords change
      if (_psd_detector_var._rotation_is_identity) {
        if(!_psd_detector_var._position_relative_is_zero) {
          coords_get(coords_add(coords_set(x,y,z), _psd_detector_var._position_relative),&x, &y, &z);
        }
      } else {
          mccoordschange(_psd_detector_var._position_relative, _psd_detector_var._rotation_relative, _particle);
      }
    }
    if (!ABSORBED && _particle->_index == 21) {
      _particle->flag_nocoordschange=0; /* Reset if we came here from a JUMP */
      _particle_save = *_particle;
      DEBUG_COMP(_psd_detector_var._name);
      DEBUG_STATE();
      class_PSD_monitor_trace(&_psd_detector_var, _particle);
      if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
      _particle->_index++;
      if (!ABSORBED) { DEBUG_STATE(); }
    } /* end component psd_detector [21] */
#ifndef NOSPLIT
    } /* end SPLIT at phonon_bvk_pg */
    } /* if (!ABSORBED) relating to SPLIT at phonon_bvk_pg */
#endif
    if (_particle->_index > 21)
      ABSORBED++; /* absorbed when passed all components */
  } /* while !ABSORBED */

  DEBUG_LEAVE()
  particle_restore(_particle, &_particle_save);
  DEBUG_STATE()

  return(_particle->_index);
} /* raytrace */

/* loop to generate events and call raytrace() propagate them */
void raytrace_all(unsigned long long ncount, unsigned long seed) {

  /* CPU-loop */
  unsigned long long loops;
  loops = ceil((double)ncount/gpu_innerloop);
  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif

  #ifdef OPENACC
  if (ncount>gpu_innerloop) {
    printf("Defining %llu CPU loops around GPU kernel and adjusting ncount\n",loops);
    mcset_ncount(loops*gpu_innerloop);
  } else {
    #endif
    loops=1;
    gpu_innerloop = ncount;
    #ifdef OPENACC
  }
    #endif

  for (unsigned long long cloop=0; cloop<loops; cloop++) {
    #ifdef OPENACC
    if (loops>1) fprintf(stdout, "%d..", (int)cloop); fflush(stdout);
    #endif

    /* if on GPU, re-nullify printf */
    #ifdef OPENACC
    #ifndef MULTICORE
    #define printf(...) noprintf()
    #endif
    #endif

    #pragma acc parallel loop num_gangs(numgangs) vector_length(vecsize)
    for (unsigned long pidx=0 ; pidx < gpu_innerloop ; pidx++) {
      _class_particle particleN = mcgenstate(); // initial particle
      _class_particle* _particle = &particleN;
      particleN._uid = pidx;

      srandom(_hash((pidx+1)*(seed+1)));
      particle_uservar_init(_particle);

      raytrace(_particle);
    } /* inner for */
    seed = seed+gpu_innerloop;
    mcncount += gpu_innerloop;
  } /* CPU for */
  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif
  MPI_MASTER(
  printf("*** TRACE end *** \n");
  );
} /* raytrace_all */

#endif //no-FUNNEL

#ifdef FUNNEL
// Alternative raytrace algorithm which iterates all particles through
// one component at the time, can remove absorbs from the next loop and
// switch between cpu/gpu.
void raytrace_all_funnel(unsigned long long ncount, unsigned long seed) {

  // set up outer (CPU) loop / particle batches
  unsigned long long loops;

  /* if on GPU, printf has been globally nullified, re-enable here */
  #ifdef OPENACC
  #ifndef MULTICORE
  #undef printf
  #endif
  #endif

  #ifdef OPENACC
  loops = ceil((double)ncount/gpu_innerloop);
  if (ncount>gpu_innerloop) {
    printf("Defining %llu CPU loops around kernel and adjusting ncount\n",loops);
    mcset_ncount(loops*gpu_innerloop);
  } else {
  #endif
    loops=1;
    gpu_innerloop = ncount;
  #ifdef OPENACC
  }
  #endif

  // create particles struct and pointer arrays (same memory used by all batches)
  _class_particle* particles = malloc(gpu_innerloop*sizeof(_class_particle));
  _class_particle* pbuffer = malloc(gpu_innerloop*sizeof(_class_particle));
  long livebatchsize = gpu_innerloop;

  #undef ABSORB0
  #undef ABSORB
  #define ABSORB0 do { DEBUG_ABSORB(); MAGNET_OFF; ABSORBED++; } while(0)
  #define ABSORB ABSORB0
  // outer loop / particle batches
  for (unsigned long long cloop=0; cloop<loops; cloop++) {
    if (loops>1) fprintf(stdout, "%d..", (int)cloop); fflush(stdout);

    // init particles
    #pragma acc parallel loop present(particles)
    for (unsigned long pidx=0 ; pidx < livebatchsize ; pidx++) {
      // generate particle state, set loop index and seed
      particles[pidx] = mcgenstate();
      _class_particle* _particle = particles + pidx;
      _particle->_uid = pidx;
      srandom(_hash((pidx+1)*(seed+1))); // _particle->state usage built into srandom macro
      particle_uservar_init(_particle);
    }

    // iterate components

    #pragma acc parallel loop present(particles)
    for (unsigned long pidx=0 ; pidx < livebatchsize ; pidx++) {
      _class_particle* _particle = &particles[pidx];
      _class_particle _particle_save;

      // origin
    if (!ABSORBED && _particle->_index == 1) {
        if (_origin_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _origin_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_origin_var._position_relative, _origin_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Progress_bar_trace(&_origin_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // source
    if (!ABSORBED && _particle->_index == 2) {
        if (_source_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _source_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_source_var._position_relative, _source_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Source_Maxwell_3_trace(&_source_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // l_monitor_before_mono
    if (!ABSORBED && _particle->_index == 3) {
        if (_l_monitor_before_mono_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _l_monitor_before_mono_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_l_monitor_before_mono_var._position_relative, _l_monitor_before_mono_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_L_monitor_trace(&_l_monitor_before_mono_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // monochromator_flat
    if (!ABSORBED && _particle->_index == 4) {
        if (_monochromator_flat_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _monochromator_flat_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_monochromator_flat_var._position_relative, _monochromator_flat_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_flat_trace(&_monochromator_flat_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // arm1
    if (!ABSORBED && _particle->_index == 5) {
        _particle->_index++;
      }

      // collimator_linear1
    if (!ABSORBED && _particle->_index == 6) {
        if (_collimator_linear1_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _collimator_linear1_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_collimator_linear1_var._position_relative, _collimator_linear1_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Collimator_linear_trace(&_collimator_linear1_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // l_monitor_before_sample
    if (!ABSORBED && _particle->_index == 7) {
        if (_l_monitor_before_sample_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _l_monitor_before_sample_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_l_monitor_before_sample_var._position_relative, _l_monitor_before_sample_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_L_monitor_trace(&_l_monitor_before_sample_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // e_monitor_before_sample
    if (!ABSORBED && _particle->_index == 8) {
        if (_e_monitor_before_sample_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_before_sample_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_e_monitor_before_sample_var._position_relative, _e_monitor_before_sample_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_E_monitor_trace(&_e_monitor_before_sample_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // PSD_monitor_before_sample
    if (!ABSORBED && _particle->_index == 9) {
        if (_PSD_monitor_before_sample_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _PSD_monitor_before_sample_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_PSD_monitor_before_sample_var._position_relative, _PSD_monitor_before_sample_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_PSD_monitor_trace(&_PSD_monitor_before_sample_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // arm2
    if (!ABSORBED && _particle->_index == 10) {
        _particle->_index++;
      }
    }

    // SPLIT with available livebatchsize 
    long mult_phonon_bvk_pg;
    livebatchsize = sort_absorb_last(particles, pbuffer, livebatchsize, gpu_innerloop, 1, &mult_phonon_bvk_pg);
    //printf("livebatchsize: %ld, split: %ld\n",  livebatchsize, mult);

    #pragma acc parallel loop present(particles)
    for (unsigned long pidx=0 ; pidx < livebatchsize ; pidx++) {
      _class_particle* _particle = &particles[pidx];
      _class_particle _particle_save;

      // phonon_bvk_pg
    if (!ABSORBED && _particle->_index == 11) {
        if (_phonon_bvk_pg_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _phonon_bvk_pg_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_phonon_bvk_pg_var._position_relative, _phonon_bvk_pg_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Phonon_BvK_PG_eigenvector_trace(&_phonon_bvk_pg_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // arm3
    if (!ABSORBED && _particle->_index == 12) {
        _particle->_index++;
      }

      // collimator_linear2
    if (!ABSORBED && _particle->_index == 13) {
        if (_collimator_linear2_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _collimator_linear2_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_collimator_linear2_var._position_relative, _collimator_linear2_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Collimator_linear_trace(&_collimator_linear2_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // psd_monitoraftersample
    if (!ABSORBED && _particle->_index == 14) {
        if (_psd_monitoraftersample_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _psd_monitoraftersample_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_psd_monitoraftersample_var._position_relative, _psd_monitoraftersample_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_PSD_monitor_trace(&_psd_monitoraftersample_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // slit1
    if (!ABSORBED && _particle->_index == 15) {
        if (_slit1_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _slit1_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_slit1_var._position_relative, _slit1_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Slit_trace(&_slit1_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // e_monitor_beforeana
    if (!ABSORBED && _particle->_index == 16) {
        if (_e_monitor_beforeana_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _e_monitor_beforeana_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_e_monitor_beforeana_var._position_relative, _e_monitor_beforeana_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_E_monitor_trace(&_e_monitor_beforeana_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // analyzer
    if (!ABSORBED && _particle->_index == 17) {
        if (_analyzer_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _analyzer_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_analyzer_var._position_relative, _analyzer_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Monochromator_flat_trace(&_analyzer_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // arm4
    if (!ABSORBED && _particle->_index == 18) {
        _particle->_index++;
      }

      // slit2
    if (!ABSORBED && _particle->_index == 19) {
        if (_slit2_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _slit2_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_slit2_var._position_relative, _slit2_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_Slit_trace(&_slit2_var, _particle); /* contains EXTEND code */
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // Emon_after_analyzer
    if (!ABSORBED && _particle->_index == 20) {
        if (_Emon_after_analyzer_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _Emon_after_analyzer_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_Emon_after_analyzer_var._position_relative, _Emon_after_analyzer_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_E_monitor_trace(&_Emon_after_analyzer_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

      // psd_detector
    if (!ABSORBED && _particle->_index == 21) {
        if (_psd_detector_var._rotation_is_identity)
          coords_get(coords_add(coords_set(x,y,z), _psd_detector_var._position_relative),&x, &y, &z);
        else
          mccoordschange(_psd_detector_var._position_relative, _psd_detector_var._rotation_relative, _particle);
        _particle_save = *_particle;
        class_PSD_monitor_trace(&_psd_detector_var, _particle);
        if (_particle->_restore)
        particle_restore(_particle, &_particle_save);
        _particle->_index++;
      }

    }

    // jump to next viable seed
    seed = seed + gpu_innerloop;
  } // outer loop / particle batches

  free(particles);
  free(pbuffer);

  printf("\n");
} /* raytrace_all_funnel */
#endif // FUNNEL

#undef x
#undef y
#undef z
#undef vx
#undef vy
#undef vz
#undef t
#undef sx
#undef sy
#undef sz
#undef p
#undef mcgravitation
#undef mcMagnet
#undef allow_backprop
#ifdef OPENACC
#ifndef MULTICORE
#undef strlen
#undef strcmp
#undef exit
#undef printf
#undef sprintf
#undef fprintf
#endif
#endif
#undef SCATTERED
#undef RESTORE
#undef RESTORE_NEUTRON
#undef STORE_NEUTRON
#undef ABSORBED
#undef ABSORB
#undef ABSORB0
/* *****************************************************************************
* instrument 'template_body_simple' and components SAVE
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_save(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_save] component origin=Progress_bar() SAVE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:120]");

  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", instrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, instrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &(instrument->counter_N[1]),&(instrument->counter_P[1]),&(instrument->counter_P2[1]),
        filename);

  }
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_save */

_class_L_monitor *class_L_monitor_save(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  SIG_MESSAGE("[_l_monitor_before_mono_save] component l_monitor_before_mono=L_monitor() SAVE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:118]");

if (!nowritefile) {
  DETECTOR_OUT_1D(
    "Wavelength monitor",
    "Wavelength [AA]",
    "Intensity",
    "L", Lmin, Lmax, nL,
    &L_N[0],&L_p[0],&L_p2[0],
    filename);
}
  #undef nL
  #undef filename
  #undef nowritefile
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_save */

_class_E_monitor *class_E_monitor_save(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_before_sample_save] component e_monitor_before_sample=E_monitor() SAVE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:131]");

if (!nowritefile) {
  DETECTOR_OUT_1D(
      "Energy monitor",
      "Energy [meV]",
      "Intensity",
      "E", Emin, Emax, nE,
      &E_N[0],&E_p[0],&E_p2[0],
      filename);
  if (S_p) printf("<E> : %g meV , E-width : %g meV \n",
   S_pE/S_p,sqrt(S_pE2/S_p - S_pE*S_pE/(S_p*S_p)) );
}
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_save */

_class_PSD_monitor *class_PSD_monitor_save(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_PSD_monitor_before_sample_save] component PSD_monitor_before_sample=PSD_monitor() SAVE [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:106]");

    if (!nowritefile) {
      DETECTOR_OUT_2D(
          "PSD monitor",
          "X position [cm]",
          "Y position [cm]",
          xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
          nx, ny,
          &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
          filename);
    }
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_save */



int save(FILE *handle) { /* called by mccode_main for template_body_simple:SAVE */
  if (!handle) siminfo_init(NULL);

  /* call iteratively all components SAVE */
  class_Progress_bar_save(&_origin_var);


  class_L_monitor_save(&_l_monitor_before_mono_var);




  class_L_monitor_save(&_l_monitor_before_sample_var);

  class_E_monitor_save(&_e_monitor_before_sample_var);

  class_PSD_monitor_save(&_PSD_monitor_before_sample_var);





  class_PSD_monitor_save(&_psd_monitoraftersample_var);


  class_E_monitor_save(&_e_monitor_beforeana_var);




  class_E_monitor_save(&_Emon_after_analyzer_var);

  class_PSD_monitor_save(&_psd_detector_var);

  if (!handle) siminfo_close(); 

  return(0);
} /* save */

/* *****************************************************************************
* instrument 'template_body_simple' and components FINALLY
***************************************************************************** */

_class_Progress_bar *class_Progress_bar_finally(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_finally] component origin=Progress_bar() FINALLY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:138]");

  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", instrument_name, dirname ? dirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3600.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_finally */

_class_L_monitor *class_L_monitor_finally(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  SIG_MESSAGE("[_l_monitor_before_mono_finally] component l_monitor_before_mono=L_monitor() FINALLY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:131]");

  destroy_darr1d(L_N);
  destroy_darr1d(L_p);
  destroy_darr1d(L_p2);
  #undef nL
  #undef filename
  #undef nowritefile
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_finally */

_class_E_monitor *class_E_monitor_finally(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_before_sample_finally] component e_monitor_before_sample=E_monitor() FINALLY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:146]");

  destroy_darr1d(E_N);
  destroy_darr1d(E_p);
  destroy_darr1d(E_p2);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_finally */

_class_PSD_monitor *class_PSD_monitor_finally(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_PSD_monitor_before_sample_finally] component PSD_monitor_before_sample=PSD_monitor() FINALLY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:119]");

  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_finally */



int finally(void) { /* called by mccode_main for template_body_simple:FINALLY */
#pragma acc update host(_origin_var)
#pragma acc update host(_source_var)
#pragma acc update host(_l_monitor_before_mono_var)
#pragma acc update host(_monochromator_flat_var)
#pragma acc update host(_arm1_var)
#pragma acc update host(_collimator_linear1_var)
#pragma acc update host(_l_monitor_before_sample_var)
#pragma acc update host(_e_monitor_before_sample_var)
#pragma acc update host(_PSD_monitor_before_sample_var)
#pragma acc update host(_arm2_var)
#pragma acc update host(_phonon_bvk_pg_var)
#pragma acc update host(_arm3_var)
#pragma acc update host(_collimator_linear2_var)
#pragma acc update host(_psd_monitoraftersample_var)
#pragma acc update host(_slit1_var)
#pragma acc update host(_e_monitor_beforeana_var)
#pragma acc update host(_analyzer_var)
#pragma acc update host(_arm4_var)
#pragma acc update host(_slit2_var)
#pragma acc update host(_Emon_after_analyzer_var)
#pragma acc update host(_psd_detector_var)
#pragma acc update host(_instrument_var)

  siminfo_init(NULL);
  save(siminfo_file); /* save data when simulation ends */

  /* call iteratively all components FINALLY */
  class_Progress_bar_finally(&_origin_var);


  class_L_monitor_finally(&_l_monitor_before_mono_var);




  class_L_monitor_finally(&_l_monitor_before_sample_var);

  class_E_monitor_finally(&_e_monitor_before_sample_var);

  class_PSD_monitor_finally(&_PSD_monitor_before_sample_var);





  class_PSD_monitor_finally(&_psd_monitoraftersample_var);


  class_E_monitor_finally(&_e_monitor_beforeana_var);




  class_E_monitor_finally(&_Emon_after_analyzer_var);

  class_PSD_monitor_finally(&_psd_detector_var);

  siminfo_close(); 

  return(0);
} /* finally */

/* *****************************************************************************
* instrument 'template_body_simple' and components DISPLAY
***************************************************************************** */

  #define magnify     mcdis_magnify
  #define line        mcdis_line
  #define dashed_line mcdis_dashed_line
  #define multiline   mcdis_multiline
  #define rectangle   mcdis_rectangle
  #define box         mcdis_box
  #define circle      mcdis_circle
  #define cylinder    mcdis_cylinder
  #define sphere      mcdis_sphere
_class_Progress_bar *class_Progress_bar_display(_class_Progress_bar *_comp
) {
  #define profile (_comp->_parameters.profile)
  #define percent (_comp->_parameters.percent)
  #define flag_save (_comp->_parameters.flag_save)
  #define minutes (_comp->_parameters.minutes)
  #define IntermediateCnts (_comp->_parameters.IntermediateCnts)
  #define StartTime (_comp->_parameters.StartTime)
  #define EndTime (_comp->_parameters.EndTime)
  #define CurrentTime (_comp->_parameters.CurrentTime)
  SIG_MESSAGE("[_origin_display] component origin=Progress_bar() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//misc/Progress_bar.comp:152]");

  printf("MCDISPLAY: component %s\n", _comp->_name);

  #undef profile
  #undef percent
  #undef flag_save
  #undef minutes
  #undef IntermediateCnts
  #undef StartTime
  #undef EndTime
  #undef CurrentTime
  return(_comp);
} /* class_Progress_bar_display */

_class_Source_Maxwell_3 *class_Source_Maxwell_3_display(_class_Source_Maxwell_3 *_comp
) {
  #define size (_comp->_parameters.size)
  #define yheight (_comp->_parameters.yheight)
  #define xwidth (_comp->_parameters.xwidth)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define dist (_comp->_parameters.dist)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define T1 (_comp->_parameters.T1)
  #define T2 (_comp->_parameters.T2)
  #define T3 (_comp->_parameters.T3)
  #define I1 (_comp->_parameters.I1)
  #define I2 (_comp->_parameters.I2)
  #define I3 (_comp->_parameters.I3)
  #define target_index (_comp->_parameters.target_index)
  #define lambda0 (_comp->_parameters.lambda0)
  #define dlambda (_comp->_parameters.dlambda)
  #define l_range (_comp->_parameters.l_range)
  #define w_mult (_comp->_parameters.w_mult)
  #define w_source (_comp->_parameters.w_source)
  #define h_source (_comp->_parameters.h_source)
  SIG_MESSAGE("[_source_display] component source=Source_Maxwell_3() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//sources/Source_Maxwell_3.comp:158]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  multiline(5, -(double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0,
                (double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0,
                (double)focus_xw/2.0,  (double)focus_yh/2.0, 0.0,
               -(double)focus_xw/2.0,  (double)focus_yh/2.0, 0.0,
               -(double)focus_xw/2.0, -(double)focus_yh/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
  #undef size
  #undef yheight
  #undef xwidth
  #undef Lmin
  #undef Lmax
  #undef dist
  #undef focus_xw
  #undef focus_yh
  #undef T1
  #undef T2
  #undef T3
  #undef I1
  #undef I2
  #undef I3
  #undef target_index
  #undef lambda0
  #undef dlambda
  #undef l_range
  #undef w_mult
  #undef w_source
  #undef h_source
  return(_comp);
} /* class_Source_Maxwell_3_display */

_class_L_monitor *class_L_monitor_display(_class_L_monitor *_comp
) {
  #define nL (_comp->_parameters.nL)
  #define filename (_comp->_parameters.filename)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Lmin (_comp->_parameters.Lmin)
  #define Lmax (_comp->_parameters.Lmax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define L_N (_comp->_parameters.L_N)
  #define L_p (_comp->_parameters.L_p)
  #define L_p2 (_comp->_parameters.L_p2)
  SIG_MESSAGE("[_l_monitor_before_mono_display] component l_monitor_before_mono=L_monitor() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/L_monitor.comp:138]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nL
  #undef filename
  #undef nowritefile
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef Lmin
  #undef Lmax
  #undef restore_neutron
  #undef L_N
  #undef L_p
  #undef L_p2
  return(_comp);
} /* class_L_monitor_display */

_class_Monochromator_flat *class_Monochromator_flat_display(_class_Monochromator_flat *_comp
) {
  #define zmin (_comp->_parameters.zmin)
  #define zmax (_comp->_parameters.zmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define zwidth (_comp->_parameters.zwidth)
  #define yheight (_comp->_parameters.yheight)
  #define mosaich (_comp->_parameters.mosaich)
  #define mosaicv (_comp->_parameters.mosaicv)
  #define r0 (_comp->_parameters.r0)
  #define Q (_comp->_parameters.Q)
  #define DM (_comp->_parameters.DM)
  #define mos_rms_y (_comp->_parameters.mos_rms_y)
  #define mos_rms_z (_comp->_parameters.mos_rms_z)
  #define mos_rms_max (_comp->_parameters.mos_rms_max)
  #define mono_Q (_comp->_parameters.mono_Q)
  SIG_MESSAGE("[_monochromator_flat_display] component monochromator_flat=Monochromator_flat() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Monochromator_flat.comp:255]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  multiline(5, 0.0, (double)ymin, (double)zmin,
               0.0, (double)ymax, (double)zmin,
               0.0, (double)ymax, (double)zmax,
               0.0, (double)ymin, (double)zmax,
               0.0, (double)ymin, (double)zmin);
  #undef zmin
  #undef zmax
  #undef ymin
  #undef ymax
  #undef zwidth
  #undef yheight
  #undef mosaich
  #undef mosaicv
  #undef r0
  #undef Q
  #undef DM
  #undef mos_rms_y
  #undef mos_rms_z
  #undef mos_rms_max
  #undef mono_Q
  return(_comp);
} /* class_Monochromator_flat_display */

_class_Arm *class_Arm_display(_class_Arm *_comp
) {
  SIG_MESSAGE("[_arm1_display] component arm1=Arm() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Arm.comp:40]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
  return(_comp);
} /* class_Arm_display */

_class_Collimator_linear *class_Collimator_linear_display(_class_Collimator_linear *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define length (_comp->_parameters.length)
  #define divergence (_comp->_parameters.divergence)
  #define transmission (_comp->_parameters.transmission)
  #define divergenceV (_comp->_parameters.divergenceV)
  #define slope (_comp->_parameters.slope)
  #define slopeV (_comp->_parameters.slopeV)
  SIG_MESSAGE("[_collimator_linear1_display] component collimator_linear1=Collimator_linear() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Collimator_linear.comp:101]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  double x;
  int i;

  
  for(x = xmin, i = 0; i <= 3; i++, x += (xmax - xmin)/3.0)
    multiline(5, x, (double)ymin, 0.0, x, (double)ymax, 0.0,
              x, (double)ymax, (double)length, x, (double)ymin, (double)length,
              x, (double)ymin, 0.0);
  line(xmin, ymin, 0,   xmax, ymin, 0);
  line(xmin, ymax, 0,   xmax, ymax, 0);
  line(xmin, ymin, length, xmax, ymin, length);
  line(xmin, ymax, length, xmax, ymax, length);
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef length
  #undef divergence
  #undef transmission
  #undef divergenceV
  #undef slope
  #undef slopeV
  return(_comp);
} /* class_Collimator_linear_display */

_class_E_monitor *class_E_monitor_display(_class_E_monitor *_comp
) {
  #define nE (_comp->_parameters.nE)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define Emin (_comp->_parameters.Emin)
  #define Emax (_comp->_parameters.Emax)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define E_N (_comp->_parameters.E_N)
  #define E_p (_comp->_parameters.E_p)
  #define E_p2 (_comp->_parameters.E_p2)
  #define S_p (_comp->_parameters.S_p)
  #define S_pE (_comp->_parameters.S_pE)
  #define S_pE2 (_comp->_parameters.S_pE2)
  SIG_MESSAGE("[_e_monitor_before_sample_display] component e_monitor_before_sample=E_monitor() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/E_monitor.comp:153]");

  printf("MCDISPLAY: component %s\n", _comp->_name);

  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nE
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef nowritefile
  #undef xwidth
  #undef yheight
  #undef Emin
  #undef Emax
  #undef restore_neutron
  #undef E_N
  #undef E_p
  #undef E_p2
  #undef S_p
  #undef S_pE
  #undef S_pE2
  return(_comp);
} /* class_E_monitor_display */

_class_PSD_monitor *class_PSD_monitor_display(_class_PSD_monitor *_comp
) {
  #define nx (_comp->_parameters.nx)
  #define ny (_comp->_parameters.ny)
  #define filename (_comp->_parameters.filename)
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define restore_neutron (_comp->_parameters.restore_neutron)
  #define nowritefile (_comp->_parameters.nowritefile)
  #define PSD_N (_comp->_parameters.PSD_N)
  #define PSD_p (_comp->_parameters.PSD_p)
  #define PSD_p2 (_comp->_parameters.PSD_p2)
  SIG_MESSAGE("[_PSD_monitor_before_sample_display] component PSD_monitor_before_sample=PSD_monitor() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//monitors/PSD_monitor.comp:126]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  multiline(5, (double)xmin, (double)ymin, 0.0,
               (double)xmax, (double)ymin, 0.0,
               (double)xmax, (double)ymax, 0.0,
               (double)xmin, (double)ymax, 0.0,
               (double)xmin, (double)ymin, 0.0);
  #undef nx
  #undef ny
  #undef filename
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef xwidth
  #undef yheight
  #undef restore_neutron
  #undef nowritefile
  #undef PSD_N
  #undef PSD_p
  #undef PSD_p2
  return(_comp);
} /* class_PSD_monitor_display */

_class_Phonon_BvK_PG_eigenvector *class_Phonon_BvK_PG_eigenvector_display(_class_Phonon_BvK_PG_eigenvector *_comp
) {
  #define hh (_comp->_parameters.hh)
  #define kk (_comp->_parameters.kk)
  #define ll (_comp->_parameters.ll)
  #define radius (_comp->_parameters.radius)
  #define yheight (_comp->_parameters.yheight)
  #define sigma_abs (_comp->_parameters.sigma_abs)
  #define sigma_inc (_comp->_parameters.sigma_inc)
  #define DW (_comp->_parameters.DW)
  #define T (_comp->_parameters.T)
  #define focus_r (_comp->_parameters.focus_r)
  #define focus_xw (_comp->_parameters.focus_xw)
  #define focus_yh (_comp->_parameters.focus_yh)
  #define focus_aw (_comp->_parameters.focus_aw)
  #define focus_ah (_comp->_parameters.focus_ah)
  #define target_x (_comp->_parameters.target_x)
  #define target_y (_comp->_parameters.target_y)
  #define target_z (_comp->_parameters.target_z)
  #define target_index (_comp->_parameters.target_index)
  #define mode_input (_comp->_parameters.mode_input)
  #define e_steps_low (_comp->_parameters.e_steps_low)
  #define e_steps_high (_comp->_parameters.e_steps_high)
  #define verbose_input (_comp->_parameters.verbose_input)
  #define dispersion (_comp->_parameters.dispersion)
  #define V_rho (_comp->_parameters.V_rho)
  #define V_my_s (_comp->_parameters.V_my_s)
  #define V_my_a_v (_comp->_parameters.V_my_a_v)
  #define Matrix (_comp->_parameters.Matrix)
  #define q (_comp->_parameters.q)
  #define qx (_comp->_parameters.qx)
  #define qy (_comp->_parameters.qy)
  #define qz (_comp->_parameters.qz)
  #define q_x (_comp->_parameters.q_x)
  #define q_y (_comp->_parameters.q_y)
  #define q_z (_comp->_parameters.q_z)
  #define eigenvectormode (_comp->_parameters.eigenvectormode)
  #define p_call (_comp->_parameters.p_call)
  SIG_MESSAGE("[_phonon_bvk_pg_display] component phonon_bvk_pg=Phonon_BvK_PG_eigenvector() DISPLAY [Phonon_BvK_PG_eigenvector.comp:1182]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  circle("xz", 0,  yheight/2.0, 0, radius);
  circle("xz", 0, -yheight/2.0, 0, radius);
  line(-radius, -yheight/2.0, 0, -radius, +yheight/2.0, 0);
  line(+radius, -yheight/2.0, 0, +radius, +yheight/2.0, 0);
  line(0, -yheight/2.0, -radius, 0, +yheight/2.0, -radius);
  line(0, -yheight/2.0, +radius, 0, +yheight/2.0, +radius);
  #undef hh
  #undef kk
  #undef ll
  #undef radius
  #undef yheight
  #undef sigma_abs
  #undef sigma_inc
  #undef DW
  #undef T
  #undef focus_r
  #undef focus_xw
  #undef focus_yh
  #undef focus_aw
  #undef focus_ah
  #undef target_x
  #undef target_y
  #undef target_z
  #undef target_index
  #undef mode_input
  #undef e_steps_low
  #undef e_steps_high
  #undef verbose_input
  #undef dispersion
  #undef V_rho
  #undef V_my_s
  #undef V_my_a_v
  #undef Matrix
  #undef q
  #undef qx
  #undef qy
  #undef qz
  #undef q_x
  #undef q_y
  #undef q_z
  #undef eigenvectormode
  #undef p_call
  return(_comp);
} /* class_Phonon_BvK_PG_eigenvector_display */

_class_Slit *class_Slit_display(_class_Slit *_comp
) {
  #define xmin (_comp->_parameters.xmin)
  #define xmax (_comp->_parameters.xmax)
  #define ymin (_comp->_parameters.ymin)
  #define ymax (_comp->_parameters.ymax)
  #define radius (_comp->_parameters.radius)
  #define xwidth (_comp->_parameters.xwidth)
  #define yheight (_comp->_parameters.yheight)
  #define isradial (_comp->_parameters.isradial)
  SIG_MESSAGE("[_slit1_display] component slit1=Slit() DISPLAY [/Applications/McStas-3.2.app/Contents/Resources/mcstas/3.2//optics/Slit.comp:107]");

  printf("MCDISPLAY: component %s\n", _comp->_name);
  
  if (is_unset(radius)) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
  #undef xmin
  #undef xmax
  #undef ymin
  #undef ymax
  #undef radius
  #undef xwidth
  #undef yheight
  #undef isradial
  return(_comp);
} /* class_Slit_display */


  #undef magnify
  #undef line
  #undef dashed_line
  #undef multiline
  #undef rectangle
  #undef box
  #undef circle
  #undef cylinder
  #undef sphere

int display(void) { /* called by mccode_main for template_body_simple:DISPLAY */
  printf("MCDISPLAY: start\n");

  /* call iteratively all components DISPLAY */
  class_Progress_bar_display(&_origin_var);

  class_Source_Maxwell_3_display(&_source_var);

  class_L_monitor_display(&_l_monitor_before_mono_var);

  class_Monochromator_flat_display(&_monochromator_flat_var);

  class_Arm_display(&_arm1_var);

  class_Collimator_linear_display(&_collimator_linear1_var);

  class_L_monitor_display(&_l_monitor_before_sample_var);

  class_E_monitor_display(&_e_monitor_before_sample_var);

  class_PSD_monitor_display(&_PSD_monitor_before_sample_var);

  class_Arm_display(&_arm2_var);

  class_Phonon_BvK_PG_eigenvector_display(&_phonon_bvk_pg_var);

  class_Arm_display(&_arm3_var);

  class_Collimator_linear_display(&_collimator_linear2_var);

  class_PSD_monitor_display(&_psd_monitoraftersample_var);

  class_Slit_display(&_slit1_var);

  class_E_monitor_display(&_e_monitor_beforeana_var);

  class_Monochromator_flat_display(&_analyzer_var);

  class_Arm_display(&_arm4_var);

  class_Slit_display(&_slit2_var);

  class_E_monitor_display(&_Emon_after_analyzer_var);

  class_PSD_monitor_display(&_psd_detector_var);

  printf("MCDISPLAY: end\n");

  return(0);
} /* display */

void* _getvar_parameters(char* compname)
/* enables settings parameters based use of the GETPAR macro */
{
  #ifdef OPENACC
    #define strcmp(a,b) str_comp(a,b)
  #endif
  if (!strcmp(compname, "origin")) return (void *) &(_origin_var._parameters);
  if (!strcmp(compname, "source")) return (void *) &(_source_var._parameters);
  if (!strcmp(compname, "l_monitor_before_mono")) return (void *) &(_l_monitor_before_mono_var._parameters);
  if (!strcmp(compname, "monochromator_flat")) return (void *) &(_monochromator_flat_var._parameters);
  if (!strcmp(compname, "arm1")) return (void *) &(_arm1_var._parameters);
  if (!strcmp(compname, "collimator_linear1")) return (void *) &(_collimator_linear1_var._parameters);
  if (!strcmp(compname, "l_monitor_before_sample")) return (void *) &(_l_monitor_before_sample_var._parameters);
  if (!strcmp(compname, "e_monitor_before_sample")) return (void *) &(_e_monitor_before_sample_var._parameters);
  if (!strcmp(compname, "PSD_monitor_before_sample")) return (void *) &(_PSD_monitor_before_sample_var._parameters);
  if (!strcmp(compname, "arm2")) return (void *) &(_arm2_var._parameters);
  if (!strcmp(compname, "phonon_bvk_pg")) return (void *) &(_phonon_bvk_pg_var._parameters);
  if (!strcmp(compname, "arm3")) return (void *) &(_arm3_var._parameters);
  if (!strcmp(compname, "collimator_linear2")) return (void *) &(_collimator_linear2_var._parameters);
  if (!strcmp(compname, "psd_monitoraftersample")) return (void *) &(_psd_monitoraftersample_var._parameters);
  if (!strcmp(compname, "slit1")) return (void *) &(_slit1_var._parameters);
  if (!strcmp(compname, "e_monitor_beforeana")) return (void *) &(_e_monitor_beforeana_var._parameters);
  if (!strcmp(compname, "analyzer")) return (void *) &(_analyzer_var._parameters);
  if (!strcmp(compname, "arm4")) return (void *) &(_arm4_var._parameters);
  if (!strcmp(compname, "slit2")) return (void *) &(_slit2_var._parameters);
  if (!strcmp(compname, "Emon_after_analyzer")) return (void *) &(_Emon_after_analyzer_var._parameters);
  if (!strcmp(compname, "psd_detector")) return (void *) &(_psd_detector_var._parameters);
  return 0;
}

void* _get_particle_var(char *token, _class_particle *p)
/* enables setpars based use of GET_PARTICLE_DVAR macro and similar */
{
  return 0;
}

int _getcomp_index(char* compname)
/* Enables retrieving the component position & rotation when the index is not known.
 * Component indexing into MACROS, e.g., POS_A_COMP_INDEX, are 1-based! */
{
  if (!strcmp(compname, "origin")) return 1;
  if (!strcmp(compname, "source")) return 2;
  if (!strcmp(compname, "l_monitor_before_mono")) return 3;
  if (!strcmp(compname, "monochromator_flat")) return 4;
  if (!strcmp(compname, "arm1")) return 5;
  if (!strcmp(compname, "collimator_linear1")) return 6;
  if (!strcmp(compname, "l_monitor_before_sample")) return 7;
  if (!strcmp(compname, "e_monitor_before_sample")) return 8;
  if (!strcmp(compname, "PSD_monitor_before_sample")) return 9;
  if (!strcmp(compname, "arm2")) return 10;
  if (!strcmp(compname, "phonon_bvk_pg")) return 11;
  if (!strcmp(compname, "arm3")) return 12;
  if (!strcmp(compname, "collimator_linear2")) return 13;
  if (!strcmp(compname, "psd_monitoraftersample")) return 14;
  if (!strcmp(compname, "slit1")) return 15;
  if (!strcmp(compname, "e_monitor_beforeana")) return 16;
  if (!strcmp(compname, "analyzer")) return 17;
  if (!strcmp(compname, "arm4")) return 18;
  if (!strcmp(compname, "slit2")) return 19;
  if (!strcmp(compname, "Emon_after_analyzer")) return 20;
  if (!strcmp(compname, "psd_detector")) return 21;
  return -1;
}

/* embedding file "mccode_main.c" */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
  /*  double run_num = 0; */
  time_t  t;
  clock_t ct;

#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, instrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

  ct = clock();

  // device and host functional RNG seed
  struct timeval tm;
  gettimeofday(&tm, NULL);
  mcseed = (long) tm.tv_sec*1000000 + tm.tv_usec;
  mcstartdate = (long)tm.tv_sec;  /* set start date before parsing options and creating sim file */
  // init global _particle.randstate for random number use
  // during init(), finally() and display(). NOTE: during trace, a local
  // "_particle" variable is present and thus used instead.
  srandom(_hash(mcseed-1));

#ifdef USE_MPI
  /* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      instrument_name, instrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per noe */
  }
#endif /* USE_MPI */

#ifdef OPENACC
#ifdef USE_MPI
  int num_devices = acc_get_num_devices(acc_device_nvidia);
  if(num_devices>0){
    int my_device = mpi_node_rank % num_devices;
    acc_set_device_num( my_device, acc_device_nvidia );
    printf("Have found %d GPU devices on rank %d. Will use device %d.\n", num_devices, mpi_node_rank, my_device);
  }else{
    printf("There was an issue probing acc_get_num_devices, fallback to host\n");
    acc_set_device_type( acc_device_host );
  }
#endif
#endif

  /* *** parse options ******************************************************* */
  SIG_MESSAGE("[" __FILE__ "] main START");
  mcformat = getenv(FLAVOR_UPPER "_FORMAT") ?
             getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  instrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */


#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif


/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */


  // init executed by master/host
  siminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("[" __FILE__ "] main INITIALISE");
  init();


#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#ifdef USE_MPI
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

// MT specific init, note that per-ray init is empty
#if RNG_ALG == 2
  mt_srandom(mcseed);
#endif


// main raytrace work loop
#ifndef FUNNEL
  // legacy version
  raytrace_all(mcncount, mcseed);
#else
  MPI_MASTER(
  // "funneled" version in which propagation is more parallelizable
  printf("\nNOTE: CPU COMPONENT grammar activated:\n 1) \"FUNNEL\" raytrace algorithm enabled.\n 2) Any SPLIT's are dynamically allocated based on available buffer size. \n");
	     );
  raytrace_all_funnel(mcncount, mcseed);
#endif


#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif


  // save/finally executed by master node/thread/host
  finally();


#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */


  return 0;
} /* mccode_main */
/* End of file "mccode_main.c". */

/* end of generated C code ./phonon_eigenvector.c */
