/*  under /usr/include  */
/* #include <getopt.h> */ 

/*  own header file (needed for some architectures) */
#include "getopt.h" 

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>

#define BAD_OPTION -100 /* taken from own local "getopt.h", which is not active now */

/* The following external Fortran90 SetHSMR routines are
 * to pass various values to the HSMR program.           */

extern void sethsmrflags_(int *,int *,int *,int *,int *,
                          int *,int *,int *,int *,int *,
                          int *,int *,int *,int *,int *,
                          int *,int *,int *,int *     );
extern void sethsmrpars_(int *,int *,int *,int *);
extern void sethsmrshiftpars_(double *,double *,double *,double *);
extern void sethsmrfiles_(char *,char *,int,int);
extern void sethsmrconvtol_(double *);
extern void sethsmrinttol_(double *);
extern void sethsmrmaxit_(int *);

/* the main F90 subroutine */
extern void hsmrcc_(void);  

/* routine for setting general print level variable - iprint */
extern void sethsmr_iprint_(int *); 

extern void qenter_c_(char *, int *);

#define SHORTARGS "hgpvsryaTqfoG:C:V:A:B:t:c:m:i:d"

#define HELP       0
#define GENINT     1
#define PRINT      2
#define VERBOSE    3
#define SILENT     4
#define RESTART    5
#define TTYRUN     6
#define AVGVAL     7
#define NOAVG      8
#define REPTIME    9
#define OUTFILE   10
#define SHFTCORE  11
#define SHFTVIRT  12
#define SHFTACTV  13
#define BEETAVAL  14
#define INTTOL    15
#define CONVTOL   16
#define MAXIT     17
#define NCMO      18
#define NHOLE     19
#define NPTCL     20
#define NEMO      21
#define NOHSMR    22
#define HEFF      23
#define TAMP      24
#define QUAD      25
#define FULL      26
#define NOAPND    27
#define RDGUESS   28

#define IPR       29   /*  parameter for general print-out  */

#define DIRAC     30   /* flag for own DIRAC input  */

#define GMSCHK    31   /* flag for checking Gamess integral files */

#define NOPTS     32   /* always update when you add new flag  */

int opts[NOPTS]={0}; /* initialize to zero */

struct option longopts[] =
{
 /* { name               has_arg             *flag        val  } */
    {"help",           no_argument,       &opts[HELP],     1   },
    {"gen-int",        no_argument,       &opts[GENINT],   1   },
    {"print-info",     no_argument,       &opts[PRINT],    1   },
    {"silent",         no_argument,       &opts[SILENT],   1   },
    {"verbose",        no_argument,       &opts[VERBOSE],  1   },
    {"restart",        no_argument,       &opts[RESTART],  1   },
    {"tty-run",        no_argument,       &opts[TTYRUN],   1   },
    {"avg-val",        no_argument,       &opts[AVGVAL],   1   },
    {"no-avg-val",     no_argument,       &opts[NOAVG],    1   },
    {"rep-time",       no_argument,       &opts[REPTIME],  1   },
    {"out-file",       required_argument, &opts[OUTFILE],  1   },
    {"core-shift",     required_argument, &opts[SHFTCORE], 1   },
    {"virt-shift",     required_argument, &opts[SHFTVIRT], 1   },
    {"actv-shift",     required_argument, &opts[SHFTACTV], 1   },
    {"beeta",          required_argument, &opts[BEETAVAL], 1   },
    {"int-tol",        required_argument, &opts[INTTOL],   1   },
    {"conv-tol",       required_argument, &opts[CONVTOL],  1   },
    {"max-it",         required_argument, &opts[MAXIT],    1   },
    {"ncmo",           required_argument, &opts[NCMO],     1   },
    {"nhole",          required_argument, &opts[NHOLE],    1   },
    {"nptcl",          required_argument, &opts[NPTCL],    1   },
    {"nemo",           required_argument, &opts[NEMO],     1   },
    {"no-hsmr",        no_argument,       &opts[NOHSMR],   1   },
    {"heff",           no_argument,       &opts[HEFF],     1   },
    {"tamp",           no_argument,       &opts[TAMP],     1   },
    {"quadratic",      no_argument,       &opts[QUAD],     1   },
    {"full",           no_argument,       &opts[FULL],     1   },
    {"no-append",      no_argument,       &opts[NOAPND],   1   },
    {"read-guess",     no_argument,       &opts[RDGUESS],  1   },
    {"ipr",            required_argument, &opts[IPR],      1   },
    {"DIRAC-input",    no_argument,       &opts[DIRAC],    1   },
    {"GMS-fi check",   no_argument,       &opts[GMSCHK],   1   },
    { 0,                   0,                 0,           0   }
};

 /* define default values for different parameters */

#define DEFAULT_SHFTCORE 0.0
#define DEFAULT_SHFTVIRT 0.0
#define DEFAULT_SHFTACTV 0.0

#define DEFAULT_BEETAVAL 0.0

#define DEFAULT_CONVTOL_TAMP  1.0e-08
#define DEFAULT_CONVTOL_HEFF  1.0e-06

#define DEFAULT_INTTOL 1.0e-12
#define DEFAULT_MAXIT  120

#define DEFAULT_NCMO  -1
#define DEFAULT_NEMO  -1
#define DEFAULT_NHOLE -1
#define DEFAULT_NPTCL -1

void usage(prog)
    char *prog;
{
    printf("\n");
    printf("NOTE : Choosing help option disables further processing of other options.\n");
    printf("*************************************************************************\n");
    printf("\n");
    printf("Program usage: %s [-hvsryTidqfagpG] [-oCVABmci args ] [--options] input-file \n",prog);
    printf("\n");
}

void help()
{
    static char  *help_msg[] = {
 "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
 "|           *** Various options for HSMRCC program ***              |",
 "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
 "|                                                                   |",
 "| General Options :                                                 |",
 "| *****************                                                 |",
 "|                                                                   |",
 "|  -h --help             Print program usage information.           |",
 "|  -v --verbose          Produce verbose output.                    |",
 "|  -s --silent           Produce non-verbose output.                |",
 "|  -r --restart          Restart an incomplete HSMRCC run.          |",
 "|  -y --tty-run          Output the results of into terminal.       |",
 "|  -T --rep-time         Report timings for diagram routines.       |",
 "|     --no-append        Do not append in case of restart run.      |",
 "|     --read-guess       Read Initial Guess for T-amplitudes.       |",
 "|  -i --ipr=[IPRINT]     Read control print-out variable.           |",
 "|  -G                    Check GAMESS integrals files               |",
 "|                                                                   |",
 "|   DIRAC read-in :                                                 |",
 "| *****************                                                 |",
 "|  -d --DIRAC-input      Read all data from DIRAC enviroment.       |",
 "|                                                                   |",
 "| HSMRCC Options :                                                  |",
 "| ****************                                                  |",
 "|                                                                   |",
 "|  -q --quadratic        Enable quadratic approximation             |",
 "|                            for renormalization terms.             |",
 "|  -f --full             All renormalization terms.                 |",
 "|  -a --avg-val          Average valence MO energy values.          |",
 "|     --no-avg-val       Do not Average valence MO energy values.   |",
 "|  -o --out-file=[FILE]  Direct HSMRCC output to FILE.              |",
 "|                                                                   |",
 "| For Energy Shifting and Extrapolation :                           |",
 "| ***************************************                           |",
 "|                                                                   |",
 "|  -C --core-shift=[SHIFT]                                          |",
 "|  -V --virt-shift=[SHIFT]                                          |",
 "|  -A --actv-shift=[SHIFT]                                          |",
 "|  -B --beeta=[BEETA]                                               |",
 "|                                                                   |",
 "|      SHIFT = Amount of shift                                      |",
 "|      BEETA = extrapolation parameter                              |",
 "|                                                                   |",
 "| For HSMRCC Parameters:                                            |",
 "| **********************                                            |",
 "|                                                                   |",
 "|     --ncmo=[NCMO]        number of core orbitals.                 |",
 "|     --nemo=[NEMO]        number of unoccupied orbitals.           |",
 "|     --nhole=[NHOLE]      number of valence holes.                 |",
 "|     --nptcl=[NPTCL]      number of valence particles.             |",
 "|                                                                   |",
 "| For Convergence :                                                 |",
 "| *****************                                                 |",
 "|                                                                   |",
 "|  -m --max-it=[MXIT]    Maximum no. of HSMRCC iterations.          |",
 "|  -c --conv-tol=[CTOL]  Use CTOL as convergence tolerance.         |",
 "|     --heff             Use HEFF for deciding convergence.         |",
 "|     --tamp             Use TAMP for deciding convergence.         |",
 "|                                                                   |",
 "| For Integral Generation :                                         |",
 "| *************************                                         |",
 "|                                                                   |",
 "|  -g --gen-int          Generate HSMRCC integral files.            |",
 "|  -t --int-tol=[ITOL]   Integral Cutoff for HSMRCC integrals.      |",
 "|  -p --print-info       Print out infomation present in            |",
 "|                            HSMRCC integral files.                 |",
 "|     --no-hsmr          Only generate HSMRCC integrals.            |",
 "|                                                                   |",
 "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
  0};

    char **p = help_msg;
    while (*p) printf("%s\n",*p++);
}

int readdble(opt,val)
    char *opt;
    double *val;

{
    if (opt == NULL) return (-1);
    if (strlen(opt) == 0) return (-1);
    if (sscanf(opt,"%lf",val) == 1) return (0);
    if (sscanf(opt,"%le",val) == 1) return (0);
    return (-1);

}

int readint(opt,val)
    char *opt;
    int *val;

{
    if (opt == NULL) return (-1);
    if (strlen(opt) == 0) return (-1);
    if (sscanf(opt,"%d",val) == 1) return (0);
    return (-1);

}

/*   Main program   */

#if defined STANDALONE_HSMRCC
int main (argc, argv) 
#else
int main_hsmrcc (argc, argv)  
#endif
   int argc;
   char **argv;

{

   int c;
   int option_index=0;
   char *inpfile;
   char *outfile;
   double scval;
   double svval;
   double saval;
   double beeta;
   double inttol, convtol;
   int maxit;
   int ncmo, nemo, nhole, nptcl;
   pid_t pid;
   int ipr; /* print level... */

   /* name for QEnter/QExit stuff  */
   char routine_name[]="hsmrmain.c:int main(argc,argv)"; 
   qenter_c_(routine_name, (int *) strlen(routine_name));

   /* short printout */
    printf("\n");
    printf("*************************************************************************\n");
    printf("                  HSMRCC starting from main C-driver.                    \n");
    printf("*************************************************************************\n");
    printf("\n");
    
   /* change argv[0] to contain only base name */

   argv[0] = basename(argv[0]);

   /* pointer initialization */

   inpfile = NULL;
   outfile = NULL;

   /* process the command line options */

   while (1) {
      c = getopt_long(argc,argv,SHORTARGS,longopts,&option_index);
      if (c == EOF || c == BAD_OPTION) break;
      switch (c) {
         case 0:
            if (longopts[option_index].has_arg == no_argument) break;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '%s' must have non-zero length\n",
                      argv[0],longopts[option_index].name); return (-1); }
            switch (option_index)
               {
               case OUTFILE:
                  outfile = strdup(optarg);
                  break;
               case SHFTCORE:
                  if ( readdble(optarg,&scval) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case SHFTVIRT:
                  if ( readdble(optarg,&svval) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case SHFTACTV:
                  if ( readdble(optarg,&saval) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case BEETAVAL:
                  if ( readdble(optarg,&beeta) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case INTTOL:
                  if ( readdble(optarg,&inttol) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case CONVTOL:
                  if ( readdble(optarg,&convtol) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case MAXIT:
                  if ( readint(optarg,&maxit) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case NCMO:
                  if ( readint(optarg,&ncmo) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case NEMO:
                  if ( readint(optarg,&nemo) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case NHOLE:
                  if ( readint(optarg,&nhole) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               case NPTCL:
                  if ( readint(optarg,&nptcl) != 0 ) {
                     printf("%s: Error reading arguments for option '%s'\n",
                            argv[0],longopts[option_index].name); return (-1); }
                  break;
               default:
                  break;
               }
            break;

         case 'h':
            opts[HELP] = 1;
            break;

         case 'g':
            opts[GENINT] = 1;
            break;

         case 'p':
            opts[PRINT] = 1;
            break;

         case 'v':
            opts[VERBOSE] = 1;
            break;

         case 's':
            opts[SILENT] = 1;
            break;

         case 'r':
            opts[RESTART] = 1;
            break;

         case 'y':
            opts[TTYRUN] = 1;
            break;

         case 'G':
            opts[GMSCHK] = 1;
            break;

         case 'd':
            opts[DIRAC] = 1;
            break;

         case 'a':
            opts[AVGVAL] = 1;
            break;

         case 'T':
            opts[REPTIME] = 1;
            break;

         case 'q':
            opts[QUAD] = 1;
            break;

         case 'f':
            opts[FULL] = 1;
            break;

         case 'o':
            opts[OUTFILE] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            outfile = strdup(optarg);
            break;

         case 'C':
            opts[SHFTCORE] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&scval) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'V':
            opts[SHFTVIRT] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&svval) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'A':
            opts[SHFTACTV] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&saval) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'B':
            opts[BEETAVAL] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&beeta) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 't':
            opts[INTTOL] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&inttol) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'c':
            opts[CONVTOL] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readdble(optarg,&convtol) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'm':
            opts[MAXIT] = 1; 
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%c' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readint(optarg,&maxit) != 0 ) {
               printf("%s: Error reading arguments for option '-%c'\n",
                      argv[0],c); return (-1); }
            break;

         case 'i':
            opts[IPR] = 1;
            if (strlen(optarg) == 0) {
               printf("%s: argument for option '-%i' must have non-zero length\n",
                      argv[0],c); return (1); }
            if ( readint(optarg,&ipr) != 0 ) {
               printf("%s: Error reading arguments for option '-%i'\n",
                      argv[0],c); return (-1); }
            break;

         case '?':
            break;
  
         default:
            break;
      }
   }

   if (c == BAD_OPTION) return (-1);

   if (opts[HELP]) { usage(argv[0]); help(); return (0);}

   /*  Don't care about gamess-like input file if this is DIRAC based run  */
   if (optind == argc && !opts[DIRAC] ) {
      printf("\n ... error in input:  optind, argc %d  %d \n", optind, argc);
      printf("%s: Error in hsmrmain.c: input file not specified\n",argv[0]); 
   return (-1); }

   if (optind < argc) {
      if (strlen(argv[optind]) == 0) {
         printf("%s: Error: input file not specified\n",argv[0]); return (-1); }
      inpfile = strdup(argv[optind]); optind++; }

   if (optind < argc) {
      printf("%s: Error: incorrect number of non-option arguments\n",argv[0]);
      return (-1); }

   /* all the options have been read
    * process the option arguments,
    * identifying Errors/inconsitencies,
    * validity of the data */

   /* if -r/--restart can not be used with *
    * -a/--avg-val, --full, --quad, --heff *
    * --tamp , --no-avg-val, --core-shift  *
    * or --virt-shift */

   if (opts[RESTART] && opts[AVGVAL]) {
      printf("%s: Error: Option '-a/--avg-val' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[NOAVG]) {
      printf("%s: Error: Option '--no-avg-val' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[HEFF]) {
      printf("%s: Error: Option '--heff' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[TAMP]) {
      printf("%s: Error: Option '--tamp' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[QUAD]) {
      printf("%s: Error: Option '--quad' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[FULL]) {
      printf("%s: Error: Option '--full' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[SHFTCORE]) {
      printf("%s: Error: Option '--core-shift' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}
   if (opts[RESTART] && opts[SHFTVIRT]) {
      printf("%s: Error: Option '--core-shift' can not be used with '-r/--restart' option.\n",argv[0]);
      return (-1);}

   /* The following options are incompatible *
    *     --tty-run    &  --out-file         *
    *     --quadratic  &  --full             *
    *     --heff       &  --tamp             *
    *     --restart    &  --gen-int          *
    *     --restart    &  --read-guess       *
    *     --verbose    &  --silent           *
    *     --silent     &  --rep-time         *
    *     --avg-val    &  --no-avg-val       */

   if (opts[TTYRUN] && opts[OUTFILE]) {
      printf("%s: Error: Incompatible option pairs => (-y,-o)/(--tty-run,--out-file)\n",argv[0]);
      return (-1); }
   if (opts[QUAD] && opts[FULL]) {
      printf("%s: Error: Incompatible option pairs => (-q,-f)/(--quadratic,--full)\n",argv[0]);
      return (-1); }
   if (opts[HEFF] && opts[TAMP]) {
      printf("%s: Error: Incompatible option pairs => (--heff,--tamp)\n",argv[0]);
      return (-1); }
   if (opts[RESTART] && opts[GENINT]) {
      printf("%s: Error: Incompatible option pairs => (-g,-r)/(--restart,--gen-int)\n",argv[0]);
      return (-1); }
   if (opts[RESTART] && opts[RDGUESS]) {
      printf("%s: Error: Incompatible option pairs => (-r/--restart,--read-guess)\n",argv[0]);
      return (-1); }
   if (opts[VERBOSE] && opts[SILENT]) {
      printf("%s: Error: Incompatible option pairs => (-v,-s)/(--verbose,--silent)\n",argv[0]);
      return (-1); }
   if (opts[SILENT] && opts[REPTIME]) {
      printf("%s: Error: Incompatible option pairs => (-T,-s)/(--rep-time,--silent)\n",argv[0]);
      return (-1); }
   if (opts[AVGVAL] && opts[NOAVG]) {
      printf("%s: Error: Incompatible option pairs => (-a/--avg-val,--no-avg-val)\n",argv[0]);
      return (-1); }

   /* the following options are meaningful *
    * only with -g/--gen-int option        *
    *
    *   --ncmo
    *   --nemp
    *   --nhole
    *   --nptcl
    *   --no-hsmr
    *   -t/--int-tol
    *                                      */

   if (!opts[GENINT] && opts[NOHSMR]) {
      printf("%s: Error: Option '--no-hsmr' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }
   if (!opts[GENINT] && opts[NCMO]) {
      printf("%s: Error: Option '--ncmo' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }
   if (!opts[GENINT] && opts[NHOLE]) {
      printf("%s: Error: Option '--nhole' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }
   if (!opts[GENINT] && opts[NPTCL]) {
      printf("%s: Error: Option '--nptcl' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }
   if (!opts[GENINT] && opts[NEMO]) {
      printf("%s: Error: Option '--nemo' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }
   if (!opts[GENINT] && opts[INTTOL]) {
      printf("%s: Error: Option '-t/--int-tol' requires '-g/--gen-int' option\n",argv[0]);
      return (-1); }

   /* --no-append is meaningful only when --restart is chosen */

   if (opts[NOAPND] && !opts[RESTART]) {
      printf("%s: Error: Option '--no-append' is meaningful only with '-r/--restart' option.\n",argv[0]);
      return (-1);
   }

   /* when --no-hsmr option is chosen,  *
    * hsmrcc options are not meaningful */

   if (opts[NOHSMR]) {
      if (opts[RESTART]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-g/--gen-int"); return (-1); }
      if (opts[SHFTCORE]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-C/--core-shift"); return (-1); }
      if (opts[SHFTVIRT]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-V/--virt-shift"); return (-1); }
      if (opts[SHFTACTV]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-A/--actv-shift"); return (-1); }
      if (opts[BEETAVAL]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-B/--beeta"); return (-1); }
      if (opts[CONVTOL]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-c/--conv-tol"); return (-1); }
      if (opts[HEFF]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"--heff"); return (-1); }
      if (opts[TAMP]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"--tamp"); return (-1); }
      if (opts[MAXIT]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-m/--max-it"); return (-1); }
      if (opts[AVGVAL]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-a/--avg-val"); return (-1); }
      if (opts[QUAD]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-q/--quad"); return (-1); }
      if (opts[FULL]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-f/--full"); return (-1); }
      if (opts[REPTIME]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"-T/--rep-time"); return (-1); }
      if (opts[RDGUESS]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"--read-guess"); return (-1); }
      if (opts[IPR]) {
         printf("%s: Error: Option '%s' not meaningful with '--no-hsmr' option.\n",
                argv[0],"--ipr"); return (-1); }
   }

   /* Check whether integer input makes sense */

   if (opts[MAXIT] && (maxit<1)) {
      printf("%s: Error: Invalid value for MAXIT.\n",argv[0]); return (-1); }
   if (opts[NCMO] && (ncmo<1)) {
      printf("%s: Error: Invalid value for NCMO.\n",argv[0]); return (-1); }
   if (opts[NHOLE] && (nhole<1)) {
      printf("%s: Error: Invalid value for NHOLE.\n",argv[0]); return (-1); }
   if (opts[NPTCL] && (nptcl<1)) {
      printf("%s: Error: Invalid value for NPTCL.\n",argv[0]); return (-1); }
   if (opts[NEMO] && (nemo<1)) {
      printf("%s: Error: Invalid value for NEMO.\n",argv[0]); return (-1); }

   /* Check whether tolerance values make sense 
    * too high and negative values are not meaningful */

   if (opts[INTTOL] && (inttol<0.0e00)) {
      printf("%s: Error: Negative INTTOL \n",argv[0]); return (-1); }
   if (opts[INTTOL] && (inttol>1.0e00)) {
      printf("%s: Error: INTTOL too high \n",argv[0]); return (-1); }
   if (opts[CONVTOL] && (convtol<0.0e00)) {
      printf("%s: Error: Negative CONVTOL \n",argv[0]); return (-1); }
   if (opts[CONVTOL] && (convtol>1.0e00)) {
      printf("%s: Error: CONVTOL too high \n",argv[0]); return (-1); }
   if (opts[IPR] && (ipr < 0 || ipr > 100 )) {
      printf("%s: Error: ipr out of range <0;100>, ipr=\n",argv[0]); return (-1); }

   /* Check whether shift parameters are meaningful. */

/* if (opts[SHFTCORE]) {
      if (scval<0.0e00) {printf("%s: Error: Negative Core SHIFT value \n",argv[0]);
         return (-1);}
   }
   if (opts[SHFTVIRT]) {
      if (svval<0.0e00) {printf("%s: Error: Negative Virtual SHIFT value \n",argv[0]);
         return (-1);}
   }
*/

/*           
  ------------------------------------------------------------------------------------------
   Set IPRINT value first so that next F90 routines will have the printing parameter ready.
  ------------------------------------------------------------------------------------------
*/
   if (opts[IPR]) { 
     printf("hsmrmain.c: entering iprint parameter, ipr=%d\n",ipr);
     sethsmr_iprint_(&ipr);
   }

   /* Set Defaults for Unselected Flags */

   if (!opts[TAMP]    && !opts[HEFF])  opts[TAMP]   = 1;
   if (!opts[FULL]    && !opts[QUAD])  opts[QUAD]   = 0;
   if (!opts[AVGVAL]  && !opts[NOAVG]) opts[AVGVAL] = 0;
   if (!opts[VERBOSE] && !opts[SILENT]) {
      opts[VERBOSE] = 1;
      if (opts[TTYRUN]) opts[VERBOSE] = 0;
   }

   if (!opts[SHFTCORE]) {
      scval = DEFAULT_SHFTCORE;
   }

   if (!opts[SHFTVIRT]) {
      svval = DEFAULT_SHFTVIRT;
   }

   if (!opts[SHFTACTV]) {
      saval = DEFAULT_SHFTACTV;
   }

   if (!opts[BEETAVAL]) {
      beeta = DEFAULT_BEETAVAL;
   }
   if (!opts[CONVTOL]) {
      convtol = DEFAULT_CONVTOL_TAMP;
      if (opts[HEFF]) convtol = DEFAULT_CONVTOL_HEFF;
   }
   if (!opts[INTTOL]) inttol = DEFAULT_INTTOL;

   if (!opts[NCMO])  ncmo  = DEFAULT_NCMO ;
   if (!opts[NEMO])  nemo  = DEFAULT_NEMO ;
   if (!opts[NHOLE]) nhole = DEFAULT_NHOLE;
   if (!opts[NPTCL]) nptcl = DEFAULT_NPTCL;
   if (!opts[MAXIT]) maxit = DEFAULT_MAXIT;

   /* Pass the options to HSMR program.   *
    * sethsmr* routines store the passed  *
    * values in corresponding module vars */

   sethsmrflags_(&opts[GENINT] , &opts[PRINT]   , &opts[VERBOSE] ,
                 &opts[RESTART], &opts[TTYRUN]  , &opts[AVGVAL]  ,
                 &opts[REPTIME], &opts[QUAD]    , &opts[TAMP]    ,
                 &opts[NOHSMR] , &opts[MAXIT]   , &opts[INTTOL]  ,
                 &opts[CONVTOL], &opts[SHFTCORE], &opts[SHFTVIRT],
                 &opts[NOAPND] , &opts[RDGUESS] , &opts[DIRAC]   ,
                 &opts[GMSCHK] );

   sethsmrpars_(&ncmo,&nhole,&nptcl,&nemo);
   sethsmrshiftpars_(&scval,&svval,&saval,&beeta);
   sethsmrconvtol_(&convtol);
   sethsmrinttol_(&inttol); /* */
   sethsmrmaxit_(&maxit); /* max.number of iterations */

   if ( !opts[DIRAC] ) { /* Gamess files */
     if (outfile==NULL) {
       sethsmrfiles_(inpfile,outfile,strlen(inpfile),0);
     }
     else {
       sethsmrfiles_(inpfile,outfile,strlen(inpfile),strlen(outfile));
     }
   }
   else /* the HSMRCC run is based solely upon the DIRAC program  */
     {
       /* sethsmrfiles_(inpfile,outfile,0,0); */
       /* inpfile = "dirac_with_HSCC"; */
        inpfile = "Dirac_HSMRCC"; /* common part of file name of DIRAC prepared files */
        sethsmrfiles_(inpfile,outfile,strlen(inpfile),0);
     }

   /* Invoke HSMRCC Driver program */

   if (opts[TTYRUN]) {
    /* run the main F90 program in the terminal window */
       hsmrcc_();
       printf(" hsmrmain.c : HSMRCC fortran90 part finished ... \n");
   }
   else { /* run the program as the backgroung job */
      if ((pid = fork()) == -1) {printf("%s: Error: Unable to fork child process.\n",argv[0]); return (-1);}
      if (pid==0) { hsmrcc_(); exit(0);}
      if (pid!=0) {
         printf("%s: Submitting Background Job with PID = %d\n",argv[0],pid);
         fclose(stdin);
         fclose(stdout);
         fclose(stderr);
      }
   }

   qexit_c_(routine_name, (int *) strlen(routine_name));
   exit(0);
}


void call_dirac_sethsmrfiles_(void)
/* call from F90 file .. not working ...  */
{
  char *inpfile;
  char *outfile;

  inpfile = NULL;
  outfile = NULL;

 inpfile = "Dirac_HSMRCC"; /* common part of file name of DIRAC prepared files */
 sethsmrfiles_(inpfile,outfile,strlen(inpfile),0);
}


