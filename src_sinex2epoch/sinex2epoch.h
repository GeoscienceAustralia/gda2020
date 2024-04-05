#define		ON				1
#define		OFF				0
#define		MAXLINE				256
#define		MAXLINE_LONG			256
#define		MAXNAME				256
#define         MAS2RAD         		(1./(6.48e8/PI))
#define		VCV_VARF			1.	
#define		DEWEIGHT_COORDINATE_FACTOR	1.
#define		DEWEIGHT_VELOCITY_FACTOR	1.
#define		VEL_SIGMA_MY			0.0005
#define		MAX_VEL_DIFF_MY  		0.0015
#define         MAX_VEL_SIGMA_MY                0.0015
#define         NTYPEB                          0.02
#define         ETYPEB                          0.02
#define         UTYPEB                          0.03
#define		VERSION				1.2 
#define 	MIN(a,b) 			(((a)<(b))?(a):(b))
#define 	MAX(a,b) 			(((a)>(b))?(a):(b))

#include "matrix.h"

  typedef struct {
    char name[MAXNAME];		/* station name                         */
    char ptcode[MAXNAME];	/* point code                           */
    char domes[MAXNAME];	/* domes number                         */
    char solid[MAXNAME];	/* solution id                          */
    char desc[MAXNAME];	        /* station description                  */
    int deweight;		/* deweight flag			*/
    double mjd;			/* epoch                                */
    double mmjd;                /* mean epoch of data                   */
    double smjd;                /* start epoch of data                  */
    double emjd;                /* end epoch of data                    */
    double x;			/* x cartesian coordinate               */
    double y;			/* y cartesian coordinate               */
    double z;			/* z cartesian coordinate               */
    double vx;			/* x cartesian velocity                 */
    double vy;			/* y cartesian velocity                 */
    double vz;			/* z cartesian velocity                 */
    double svx;			/* x cartesian velocity sigma           */
    double svy;			/* y cartesian velocity sigma           */
    double svz;			/* z cartesian velocity sigma           */
    int typeb;
    double n_typeb;		/* type B uncertanity			*/
    double e_typeb;		/* type B uncertanity			*/
    double h_typeb;		/* type B uncertanity			*/
  } STATION;

  typedef struct {
    int nsta;			/* # of sites input                     */
    int dof;			/* degrees of freedom input             */
    int npar;			/* # of parameters input                */
    int velocity_found;		/* are VEL parameters present           */
    double varf;		/* variance factor of input             */
    double smjd;                /*                                      */
    double emjd;                /*                                      */
    char snxVersion[4]; /* version of input SINEX file          */   
    STATION *Station;		/* station information                  */
    MATRIX Q;			/* cofactor matrix of the solution      */
    char snxComments[10][MAXLINE]; /* list of comments */
  } SOLUTION;

  typedef struct {
    int verbose;
    int fast;
    int euler;
    int vel_reset_then_end;
    int reference_epoch_year;
    int reference_epoch_doy;
    char infile[MAXLINE];
    double ntypeb;
    double etypeb;
    double utypeb;
    double rx;
    double ry;
    double rz;
    double mjd;
    double varf;
    int euler_vcv;
    double q11;
    double q12;
    double q13;
    double q22;
    double q23;
    double q33;
    SOLUTION sol;
  } GLOBAL;

SOLUTION read_sinex(char *filename);
void write_sinex(SOLUTION sol, char *outfile);
void fatal_error(char *str);
void project2epoch(SOLUTION sol,int yy, int doy, int vcv_flag);
void project2meanepoch(SOLUTION sol, int vcv_flag);
void project2endepoch(SOLUTION sol, int vcv_flag);
void velocity_reset(SOLUTION sol, double maxdiffmy, double max_vel_sig_my);
int read_commandline(int argc, char **argv);
void usage(void);
void scale_vcv(SOLUTION sol, double varf);
void deweight(SOLUTION sol, double deweight_coordinate_factor, double deweight_velocity_factor);
void apply_typeb_uncert(SOLUTION sol, char *infile);
void project_coordinate_only_2epoch(SOLUTION sol,int yy, int doy);
