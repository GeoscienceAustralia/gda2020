#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include "geodetic.h"
#include "time.h"
#include "sinex2epoch.h"

GLOBAL global;

int main(int argc, char **argv)
{
    SOLUTION sol;

    fprintf(stderr,"+ sinex2epoch version %lg\n", VERSION);
 
    read_commandline(argc, argv);

    sol = read_sinex(global.infile);

    scale_vcv(sol,global.varf);

    if (sol.velocity_found == ON) {
	    project2meanepoch(sol,!global.fast);

	    if (global.euler == ON) 
		velocity_reset(sol,MAX_VEL_DIFF_MY,MAX_VEL_SIGMA_MY); 

	    if (global.vel_reset_then_end == ON)
		exit(1);

	    project2epoch(sol,global.reference_epoch_year,global.reference_epoch_doy,!global.fast);

	    if (global.euler == ON) 
		velocity_reset(sol, 0., MAX_VEL_SIGMA_MY); 
    }
    else {
           velocity_reset(sol, 0., MAX_VEL_SIGMA_MY); 
           project_coordinate_only_2epoch(sol,global.reference_epoch_year,global.reference_epoch_doy);
    }

    apply_typeb_uncert(sol,"typeb.dat");

    write_sinex(sol,"SNXEPO.SNX");

    return 1;
}


/*
 * apply_typeb_uncert()
 * 	n, e, h, site, solid
 */
void apply_typeb_uncert(SOLUTION sol, char *filename) {
  int i,j,ii;
  char line[MAXLINE], name[MAXLINE], solid[MAXLINE];
  FILE *fp;
  double e, n, h, x, y, z, Qneh[3][3], Qxyz[3][3];

  fprintf(stderr,"+ reading %s\n", filename);

  if ((fp = fopen(filename, "r")) == NULL) {
      fprintf(stderr,"+ no TYPE B file found - default uncertainities applied: %lg %lg %lg m NEU\n", 
		global.ntypeb, global.etypeb, global.utypeb);
      for (i=0; i<sol.nsta; i++) {
		sol.Station[i].typeb = ON;
		sol.Station[i].n_typeb = global.ntypeb;
		sol.Station[i].e_typeb = global.etypeb;
		sol.Station[i].h_typeb = global.utypeb;
      }
  } else {
	  while (fgets(line, MAXLINE, fp) != NULL) {
		ii = sscanf(line, "%lf %lf %lf %s %*s %*s %s", &n, &e, &h, name, solid);
		if (ii == 5) {
			for (i=0; i<sol.nsta; i++) {
				if (strncmp(name,sol.Station[i].name,4) == 0 && strncmp(solid,sol.Station[i].solid,2) == 0) {
					sol.Station[i].typeb = ON;
					sol.Station[i].n_typeb = n;
					sol.Station[i].e_typeb = e;
					sol.Station[i].h_typeb = h;
				}
			}
		}
	  } 
  } 

  fclose(fp);

  for (i=0; i<sol.nsta; i++) {
	if (sol.Station[i].typeb == ON) {

		for (j=0; j<3; j++) {
			for (ii=0; ii<3; ii++) {
				Qneh[ii][j] = 0;
				Qxyz[ii][j] = 0;
			}
		}

		n = sol.Station[i].n_typeb;
		e = sol.Station[i].e_typeb;
		h = sol.Station[i].h_typeb;

		fprintf(stderr, "+ add type B uncertainty (NEH): %.4lf %.4lf %.4lf --> %s %s\n", n, e, h, sol.Station[i].name, sol.Station[i].solid);

		Qneh[0][0] = n * n;
		Qneh[1][1] = e * e;
		Qneh[2][2] = h * h;

		x = sol.Station[i].x;
		y = sol.Station[i].y;
		z = sol.Station[i].z;

		sigma_neh2uvw(x,y,z,Qxyz,Qneh);

		if (sol.velocity_found == ON)
			ii = i*6;
		else
			ii = i*3;

		sol.Q[ii+0][ii+0] += Qxyz[0][0];
		sol.Q[ii+0][ii+1] += Qxyz[0][1];
		sol.Q[ii+0][ii+2] += Qxyz[0][2];
		sol.Q[ii+1][ii+0] += Qxyz[1][0];
		sol.Q[ii+1][ii+1] += Qxyz[1][1];
		sol.Q[ii+1][ii+2] += Qxyz[1][2];
		sol.Q[ii+2][ii+0] += Qxyz[2][0];
		sol.Q[ii+2][ii+1] += Qxyz[2][1];
		sol.Q[ii+2][ii+2] += Qxyz[2][2];
	}
  }
}

/*
 * deweight()
 */
void deweight(SOLUTION sol, double deweight_coordinate_factor, double deweight_velocity_factor) {
int i;

	fprintf(stderr,"+ deweight flagged stations, coordinates x %lg, velocities x %lg\n", DEWEIGHT_COORDINATE_FACTOR, DEWEIGHT_VELOCITY_FACTOR);

	for (i=0; i<sol.nsta; i++) {
		
		if (sol.Station[i].deweight == ON) {

			sol.Q[i*6+0][i*6+0] *= deweight_coordinate_factor;
			sol.Q[i*6+1][i*6+1] *= deweight_coordinate_factor;
			sol.Q[i*6+2][i*6+2] *= deweight_coordinate_factor;
			sol.Q[i*6+3][i*6+3] *= deweight_velocity_factor;
			sol.Q[i*6+4][i*6+4] *= deweight_velocity_factor;
			sol.Q[i*6+5][i*6+5] *= deweight_velocity_factor;
		}	
	}
}

/*
 * scale_vcv()
 */
void scale_vcv(SOLUTION sol, double varf) {
int i,j;

	fprintf(stderr,"+ scale vcv by %lg (sqrt(vcv) x %lg)\n", varf, sqrt(varf));

        for (i=0; i<sol.nsta; i++) {
		sol.Station[i].svx *= varf;
		sol.Station[i].svy *= varf;
		sol.Station[i].svz *= varf;
	}

	for (i=0; i<MatRow(sol.Q); i++)
		for (j=0; j<MatCol(sol.Q); j++)
			sol.Q[i][j] *= varf;

}


/*
 * velocity_reset()
 * 	- reset velocity estimates of stations exceeding threshold (maxdiffmy)
 * 	- reset velocity estimates of stations which velocity uncertainy greater than threshold (max_vel_sig_my)
 * 	- reset velocity variances and covariances too for those stations to vel_sigma_my
 */
void velocity_reset(SOLUTION sol, double maxdiffmy, double max_vel_sig_my) {
int i,j,ii,jj;
double rx, ry, rz, x, y, z, vx, vy, vz, dvx, dvy, dvz, svx, svy, svz, n, e, h, svn, sve, svh, se, sn, sh, ne, nn, nh;
MATRIX A,At,Qv,Q,QAt;

	fprintf(stderr,"+ reset velocities to plate model for those stations with a horizontal velocity difference > %lg m/y\n", maxdiffmy);
	fprintf(stderr,"+ reset velocities to plate model for those stations with a horizontal velocity uncertainty > %lg m/y\n", max_vel_sig_my);
	fprintf(stderr,"+ reset velocities uncertainties to %lg m/y XYZ\n", VEL_SIGMA_MY);

	rx = global.rx * MAS2RAD;
	ry = global.ry * MAS2RAD;
	rz = global.rz * MAS2RAD;

	for (i=0; i<sol.nsta; i++) {
		x = sol.Station[i].x;
		y = sol.Station[i].y;
		z = sol.Station[i].z;

		vx =  z * ry - y * rz;
		vy = -z * rx + x * rz;
		vz =  y * rx - x * ry;

		dvx = sol.Station[i].vx - vx;
		dvy = sol.Station[i].vy - vy;
		dvz = sol.Station[i].vz - vz;

		svx = sol.Station[i].svx;
		svy = sol.Station[i].svy;
		svz = sol.Station[i].svz;

		uvw2neh(x, y, z,  dvx, dvy, dvz, &n, &e, &h);
		uvw2neh(x, y, z,  svx, svy, svz, &svn, &sve, &svh);

		sn = 0.002;
		se = 0.002;
		sh = 0.004;

		// compute the normalised residual using expected consistency with model
		// and expected uncertainty level of uncertainty
		nn = fabs(n / maxdiffmy) + fabs(svn / max_vel_sig_my);
		ne = fabs(e / maxdiffmy) + fabs(sve / max_vel_sig_my);
		nh = fabs(h / maxdiffmy) + fabs(svh / max_vel_sig_my);
		
		// limit normalisation factor to be no greater than a threshold value but greater than 1
		nn = MIN(10.,nn);
		ne = MIN(10.,ne);
		nh = MIN(10.,nh);
		nn = MAX(1.,nn);
		ne = MAX(1.,ne);
		nh = MAX(1.,nh);

		// apply normalisation 
		sn += maxdiffmy * nn;
		se += maxdiffmy * ne;
		sh += 2. * maxdiffmy * nh;

		if (global.verbose) {
			fprintf(stdout, "%s %4.1lf %4.1lf %4.1lf NEH NRESID DIFF\n", sol.Station[i].name, 
				fabs(n / maxdiffmy), fabs(e / maxdiffmy), fabs(h / maxdiffmy));
			fprintf(stdout, "%s %4.1lf %4.1lf %4.1lf NEH NRESID SIGMA\n", sol.Station[i].name, 
				fabs(svn / max_vel_sig_my), fabs(sve / max_vel_sig_my), fabs(svh / max_vel_sig_my));
			fprintf(stdout, "%s %4.1lf %4.1lf %4.1lf NEH NRESID TOTAL\n", sol.Station[i].name, nn, ne, nh);
		}

		// 	
		// if horizontal velocity different by more than maxdiffmy mm/y then replace with plate model
		// OR 
		// if horizontal velocity uncertainty greater than max_vel_sig_my mm/y then replace with plate model
		// 	
		if (fabs(n) > maxdiffmy || fabs(e) > maxdiffmy || fabs(svn) > max_vel_sig_my || fabs(sve) > max_vel_sig_my) {
		
			if (fabs(n) > maxdiffmy || fabs(e) > maxdiffmy) {

			fprintf(stdout,
			"# %6.3lf %6.3lf %6.3lf %s %s %s %s %6.3lf %6.3lf %6.3lf m/y VEL DIFF NEH %8.4lf %8.4lf %8.4lf VEL SIG NEH** MFAIL\n",
			sn, se, sh, sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].domes, sol.Station[i].solid, n, e, h, svn, sve, svh);

			}

			else if (fabs(svn) > max_vel_sig_my || fabs(sve) > max_vel_sig_my) {

			fprintf(stdout,
			"# %6.3lf %6.3lf %6.3lf %s %s %s %s %6.3lf %6.3lf %6.3lf m/y VEL DIFF NEH %8.4lf %8.4lf %8.4lf VEL SIG NEH** SFAIL\n",
			sn, se, sh, sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].domes, sol.Station[i].solid, n, e, h, svn, sve, svh);

			}

        		A = mat_creat(3, 3, ZERO_MATRIX);
        		Q = mat_creat(3, 3, ZERO_MATRIX);

			A[0][0] = 0.; A[0][1] =  z;  A[0][2] = y;
			A[1][0] = -z; A[1][1] =  0.; A[1][2] = x;
			A[2][0] =  y; A[2][1] = -x;  A[2][2] = 0.;

			Q[0][0] = global.q11; Q[0][1] = global.q12; Q[0][2] = global.q13;
			Q[1][0] = global.q12; Q[1][1] = global.q22; Q[1][2] = global.q23;
			Q[2][0] = global.q13; Q[2][1] = global.q23; Q[2][2] = global.q33;

			At = mat_tran(A);
			QAt = mat_mul(Q,At);
			Qv = mat_mul(A,QAt);

			sol.Station[i].vx = vx;
			sol.Station[i].vy = vy;
			sol.Station[i].vz = vz;

			if (global.euler_vcv == ON) {
				sol.Station[i].svx = sqrt(Qv[0][0]);
				sol.Station[i].svy = sqrt(Qv[1][1]);
				sol.Station[i].svz = sqrt(Qv[2][2]);
			} else {
				sol.Station[i].svx = VEL_SIGMA_MY;
				sol.Station[i].svy = VEL_SIGMA_MY;
				sol.Station[i].svz = VEL_SIGMA_MY;
			}
			sol.Station[i].deweight = ON;

			mat_free(A);
			mat_free(At);
			mat_free(QAt);
			mat_free(Q);
			mat_free(Qv);

		}

		// the best stations - nothing has been flagged
		else {
			fprintf(stdout,
			"# %6.3lf %6.3lf %6.3lf %s %s %s %s %6.3lf %6.3lf %6.3lf m/y VEL DIFF NEH %8.4lf %8.4lf %8.4lf VEL SIG NEH\n",
			sn, se, sh, sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].domes, sol.Station[i].solid, n, e, h, svn, sve, svh);
		}
	}
        fprintf(stdout,"+ \n");

	if (global.fast == OFF && sol.velocity_found == ON) {
		// reset the variances and covariances of the reset velocity parameters
		for (i=0; i<sol.nsta; i++) {
			if (sol.Station[i].deweight == ON) {
				j = i*6 + 3;
				for (j = i*6+3; j<i*6+6; j++) {
					for (ii=0; ii<MatRow(sol.Q); ii++) {
						for (jj=0; jj<MatCol(sol.Q); jj++) {
							if (j == ii && j == jj) {
								sol.Q[ii][jj] = VEL_SIGMA_MY * VEL_SIGMA_MY;
							}
							else if (j == ii || j == jj) {
								sol.Q[ii][jj] = 0.;
								sol.Q[jj][ii] = 0.;
							}
						}
					}
				}
			}
		}

		// use variances and covariances from model if available (but not non-diagonal terms)
		if (global.euler_vcv == ON) {
			for (i=0; i<sol.nsta; i++) {
				sol.Q[i*6+3][i*6+3] = sol.Station[i].svx * sol.Station[i].svx;
				sol.Q[i*6+4][i*6+4] = sol.Station[i].svy * sol.Station[i].svy;
				sol.Q[i*6+5][i*6+5] = sol.Station[i].svz * sol.Station[i].svz;
			}
		}
	}
}


/*
 * project2epoch();
 *    - project coordinates and VCV to new reference epoch
 */
void project2epoch(SOLUTION sol,int yy, int doy, int vcv_flag)
{
	int i,j;
	double mjd, dt;
        MATRIX A, At, Q, QAt;

	fprintf(stderr,"+ project coordinates to the epoch of %d:%3.3d\n",yy,doy);
	if (vcv_flag == OFF) { fprintf(stderr,"+ WARNING: VCV not transformed\n"); }

	mjd = ymd2mjd(yy2yyyy(yy),1,0) + doy;
        A = mat_creat(sol.nsta*6, sol.nsta*6, ZERO_MATRIX);

	for (i=0; i<sol.nsta; i++) {
		dt = (mjd - sol.Station[i].mjd) / 365.25;
		sol.Station[i].x += dt * sol.Station[i].vx;
		sol.Station[i].y += dt * sol.Station[i].vy;
		sol.Station[i].z += dt * sol.Station[i].vz;
		A[i*6+0][i*6+0] = 1.;
		A[i*6+1][i*6+1] = 1.;
		A[i*6+2][i*6+2] = 1.;
		A[i*6+0][i*6+3] = dt;
		A[i*6+1][i*6+4] = dt;
		A[i*6+2][i*6+5] = dt;
		A[i*6+3][i*6+3] = 1.;
		A[i*6+4][i*6+4] = 1.;
		A[i*6+5][i*6+5] = 1.;
		sol.Station[i].mjd = mjd;
	}

	if (vcv_flag == ON) {

		fprintf(stderr, "+ making At(%d x %d)\n",MatCol(A),MatRow(A));
      	  	At = mat_tran(A);
		fprintf(stderr, "+ making QAt\n");
		QAt = mat_mul(sol.Q,At);
		fprintf(stderr,"+ making AQAt\n");
		Q = mat_mul(A,QAt);
		fprintf(stderr,"+ updating Q\n");
		for (i=0; i<MatRow(Q); i++) {
			for (j=0; j<MatCol(Q); j++) {
				sol.Q[i][j] = Q[i][j];
			}
		}	

		mat_free(Q);
		mat_free(At);
		mat_free(QAt);
	}

	mat_free(A);
}

/*
 * project_coordinate_only_2epoch();
 *    - project coordinates to new reference epoch
 */
void project_coordinate_only_2epoch(SOLUTION sol,int yy, int doy)
{
	int i;
	double mjd, dt;

	fprintf(stderr,"+ project coordinates to the epoch of %d:%3.3d\n",yy,doy);

	mjd = ymd2mjd(yy2yyyy(yy),1,0) + doy;

	for (i=0; i<sol.nsta; i++) {
		dt = (mjd - sol.Station[i].mjd) / 365.25;
		sol.Station[i].x += dt * sol.Station[i].vx;
		sol.Station[i].y += dt * sol.Station[i].vy;
		sol.Station[i].z += dt * sol.Station[i].vz;
		sol.Station[i].mjd = mjd;
	}
}

/*
 * project2meanepoch();
 *    - project coordinates and VCV to the individual mean epochs of each station
 */
void project2meanepoch(SOLUTION sol, int vcv_flag)
{
	int i,j;
	double dt;
        MATRIX A, At, Q, QAt;

	fprintf(stderr,"+ project station coordinates to their individual mean epoch\n");
	if (vcv_flag == OFF) { fprintf(stderr,"+ WARNING: VCV not transformed\n"); }

        A = mat_creat(sol.nsta*6, sol.nsta*6, ZERO_MATRIX);

	for (i=0; i<sol.nsta; i++) {
		dt = (sol.Station[i].mmjd - sol.Station[i].mjd) / 365.25;
		sol.Station[i].x += dt * sol.Station[i].vx;
		sol.Station[i].y += dt * sol.Station[i].vy;
		sol.Station[i].z += dt * sol.Station[i].vz;
		A[i*6+0][i*6+0] = 1.;
		A[i*6+1][i*6+1] = 1.;
		A[i*6+2][i*6+2] = 1.;
		A[i*6+0][i*6+3] = dt;
		A[i*6+1][i*6+4] = dt;
		A[i*6+2][i*6+5] = dt;
		A[i*6+3][i*6+3] = 1.;
		A[i*6+4][i*6+4] = 1.;
		A[i*6+5][i*6+5] = 1.;
		sol.Station[i].mjd = sol.Station[i].mmjd;
	}

	if (vcv_flag == ON) {

		fprintf(stderr, "+ making At(%d x %d)\n",MatCol(A),MatRow(A));
      	  	At = mat_tran(A);
		fprintf(stderr,"+ making QAt\n");
		QAt = mat_mul(sol.Q,At);
		fprintf(stderr,"+ making AQAt\n");
		Q = mat_mul(A,QAt);
		fprintf(stderr,"+ updating Q\n");
		for (i=0; i<MatRow(Q); i++) {
			for (j=0; j<MatCol(Q); j++) {
				sol.Q[i][j] = Q[i][j];
			}
		}	

		mat_free(Q);
		mat_free(At);
		mat_free(QAt);
	}

	mat_free(A);
}

/*
 * project2endepoch();
 *    - project coordinates and VCV to the individual end epochs of each station
 */
void project2endepoch(SOLUTION sol, int vcv_flag)
{
	int i,j;
	double dt;
        MATRIX A, At, Q, QAt;

	printf("+ project station coordinates to their individual end epoch\n");
	if (vcv_flag == OFF) { printf("+ WARNING: VCV not transformed\n"); }

        A = mat_creat(sol.nsta*6, sol.nsta*6, ZERO_MATRIX);

	for (i=0; i<sol.nsta; i++) {
		dt = (sol.Station[i].emjd - sol.Station[i].mjd) / 365.25;
		sol.Station[i].x += dt * sol.Station[i].vx;
		sol.Station[i].y += dt * sol.Station[i].vy;
		sol.Station[i].z += dt * sol.Station[i].vz;
		A[i*6+0][i*6+0] = 1.;
		A[i*6+1][i*6+1] = 1.;
		A[i*6+2][i*6+2] = 1.;
		A[i*6+0][i*6+3] = dt;
		A[i*6+1][i*6+4] = dt;
		A[i*6+2][i*6+5] = dt;
		A[i*6+3][i*6+3] = 1.;
		A[i*6+4][i*6+4] = 1.;
		A[i*6+5][i*6+5] = 1.;
		sol.Station[i].mjd = sol.Station[i].emjd;
	}

	if (vcv_flag == ON) {

		fprintf(stderr, "+ making At(%d x %d)\n",MatCol(A),MatRow(A));
      	  	At = mat_tran(A);
		fprintf(stderr, "+ making QAt\n");
		QAt = mat_mul(sol.Q,At);
		fprintf(stderr,"+ making AQAt\n");
		Q = mat_mul(A,QAt);
		fprintf(stderr,"+ updating Q\n");
		for (i=0; i<MatRow(Q); i++) {
			for (j=0; j<MatCol(Q); j++) {
				sol.Q[i][j] = Q[i][j];
			}
		}	

		mat_free(Q);
		mat_free(At);
		mat_free(QAt);
	}

	mat_free(A);
}

/*
 * read_sinex()
 *	- reads a SINEX file
 *	- reads the raw parameter estimates and their variance covariance information
 */
SOLUTION read_sinex(char *filename)
{
  int i, done, n, y, doy, sec, row, col;
  char line[MAXLINE_LONG], errstr[MAXLINE], type[7], name[5], ptcode[5], s1[10], domes[12], solid[2];
  double val, sd, mjd, q1, q2, q3;
  FILE *fp;
  SOLUTION sol;

  fprintf(stderr, "+ reading SINEX format variance covariance solution: %s\n", filename);

  if ((fp = fopen(filename, "r")) == NULL) 
    {
      sprintf(errstr, "UNABLE TO OPEN OUTPUT FILE %s", filename);
      fatal_error(errstr);
    }

  done = 0;
  fgets(line, MAXLINE, fp);
  sscanf(line, "%*s %*s %*s %*s %*s %d:%d:%d %*s %*s %d", &y, &doy, &sec, &sol.npar);
  sol.emjd = ymd2mjd(yy2yyyy(y),1,0) + (double) doy + (double) sec/86400.;

  sscanf(line, "%*s %*s %*s %*s %*s %*s %d:%d:%d", &y, &doy, &sec);
  sol.smjd = ymd2mjd(yy2yyyy(y),1,0) + (double) doy + (double) sec/86400.;

  sscanf(line, "%*s %s", sol.snxVersion);

  sol.velocity_found = OFF;

  // save file comments
  done = 0;
  i = 0;
  while (fgets(line, MAXLINE, fp) != NULL && !done)
    {
      if (strncmp(line, "+FILE/COMMENT", 13) == 0)
    {
      while (fgets(line, MAXLINE, fp) != NULL)
        {
          if (strncmp(line, "-FILE/COMMENT", 13) == 0)
        {
          done = 1;
          break;
        }
          sscanf(line, "%[^\n]", sol.snxComments[i]);
          i++;
        }
    }
    }

  // check for VEL parameters
  done = 0;
  while (fgets(line, MAXLINE, fp) != NULL && !done) 
    {
      if (strncmp(line, "+SOLUTION/ESTIMATE", 18) == 0) 
	{
	  while (fgets(line, MAXLINE, fp) != NULL) 
	    {
	      if (strncmp(line, "-SOLUTION/ESTIMATE", 18) == 0) 
		{
		  done = 1;
		  break;
		}

	      if (line[0] != '*') 
		{
		  sscanf(line, "%d%s", &n, type);
		  if (strstr(type, "VELX") != NULL) 
			sol.velocity_found = ON;
                }
             }
         }
     }

  if (sol.velocity_found == ON) {
  	fprintf(stderr, "+ velocity estimates found in the input SINEX\n");
  	sol.nsta = sol.npar/6; 
  } else {
  	fprintf(stderr, "+ velocity estimates NOT found in the input SINEX\n");
  	sol.nsta = sol.npar/3; 
  }

  if ((sol.Station = (STATION *) calloc(sol.nsta, sizeof(STATION))) == NULL)
    {
      fatal_error("MEMORY ALLOCATION FAILURE");
    }

  if ((sol.Q = mat_creat(sol.npar, sol.npar, UNIT_MATRIX)) == NULL)
    {
      fatal_error("MEMORY ALLOCATION FAILURE");
    }

  rewind(fp);

  i = 0;
  done = 0;
  while (fgets(line, MAXLINE, fp) != NULL && !done) 
    {
      if (strncmp(line, "+SOLUTION/ESTIMATE", 18) == 0) 
	{
	  while (fgets(line, MAXLINE, fp) != NULL) 
	    {
	      if (strncmp(line, "-SOLUTION/ESTIMATE", 18) == 0) 
		{
		  done = 1;
		  break;
		}

	      if (line[0] != '*') 
		{
		  sscanf(line, "%d%s%s%s%s%d:%d:%d", &n, type, name, ptcode, s1, &y, &doy, &sec);
		  sscanf(line + 47, "%le%le", &val, &sd);
		  y = yy2yyyy(y);
		  mjd = ymd2mjd(y, 01, 00) + doy + ((double) sec) / 86400.0;

		  if (strstr(type, "STAX") != NULL) 
		    {
		      strcpy(sol.Station[i].name, name);
		      strcpy(sol.Station[i].ptcode, ptcode);
		      strcpy(sol.Station[i].solid, s1);
		      sol.Station[i].deweight = OFF;
                      sol.Station[i].mjd = mjd;
		      sol.Station[i].x = val;
		      sol.Station[i].typeb = OFF;
		      sol.Station[i].n_typeb = 0.;
		      sol.Station[i].e_typeb = 0.;
		      sol.Station[i].h_typeb = 0.;
		    } 
		  else if (strstr(type, "STAY") != NULL) 
		    {
		      sol.Station[i].y = val;
		    } 
		  else if (strstr(type, "STAZ") != NULL) 
		    {
		      sol.Station[i].z = val;
		      if (sol.velocity_found == OFF)
			i++;
		    }
		  else if (strstr(type, "VELX") != NULL) 
		    {
		      sol.velocity_found = ON;
		      sol.Station[i].vx = val;
		      sol.Station[i].svx = sd;
                    }
		  else if (strstr(type, "VELY") != NULL) 
		    {
		      sol.Station[i].vy = val;
		      sol.Station[i].svy = sd;
                    }
		  else if (strstr(type, "VELZ") != NULL) 
		    {
		      sol.Station[i].vz = val;
		      sol.Station[i].svz = sd;
		      i++;
                    }
		}
	    }
	}
    }

  rewind(fp);

  done = 0;
  while (fgets(line, MAXLINE, fp) != NULL && !done) 
    {
      if (strncmp(line, "+SITE/ID", 8) == 0) 
	{
	  while (fgets(line, MAXLINE, fp) != NULL) 
	    {
	      if (strncmp(line, "-SITE/ID", 8) == 0) 
		{
		  done = 1;
		  break;
		}

	      if (line[0] != '*') 
		{
		  sscanf(line, "%s%s%s", name, ptcode, domes);
                  for (i=0; i<sol.nsta; i++) {
			if (strncmp(name,sol.Station[i].name,4) == 0) {
		             strcpy(sol.Station[i].domes,domes);
		             sscanf(line, "%*21c%22c",sol.Station[i].desc);
			}
		   }
		}
	    }
	}
    }

  rewind(fp);

  done = 0;
  while (fgets(line, MAXLINE, fp) != NULL && !done) 
    {
      if (strncmp(line, "+SOLUTION/EPOCHS", 16) == 0) 
	{
	  while (fgets(line, MAXLINE, fp) != NULL) 
	    {
	      if (strncmp(line, "-SOLUTION/EPOCHS", 16) == 0) 
		{
		  done = 1;
		  break;
		}

	      if (line[0] != '*') 
		{
		  sscanf(line, "%s%s%s", name, ptcode, solid);

                  for (i=0; i<sol.nsta; i++) {
			if (strncmp(name,sol.Station[i].name,4) == 0 && strncmp(solid,sol.Station[i].solid,2) == 0) {
		            sscanf(line, "%*s%*s%*s%*s %d:%d:%d", &y, &doy, &sec);
			    sol.Station[i].smjd = ymd2mjd(yy2yyyy(y),1,0) + (double) doy + (double) sec/86400.;
		            sscanf(line, "%*s%*s%*s%*s%*s %d:%d:%d", &y, &doy, &sec);
			    sol.Station[i].emjd = ymd2mjd(yy2yyyy(y),1,0) + (double) doy + (double) sec/86400.;
		            sscanf(line, "%*s%*s%*s%*s%*s%*s %d:%d:%d", &y, &doy, &sec);
			    sol.Station[i].mmjd = ymd2mjd(yy2yyyy(y),1,0) + (double) doy + (double) sec/86400.;
			}
                    }
		}
	    }
	}
    }

  rewind(fp);

  if (global.fast == OFF) {
	  done = 0;
	  while (fgets(line, MAXLINE, fp) != NULL && !done) 
	    {
	      if (strncmp(line, "+SOLUTION/MATRIX_ESTIMATE", 25) == 0) 
		{
		  while (fgets(line, MAXLINE, fp) != NULL) 
		    {
		      if (strncmp(line, "-SOLUTION/MATRIX_ESTIMATE", 25) == 0) 
			{
			  done = 1;
			  break;
			}

		      if (line[0] != '*') 
			{
			  i = sscanf(line, "%d %d %lf %lf %lf", &row, &col, &q1, &q2, &q3);
			  if (i >= 3) 
			    {
			      sol.Q[row - 1][col - 1] = q1;
			      sol.Q[col - 1][row - 1] = q1;
			    }

			  if (i >= 4) 
			    {
			      sol.Q[row - 1][col] = q2;
			      sol.Q[col][row - 1] = q2;
			    }

			  if (i >= 5) 
			    {
			      sol.Q[row - 1][col + 1] = q3;
			      sol.Q[col + 1][row - 1] = q3;
			    }
			}
		    }
		}
	    }
        }

  return sol;
}


/*
 * write_sinex()
 *	- assumes solution has been reduced to those station required for SINEX only
 */
void write_sinex(SOLUTION sol, char *outfile)
{
  int i, ii, j, k, y, doy, y2, doy2, nep, deg_lat, min_lat, deg_lon, min_lon;
  double lat, lon, hgt, sec_lat, sec_lon;
  char errstr[MAXLINE];
  FILE *fp;

  fprintf(stderr, "+ writing SINEX file: %s\n", outfile);

  mjd2ydoy(&y, &doy, global.mjd);

  if ((fp = fopen(outfile, "w")) == NULL) 
    {
      sprintf(errstr, "UNABLE TO OPEN %s", outfile);
      fatal_error(errstr);
    }

  nep = sol.npar;
  mjd2ydoy(&y,&doy,sol.smjd);
  mjd2ydoy(&y2,&doy2,sol.emjd);

  if (sol.velocity_found == ON) 
	ii = 6;
  else 
        ii = 3;

  //
  // SINEX header
  //
  if (sol.velocity_found == ON) {
    fprintf(fp, "%%=SNX %s AUS %2.2d:%3.3d:%5.5d AUS %2.2d:%3.3d:%5.5d %2.2d:%3.3d:%5.5d %s %5.5d %d %s\n",
            sol.snxVersion, 0, 0, 0, y, doy, 0, y2, doy2, 0, "C", nep, 2, "X V");
  } else {
    fprintf(fp, "%%=SNX %s AUS %2.2d:%3.3d:%5.5d AUS %2.2d:%3.3d:%5.5d %2.2d:%3.3d:%5.5d %s %5.5d %d %s\n",
            sol.snxVersion, 0, 0, 0, y, doy, 0, y2, doy2, 0, "C", nep, 2, "X  ");
  }
  
  fprintf(fp, "*-------------------------------------------------------------------------------\n");

  //
  // FILE/REFERENCE
  ///
  //fprintf(fp, "+FILE/REFERENCE                                                     \n");
  //fprintf(fp, " DESCRIPTION                                                        \n");
  //fprintf(fp, " OUTPUT             SSC SINEX                                       \n");
  //fprintf(fp, " CONTACT            ________________________________________________\n");
  //fprintf(fp, " SOFTWARE           sinex2epoch version %.2lf                       \n", VERSION);
  //fprintf(fp, " HARDWARE           ________________________________________________\n");
  //fprintf(fp, " INPUT                                                              \n");
  //fprintf(fp, "-FILE/REFERENCE                                                     \n");

  //fprintf(fp, "*-------------------------------------------------------------------------------\n");

  //
  // FILE/COMMENT 
  ///
  fprintf(fp, "+FILE/COMMENT                                                       \n");
  for (i = 0; i < sizeof(sol.snxComments)/sizeof(sol.snxComments[0]); i++)
  {
     if(strncmp(sol.snxComments[i], "*", 1) == 0)
     {
        fprintf(fp, "%s\n", sol.snxComments[i]);
     }
  }
  fprintf(fp, "* File created by sinex2epoch software by John Dawson Geoscience Australia          \n");
  fprintf(fp, "-FILE/COMMENT                                                       \n");

  fprintf(fp, "*-------------------------------------------------------------------------------\n");

  //
  // SITE/ID
  ///
  fprintf(fp, "+SITE/ID\n");
  fprintf(fp, "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_\n");
  for (i = 0; i < sol.nsta; i++) 
    {
	      xyz2plh(&lat, &lon, &hgt, sol.Station[i].x, sol.Station[i].y, sol.Station[i].z, AE, FLAT);
	      deg2dms(lat, &deg_lat, &min_lat, &sec_lat);
	      deg2dms(lon, &deg_lon, &min_lon, &sec_lon);

	      if (min_lat < 0)
		{
		  min_lat *= -1;
		}

	      if (min_lon < 0)
		{
		  min_lon *= -1;
		}

	      if (sec_lat < 0.)
		{
		  sec_lat *= -1.;
		}

	      if (sec_lon < 0.)
		{
		  sec_lon *= -1.;
		}


	      fprintf(fp, " %4s %2s %9s %s %-22s ",
		      sol.Station[i].name, sol.Station[i].ptcode,
		      sol.Station[i].domes, "P", sol.Station[i].desc);

	      fprintf(fp, "%3d %2d %4.1lf %3d %2d %4.1lf %7.1lf\n",
		      deg_lon, min_lon, sec_lon, deg_lat, min_lat, sec_lat, hgt);
    }
  fprintf(fp, "-SITE/ID\n");
  
  fprintf(fp, "*-------------------------------------------------------------------------------\n");

  //
  // SOLUTION/EPOCHS
  ///
  fprintf(fp, "+SOLUTION/EPOCHS\n");
  fprintf(fp, "*Code PT SOLN T Data_start__ Data_end____ Mean_epoch__\n");
  for (i = 0; i < sol.nsta; i++) 
    {
      mjd2ydoy(&y,&doy,sol.Station[i].smjd);
      fprintf(fp, " %4s %2s %4s %1s %2.2d:%3.3d:%5.5d", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid, "C", y, doy, 0);
      mjd2ydoy(&y,&doy,sol.Station[i].emjd);
      fprintf(fp, " %2.2d:%3.3d:%5.5d",y,doy,0);
      mjd2ydoy(&y,&doy,sol.Station[i].mmjd);
      fprintf(fp, " %2.2d:%3.3d:%5.5d\n",y,doy,0);
    }
  fprintf(fp, "-SOLUTION/EPOCHS\n");

  fprintf(fp, "*-------------------------------------------------------------------------------\n");

  //
  // SOLUTION/ESTIMATE
  ///
  fprintf(fp, "+SOLUTION/ESTIMATE\n");
  fprintf(fp, "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___\n");
  for (i = 0, k = 1; i < sol.nsta; i++) 
    {

              mjd2ydoy(&y,&doy,sol.Station[i].mjd);              
	      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
		      k, "STAX", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
		      y, doy, 0, "m", 2, sol.Station[i].x, sqrt(sol.Q[i*ii][i*ii]));
	      k++;

	      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
		      k, "STAY", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
		      y, doy, 0, "m", 2, sol.Station[i].y, sqrt(sol.Q[i*ii+1][i*ii+1]));
	      k++;

	      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
		      k, "STAZ", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
		      y, doy, 0, "m", 2, sol.Station[i].z, sqrt(sol.Q[i*ii+2][i*ii+2]));
	      k++;

	      if (sol.velocity_found == ON) {
		      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
			      k, "VELX", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
			      y, doy, 0, "m/y", 2, sol.Station[i].vx, sqrt(sol.Q[i*6+3][i*6+3]));
		      k++;

		      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
			      k, "VELY", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
			      y, doy, 0, "m/y", 2, sol.Station[i].vy, sqrt(sol.Q[i*6+4][i*6+4]));
		      k++;

		      fprintf(fp, " %5d %-6s %4s %2s %4s %2.2d:%3.3d:%5.5d %-4s %1d %21.14le %10.5le\n",
			      k, "VELZ", sol.Station[i].name, sol.Station[i].ptcode, sol.Station[i].solid,
			      y, doy, 0, "m/y", 2, sol.Station[i].vz, sqrt(sol.Q[i*6+5][i*6+5]));
		      k++;
              }
    }
  fprintf(fp, "-SOLUTION/ESTIMATE\n");

  fprintf(fp, "*-------------------------------------------------------------------------------\n");

  if (global.fast == OFF) {

	  //
	  // +SOLUTION/MATRIX_ESTIMATE U COVA             
	  ///
	  fprintf(fp, "+SOLUTION/MATRIX_ESTIMATE U COVA\n");
      fprintf(fp, "*PARA1 PARA2 ____PARA2+0__________ ____PARA2+1__________ ____PARA2+2__________\n");
	  for (i = 0; i < nep; i++) 
	    {
	      for (j = i; j < nep; j++) 
		{
		  fprintf(fp, " %5d %5d %21.14le ", i + 1, j + 1, sol.Q[i][j]);

		  j++;
		  if (j < nep)
		    {
		      fprintf(fp, "%21.14le ", sol.Q[i][j]);
		    }

		  j++;
		  if (j < nep)
		    {
		      fprintf(fp, "%21.14le", sol.Q[i][j]);
		    }

		  fprintf(fp, " \n");
		}
	    }
        }

  fprintf(fp, "-SOLUTION/MATRIX_ESTIMATE U COVA\n");

  fprintf(fp, "%sENDSNX\n", "%");
}

void fatal_error(char *errstr){
	fprintf(stderr,"Error: %s", errstr);
	exit(0);
}


/*
 * read_commandline()
 */
int read_commandline(int argc, char **argv)
{
    int i, snxfile=OFF;

    global.verbose = OFF;
    global.euler = OFF;
    global.fast = OFF;
    global.vel_reset_then_end = OFF;
    global.varf = VCV_VARF;
    global.reference_epoch_year = 20;
    global.reference_epoch_doy = 1;
    global.ntypeb = NTYPEB;
    global.etypeb = ETYPEB;
    global.utypeb = UTYPEB;
    global.euler_vcv = OFF;
    global.q11 = 0.0;
    global.q12 = 0.0;
    global.q13 = 0.0;
    global.q22 = 0.0;
    global.q23 = 0.0;
    global.q33 = 0.0;

    while (--argc > 0) {

	if ((*++argv)[0] == '-') {
	    while (*++argv[0]) {
		switch (*argv[0]) {

		case 'b':
		    ++argv[0];
                    i = sscanf(*argv, "%lf:%lf:%lf", &global.ntypeb, &global.etypeb, &global.utypeb);
                    if (i != 3) usage();
                    *argv += strlen(*argv) - 1;
		    break;

		case 'e':
		    ++argv[0];
                    i = sscanf(*argv, "%d:%d", &global.reference_epoch_year,&global.reference_epoch_doy);
                    if (i != 2) usage();
                    *argv += strlen(*argv) - 1;
		    break;

		case 'E':
		    global.euler = ON;
		    ++argv[0];
                    i = sscanf(*argv, "%lf:%lf:%lf", &global.rx, &global.ry, &global.rz);
                    if (i != 3) usage();
                    *argv += strlen(*argv) - 1;
		    break;

		case 'f':
		    global.fast = ON;
		    break;

		case 'h':
		    usage();
		    break;

		case 's':
		    ++argv[0];
                    i = sscanf(*argv, "%lf", &global.varf);
                    if (i != 1) usage();
                    *argv += strlen(*argv) - 1;
		    break;

		case 'q':
		    ++argv[0];
                    i = sscanf(*argv, "%le:%le:%le:%le:%le:%le", 
				    &global.q11, &global.q12, &global.q13,
				    &global.q22, &global.q23,
				    &global.q33);
		    global.euler_vcv = ON;
                    if (i != 6) usage();
                    *argv += strlen(*argv) - 1;
		    break;


		case 'v':
		    global.verbose++;
		    break;

		case 'x':
   		    global.vel_reset_then_end = ON;
		    global.fast = ON;
		    break;

		case '?':
		    usage();
		    break;

		default:
		    usage();
		}
	    }
	}

	else {
	    snxfile = ON;
	    strcpy(global.infile, *argv);
        }
    }

    if (snxfile == OFF) { usage(); }

    if (global.euler == ON) { fprintf(stderr,"+ Euler Pole rx,ry,tz mas = %lg %lg %lg\n", global.rx, global.ry, global.rz); }
    if (global.fast == ON) { fprintf(stderr,"+ Executing in fast mode (non-rigorous)\n"); }

    if (global.euler_vcv == ON) {
	    fprintf(stderr,"+ Euler Pole Model VCV\n");
	    fprintf(stderr,"+       %12.5le %12.5le %12.5le\n", global.q11,global.q12,global.q13);
	    fprintf(stderr,"+       %12.5le %12.5le %12.5le\n", global.q12,global.q22,global.q23);
	    fprintf(stderr,"+       %12.5le %12.5le %12.5le\n", global.q13,global.q23,global.q33);
    }

    return 1;
}

/* 
 * usage() 
 */
void usage(void)
{
    fprintf(stderr, "                                                                                       \n");
    fprintf(stderr, "usage: sinex2epoch [OPTIONS] sinex.snx                                                 \n");
    fprintf(stderr, "                                                                                       \n");
    fprintf(stderr, "OPTIONS :-                                                                             \n");
    fprintf(stderr, "                                                                                       \n");
    fprintf(stderr, "   -bN:E:U global type B error to be added to coordinates (typeb.dat takes precident)  \n");
    fprintf(stderr, "   -eYY:DDD allows a new reference epoch to be specified                               \n");
    fprintf(stderr, "   -ERX:RY:RZ user specified Euler Pole                                                \n");
    fprintf(stderr, "   -f execute in fast mode which is non-rigorous                                       \n");
    fprintf(stderr, "   -sVARF apply this global variance factor                                            \n");
    fprintf(stderr, "   -qQ11:Q12:Q13:Q21:Q22:Q23:Q31:Q32:Q33 allows input of plate model VCV               \n");
    fprintf(stderr, "   -x execute to velocity reset step then terminate                                    \n");
    fprintf(stderr, "                                                                                       \n");
    exit(0);
}
