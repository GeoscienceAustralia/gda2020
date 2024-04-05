/*
 * GRS80
 */
#define   AE  		6378137.0
#define   FLAT  	(1./298.257222101)

/*
 * International
 */
/*
#define   AE  		6378388.0
#define   FLAT  	(1./297.)
*/

#define	  PI    	(4. * atan(1.))
#define	  DEGRAD  	(PI / 180.0)

void vincentypl(double *lat, double *lon, double lat1, double lon1, double azimuth, double dist);
double vincentyd(double x1, double y1, double z1, double x2, double y2, double z2);

void plh2xyz(double *x, double *y, double *z, double glat, double elon, double ht);
void xyz2plh(double *glat, double *elon, double *ht, double x, double y, double z, double ae, double flat);
void uvw2plh(double x, double y, double z, double u, double v, double w, double *lat, double *lon, double *hgt);
void uvw2neh(double x, double y, double z, double u, double v, double w, double *n, double *e, double *h);
void neh2uvw(double x, double y, double z, double n, double e, double h, double *u, double *v, double *w);

void deg2dms(double deg, int *hr, int *mn, double *sc);
void sigma_uvw2neh(double x, double y, double z, double Qxyz[3][3], double Qneh[3][3]);
void sigma_uvw2plh(double x, double y, double z, double Qxyz[3][3], double Qplh[3][3]);
void sigma_neh2uvw(double x, double y, double z, double Qneh[3][3], double Qxyz[3][3]);

