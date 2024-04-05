#include <math.h>
#include "geodetic.h"


/*
 * vincentypl()
 *	- input: latitude (deg), longitude (deg), geodetic azimuth (deg), ellipsoidal distance
 *	- output: latitude, longitude
 * 	- Vincenty's direct formulae
 */
void vincentypl(double *lat, double *lon, double lat1, double lon1, double azimuth, double dist)
{
    int iter = 9;
    double tanu1, u_2, tans1, cosa, sina, sins, a, b, be, sigma, sigma0;
    double dsig, cos2sm, tsm, coss, sinu1, cosu1, tans2, t1, tanl, c, omeg;

    be = AE * (1. - FLAT);

    lat1 *= DEGRAD;
    lon1 *= DEGRAD;
    azimuth *= DEGRAD;

    tanu1 = (1 - FLAT) * tan(lat1);
    tans1 = tanu1 / cos(azimuth);
    sina = cos(atan(tanu1)) * sin(azimuth);
    cosa = cos(asin(sina));
    u_2 = cosa * cosa * ((AE * AE) / (be * be) - 1.0);
    a = 1. + (u_2 / 16384.) * (4096. + u_2 * (-768. + u_2 * (320 - 175 * u_2)));
    b = (u_2 / 1024.) * (256. + u_2 * (-128. + u_2 * (74. - 47 * u_2)));

    sigma = (dist / (be * a));

    while (iter--) {
	sigma0 = sigma;
	tsm = 2. * atan(tans1) + sigma;
	cos2sm = cos(tsm);
	sins = sin(sigma);
	coss = cos(sigma);
	dsig = b * sins * (cos2sm + (b / 4.) * (coss * (-1. + 2. * cos2sm * cos2sm) -
						(b / 6.) * cos2sm * (-3. + 4. * sins * sins) * (-3. +
												4. * cos2sm * cos2sm)));
	sigma = (dist / (be * a)) + dsig;

	if (fabs(sigma - sigma0) < 1.e-15)
	    iter = 0;
    }

    sins = sin(sigma);
    coss = cos(sigma);
    sinu1 = sin(atan(tanu1));
    cosu1 = cos(atan(tanu1));
    t1 = (sinu1 * sins - cosu1 * coss * cos(azimuth));
    tans2 = (sinu1 * coss + cosu1 * sins * cos(azimuth)) / ((1. - FLAT) * sqrt(sina * sina + t1 * t1));
    tanl = (sins * sin(azimuth)) / (cosu1 * coss - sinu1 * sins * cos(azimuth));
    c = (FLAT / 16.) * cosa * cosa * (4. + FLAT * (4. - 3. * cosa * cosa));
    omeg =
	atan(tanl) - (1. - c) * FLAT * sina * (sigma + c * sins * (cos(tsm) + c * coss * (-1. + 2. * cos(tsm) * cos(tsm))));

    *lat = atan(tans2) / DEGRAD;
    *lon = (lon1 + omeg) / DEGRAD;

    /*
     * for completeness compute the reverse azimuth
     *
     * tanra = sina / (-sinu1*sins + cosu1*coss*cos(azimuth));
     */
}


/*
 * vincentyd()
 *	- input: cartesian coordinates of two stations
 *	- output: inter-station distance
 * 	- Vincenty's inverse formulae
 */
double vincentyd(double x1, double y1, double z1, double x2, double y2, double z2)
{
    int iter = 9;
    double lam, lam0, sinlam, coslam, tanu1, tanu2, sinu1, cosu1, sinu2, cosu2, u1, u2;
    double sins, coss, sig, sina, cosa, alph, cos2sm, c, be, u_2, a, b, dsig, s;
    double lat1, lon1, ht1, lat2, lon2, ht2;

    if (x1 == x2 && y1 == y2 && z1 == z2)
	return 0.;

    xyz2plh(&lat1, &lon1, &ht1, x1, y1, z1, AE, FLAT);
    xyz2plh(&lat2, &lon2, &ht2, x2, y2, z2, AE, FLAT);

    lat1 *= DEGRAD;
    lon1 *= DEGRAD;
    lat2 *= DEGRAD;
    lon2 *= DEGRAD;

    be = AE * (1. - FLAT);
    tanu1 = (1 - FLAT) * tan(lat1);
    tanu2 = (1 - FLAT) * tan(lat2);
    u1 = atan(tanu1);
    u2 = atan(tanu2);
    sinu1 = sin(u1);
    cosu1 = cos(u1);
    sinu2 = sin(u2);
    cosu2 = cos(u2);
    lam = lon2 - lon1;

    while (iter--) {
	lam0 = lam;
	sinlam = sin(lam);
	coslam = cos(lam);
	coss = sinu1 * sinu2 + cosu1 * cosu2 * coslam;
	sig = acos(coss);
	sins = sin(sig);
	sina = cosu1 * cosu2 * sinlam / sins;
	alph = asin(sina);
	cosa = cos(alph);
	cos2sm = coss - (2. * sinu1 * sinu2 / (cosa * cosa));
	c = (FLAT / 16.) * cosa * cosa * (4. + FLAT * (4. - 3. * cosa * cosa));
	lam = (lon2 - lon1) + (1. - c) * FLAT * sina * (sig + c * sins * (cos2sm + c * coss * (-1. + 2. * cos2sm * cos2sm)));

	if (fabs(lam - lam0) < 1.e-15)
	    iter = 0;
    }

    u_2 = cosa * cosa * ((AE * AE) / (be * be) - 1.0);
    a = 1. + (u_2 / 16384.) * (4096. + u_2 * (-768. + u_2 * (320 - 175 * u_2)));
    b = (u_2 / 1024.) * (256. + u_2 * (-128. + u_2 * (74. - 47 * u_2)));
    dsig =
	b * sins * (cos2sm +
		    (b / 4.) * (coss * (-1. + 2. * cos2sm * cos2sm) -
				(b / 6.) * cos2sm * (-3. + 4. * sins * sins) * (-3. + 4. * cos2sm * cos2sm)));
    s = be * a * (sig - dsig);

    return s;
}

/* 
 * plh2xyz()
 * 	- convert geographic coordinates to geocentric coordinates 
 */
void plh2xyz(double *x, double *y, double *z, double glat, double elon, double ht)
{
    double N, s2, e2;

    glat *= DEGRAD;
    elon *= DEGRAD;
    e2 = 2.0 * FLAT - FLAT * FLAT;
    s2 = sin(glat) * sin(glat);
    N = AE / (sqrt(1.0 - (e2 * s2)));
    *x = cos(glat) * (ht + N) * cos(elon);
    *y = cos(glat) * (ht + N) * sin(elon);
    *z = sin(glat) * (ht + N * (1.0 - e2));
}


/* 
 * xyz2plh()
 *	- convert geocentric coordinates to geographic coordinates 
 */
void xyz2plh(double *glat, double *elon, double *ht, double x, double y, double z, double ae, double flat)
{

    double tol = 1.E-14;
    int max = 10, iter = 0, done = 0;
    double flatfn, funsq, r, e, rsq, rho, sphi, glatr, dht, dlatr, cphi, g1, g2, dr, dz;

    flatfn = (2. - flat) * flat;
    funsq = pow((1. - flat), 2);

    rsq = x * x + y * y;
    r = sqrt(rsq);
    e = atan2(y, x);
    *elon = e / DEGRAD;
    rho = sqrt(z * z + rsq);
    sphi = z / rho;
    glatr = asin(sphi);
    *ht = rho - ae * (1. - flat * sphi * sphi);

    while (!done) {
	sphi = sin(glatr);
	cphi = cos(glatr);
	g1 = ae / sqrt(1. - flatfn * sphi * sphi);
	g2 = g1 * funsq + (*ht);
	g1 += (*ht);
	dr = r - g1 * cphi;
	dz = z - g2 * sphi;
	dht = dr * cphi + dz * sphi;
	*ht += dht;
	dlatr = (dz * cphi - dr * sphi) / (ae + *ht);
	glatr += dlatr;
	iter++;
	if (iter > max)
	    done++;
	if (fabs(dht) / (ae + *ht) <= tol && fabs(dlatr) <= tol)
	    done++;
    }
    *glat = glatr / DEGRAD;
}


/* 
 * uvw2neh()
 * 	- converts a local geocentric vector into E,N,U local vector 
 *	
 *	|n|             -1          |u|
 *	|e| = H(phi) x J(phi,lam) x |v|
 *	|h|                         |w|
 *	
 */
void uvw2neh(double x, double y, double z, double u, double v, double w, double *n, double *e, double *h)
{
    int i, j, k;
    double amat[3][3], coslam, cosphi, eccsq, elon1, glat1, ht1, hmat[3][3], jinv[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);
    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    hmat[i][j] = 0.;

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    hmat[0][0] = radius_m;
    hmat[1][1] = radius_n * cosphi;
    hmat[2][2] = 1.0;

    jinv[0][0] = -sinphi * coslam / radius_m;
    jinv[0][1] = -sinphi * sinlam / radius_m;
    jinv[0][2] = cosphi / radius_m;
    jinv[1][0] = -sinlam / (radius_n * cosphi);
    jinv[1][1] = coslam / (radius_n * cosphi);
    jinv[1][2] = 0.0;
    jinv[2][0] = cosphi * coslam;
    jinv[2][1] = cosphi * sinlam;
    jinv[2][2] = sinphi;

    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    amat[i][j] = 0.;
	    for (k = 0; k < 3; k++)
		amat[i][j] = amat[i][j] + hmat[i][k] * jinv[k][j];
	}
    }

    *n = amat[0][0] * u + amat[0][1] * v + amat[0][2] * w;
    *e = amat[1][0] * u + amat[1][1] * v + amat[1][2] * w;
    *h = amat[2][0] * u + amat[2][1] * v + amat[2][2] * w;
}


/* 
 * uvw2plh()
 * 	- converts a local geocentric vector into latitude, longitude, height vector 
 *	
 *	|p|    -1          |u|
 *	|l| = J(phi,lam) x |v|
 *	|h|                |w|
 *	
 */
void uvw2plh(double x, double y, double z, double u, double v, double w, double *lat, double *lon, double *hgt)
{
    double coslam, cosphi, eccsq, elon1, glat1, ht1, jinv[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);
    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    jinv[0][0] = -sinphi * coslam / radius_m;
    jinv[0][1] = -sinphi * sinlam / radius_m;
    jinv[0][2] = cosphi / radius_m;
    jinv[1][0] = -sinlam / (radius_n * cosphi);
    jinv[1][1] = coslam / (radius_n * cosphi);
    jinv[1][2] = 0.0;
    jinv[2][0] = cosphi * coslam;
    jinv[2][1] = cosphi * sinlam;
    jinv[2][2] = sinphi;

    *lat = jinv[0][0] * u + jinv[0][1] * v + jinv[0][2] * w;
    *lon = jinv[1][0] * u + jinv[1][1] * v + jinv[1][2] * w;
    *hgt = jinv[2][0] * u + jinv[2][1] * v + jinv[2][2] * w;
}


/* 
 * sigma_uvw2neh()
 *
 *	- transform cartesian variance matrix into local system variance matrix
 *	
 *	|n|             -1          |u|
 *	|e| = H(phi) x J(phi,lam) x |v|
 *	|h|                         |w|
 *
 *                         -1                            -1        T
 *      Qneh = (H(phi) x J(phi,lam)) x Qxyz x (H(phi) x J(phi,lam))      
 *
 */
void sigma_uvw2neh(double x, double y, double z, double Qxyz[3][3], double Qneh[3][3])
{
    int i, j, k;
    double tmp[3][3], amat[3][3], coslam, cosphi, eccsq, elon1, glat1, ht1, hmat[3][3], jinv[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);
    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    hmat[i][j] = 0.;

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    hmat[0][0] = radius_m;
    hmat[1][1] = radius_n * cosphi;
    hmat[2][2] = 1.0;

    jinv[0][0] = -sinphi * coslam / radius_m;
    jinv[0][1] = -sinphi * sinlam / radius_m;
    jinv[0][2] = cosphi / radius_m;
    jinv[1][0] = -sinlam / (radius_n * cosphi);
    jinv[1][1] = coslam / (radius_n * cosphi);
    jinv[1][2] = 0.0;
    jinv[2][0] = cosphi * coslam;
    jinv[2][1] = cosphi * sinlam;
    jinv[2][2] = sinphi;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, amat[i][j] = 0.; k < 3; k++)
		amat[i][j] += hmat[i][k] * jinv[k][j];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, tmp[i][j] = 0.0; k < 3; k++)
		tmp[i][j] += amat[i][k] * Qxyz[k][j];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, Qneh[i][j] = 0.; k < 3; k++)
		Qneh[i][j] += tmp[i][k] * amat[j][k];
}


/* 
 * sigma_uvw2plh()
 *
 *	- transform cartesian variance matrix into geodetic system variance matrix
 *	
 *	|p|    -1          |u|
 *	|l| = J(phi,lam) x |v|
 *	|h|                |w|
 *
 *              -1                  -1       T
 *      Qplh = J(phi,lam) x Qxyz x J(phi,lam)      
 *
 */
void sigma_uvw2plh(double x, double y, double z, double Qxyz[3][3], double Qplh[3][3])
{
    int i, j, k;
    double tmp[3][3], coslam, cosphi, eccsq, elon1, glat1, ht1, jinv[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);
    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    jinv[0][0] = -sinphi * coslam / radius_m;
    jinv[0][1] = -sinphi * sinlam / radius_m;
    jinv[0][2] = cosphi / radius_m;
    jinv[1][0] = -sinlam / (radius_n * cosphi);
    jinv[1][1] = coslam / (radius_n * cosphi);
    jinv[1][2] = 0.0;
    jinv[2][0] = cosphi * coslam;
    jinv[2][1] = cosphi * sinlam;
    jinv[2][2] = sinphi;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, tmp[i][j] = 0.0; k < 3; k++)
		tmp[i][j] += jinv[i][k] * Qxyz[k][j];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, Qplh[i][j] = 0.; k < 3; k++)
		Qplh[i][j] += tmp[i][k] * jinv[j][k];
}


/* 
 * neh2uvw()
 * 	- converts a e,n,u local vector into a geocentric vector 
 *	
 *	|u|                 -1      |n|
 *	|v| = J(phi,lam) x H(phi) x |e|
 *	|w|                         |h|
 *
 */
void neh2uvw(double x, double y, double z, double n, double e, double h, double *u, double *v, double *w)
{
    int i, j, k;
    double amat[3][3], coslam, cosphi, eccsq, elon1, glat1, ht1, hinv[3][3], jmat[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);

    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    hinv[i][j] = 0.;

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    hinv[0][0] = 1.0 / radius_m;
    hinv[1][1] = 1.0 / (radius_n * cosphi);
    hinv[2][2] = 1.0;

    jmat[0][0] = -radius_m * coslam * sinphi;
    jmat[0][1] = -radius_n * cosphi * sinlam;
    jmat[0][2] = cosphi * coslam;
    jmat[1][0] = -radius_m * sinphi * sinlam;
    jmat[1][1] = radius_n * coslam * cosphi;
    jmat[1][2] = cosphi * sinlam;
    jmat[2][0] = radius_m * cosphi;
    jmat[2][1] = 0.0;
    jmat[2][2] = sinphi;

    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    amat[i][j] = 0.;
	    for (k = 0; k < 3; k++)
		amat[i][j] += jmat[i][k] * hinv[k][j];
	}
    }

    *u = amat[0][0] * n + amat[0][1] * e + amat[0][2] * h;
    *v = amat[1][0] * n + amat[1][1] * e + amat[1][2] * h;
    *w = amat[2][0] * n + amat[2][1] * e + amat[2][2] * h;
}


/* 
 * sigma_neh2uvw()
 *	
 *	- transform local system variance matrix into cartesian system variance matrix
 *
 *	|u|                 -1      |n|
 *	|v| = J(phi,lam) x H(phi) x |e|
 *	|w|                         |h|
 *                           -1                            -1    T
 *      Quvw = J(phi,lam) x H(phi) x Qneh x (J(phi,lam) x H(phi))     
 */
void sigma_neh2uvw(double x, double y, double z, double Qxyz[3][3], double Qneh[3][3])
{
    int i, j, k;
    double tmp[3][3], amat[3][3], coslam, cosphi, eccsq, elon1, glat1, ht1, hinv[3][3], jmat[3][3];
    double lam1, one_minus_e2sin2phi, phi1, radius_m, radius_n, sinlam, sinphi, sinsqphi;

    eccsq = 2 * FLAT - (FLAT * FLAT);

    xyz2plh(&glat1, &elon1, &ht1, x, y, z, AE, FLAT);

    phi1 = glat1 * DEGRAD;
    lam1 = elon1 * DEGRAD;
    sinphi = sin(phi1);
    cosphi = cos(phi1);
    sinsqphi = sinphi * sinphi;
    sinlam = sin(lam1);
    coslam = cos(lam1);
    one_minus_e2sin2phi = (1. - eccsq * sinsqphi);

    radius_n = AE / sqrt(one_minus_e2sin2phi);
    radius_m = radius_n * (1. - eccsq) / (one_minus_e2sin2phi);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    hinv[i][j] = 0.;

    radius_m = radius_m + ht1;
    radius_n = radius_n + ht1;

    hinv[0][0] = 1.0 / radius_m;
    hinv[1][1] = 1.0 / (radius_n * cosphi);
    hinv[2][2] = 1.0;

    jmat[0][0] = -radius_m * coslam * sinphi;
    jmat[0][1] = -radius_n * cosphi * sinlam;
    jmat[0][2] = cosphi * coslam;
    jmat[1][0] = -radius_m * sinphi * sinlam;
    jmat[1][1] = radius_n * coslam * cosphi;
    jmat[1][2] = cosphi * sinlam;
    jmat[2][0] = radius_m * cosphi;
    jmat[2][1] = 0.0;
    jmat[2][2] = sinphi;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, amat[i][j] = 0.; k < 3; k++)
		amat[i][j] += jmat[i][k] * hinv[k][j];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, tmp[i][j] = 0.0; k < 3; k++)
		tmp[i][j] += amat[i][k] * Qneh[k][j];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    for (k = 0, Qxyz[i][j] = 0.; k < 3; k++)
		Qxyz[i][j] += tmp[i][k] * amat[j][k];
}

/* 
 * deg2dms()
 * 	- convert decimal degrees into degrees, minutes and seconds 
 */
void deg2dms(double deg, int *hr, int *mn, double *sc)
{
    double h, m;
    deg = modf(deg, &h) * 60.0;
    *sc = modf(deg, &m) * 60.0;
    *hr = (int) h;
    *mn = (int) m;
}

