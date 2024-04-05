#include "time.h"
#include <stdio.h>


/*
 * mjd2ydoy()
 *	- convert a mjd into day of year and year
 */
void mjd2ydoy(int *y, int *doy, double mjd)
{
  int m, d;

  mjd2ymd(y, &m, &d, mjd);
  *doy = mjd - ymd2mjd(*y, 1, 0);
  *y = yyyy2yy(*y);
}


/*
 * mjd2gpsweek()
 */
void mjd2gpsweek(int *gpsweek, int *dow, double mjd)
{
  int y, m, d;
  mjd2ymd(&y, &m, &d, mjd);
  *gpsweek = (mjd - ymd2mjd(1980, 01, 06)) / 7;
  *dow = mjd - ymd2mjd(1980, 01, 06) - *gpsweek * 7;
}


/*
 * ydoy2gpsweek()
 */
void ydoy2gpsweek(int *gpsweek, int *dow, int y, int doy)
{
  int gweek;
  if (y < 80)
    {
      y += 2000;
    }
  else
    {
      y += 1900;
    }

  gweek = (ymd2mjd(y, 01, 00) + doy - ymd2mjd(1980, 01, 06)) / 7;
  *dow = ymd2mjd(y, 01, 00) + doy - ymd2mjd(1980, 01, 06) - gweek * 7;
  *gpsweek = gweek;
}


/*
 * ydoy2ymd()
 */
void ydoy2ymd(int *y, int *m, int *d, int year, int doy)
{
  double mjd;

  if (year < 80)
    {
      year += 2000;
    }
  else
    {
      year += 1900;
    }

  mjd = ymd2mjd(year, 01, 00) + doy;
  mjd2ymd(y, m, d, mjd);

  if (*y >= 2000)
    {
      *y -= 2000;
    }
  else
    {
      *y -= 1900;
    }
}


/*	
 * ymd2jd()
 *	converts between gregorian calender date and julian day number (jd)
 *	note: "julian day numbers run from noon to noon.  thus a calculated
 *	julian day number pertains to the noon in the corresponding calender
 *	date" (explantory supplement).  jd returned as a integer*4 
 */
double ymd2jd(int y, int m, int d)
{
  return (1461 * (y + 4800 + (m - 14) / 12)) / 4 + (367 * (m - 2 - 12 * ((m - 14) / 12))) / 12 -
    (3 * ((y + 4900 + (m - 14) / 12) / 100)) / 4 + d - 32075;
}


/*
 * jd2ymd()
 *	- converts a julian day number into a gregorian calender date.   
 */
void jd2ymd(int *y, int *m, int *d, int jd)
{

  int l, n, i, j;

  l = jd + 68569;
  n = (4 * l) / 146097;
  l = l - (146097 * n + 3) / 4;
  i = (4000 * (l + 1)) / 1461001;
  l = l - (1461 * i) / 4 + 31;
  j = (80 * l) / 2447;
  *d = l - (2447 * j) / 80;
  l = j / 11;
  *m = j + 2 - 12 * l;
  *y = 100 * (n - 49) + i + l;
}


/* 	
 * ymd2mjd()
 *	- converts between gregorian calender date and modified julian day number. 
 */
double ymd2mjd(int y, int m, int d)
{
  int jd;
  double mjd;

  jd = (1461 * (y + 4800 + (m - 14) / 12)) / 4 + (367 * (m - 2 - 12 * ((m - 14) / 12))) / 12 -
    (3 * ((y + 4900 + (m - 14) / 12) / 100)) / 4 + d - 32075;

  mjd = (double) (jd - 2400001.0);	/*
					 * extra 12 hours so day starts at 0:00ut 
					 */

  return (double) mjd;
}


/*	
 * mjd2ymd()
 *	- converts a modified julian day number into a gregorian calender date. 
 */
void mjd2ymd(int *y, int *m, int *d, double mjd)
{

  int jd, l, n, i, j;

  jd = (int) (mjd + 2400001.0);	/*
				 * extra 12 hours so day starts at 0:00ut 
				 */

  l = jd + 68569;
  n = (4 * l) / 146097;
  l = l - (146097 * n + 3) / 4;
  i = (4000 * (l + 1)) / 1461001;
  l = l - (1461 * i) / 4 + 31;
  j = (80 * l) / 2447;
  *d = l - (2447 * j) / 80;
  l = j / 11;
  *m = j + 2 - 12 * l;
  *y = 100 * (n - 49) + i + l;
}


/*	
 * mjd2ymdhms()
 * 	-converts a modified julian day number into a gregorian calender date with hours, minutes, seconds. 
 */
void mjd2ymdhms(int *y, int *m, int *d, int *hr, int *mn, double *sc, double mjd)
{
  mjd2ymd(y, m, d, mjd);
  mjd = (mjd - ymd2mjd(*y, *m, *d));
  *hr = (int) 24.0 *mjd;
  *mn = (int) 60.0 *(24.0 * mjd - *hr);
  *sc = 3600.0 * (24.0 * mjd - *hr - *mn / 60.0);
}


/* 
 * yyyy2yy()
 *	- reverse y2k conversion ie. 1998 >> 98
 */
int yyyy2yy(int yyyy)
{
  int yy;

  if (yyyy >= 2000)
    {
      yy = yyyy - 2000;
    }
  else
    {
      yy = yyyy - 1900;
    }

  return yy;
}


/* 
 * yy2yyyy()
 *	- y2k conversion ie. 98 >> 1998
 */
int yy2yyyy(int yy)
{
  int yyyy;

  if (yy < 41)
    {
      yyyy = yy + 2000;
    }
  else
    {
      yyyy = yy + 1900;
    }

  return yyyy;
}


/*
 * yds2yyyy()
 */
int yds2yyyy(int yy, int d, int s)
{
  int yyyy;

  if (yy == 0 && d == 0 && s == 0)
    {
      yy = 40;
    }

  if (yy < 41)
    {
      yyyy = yy + 2000;
    }
  else
    {
      yyyy = yy + 1900;
    }

  return yyyy;
}
