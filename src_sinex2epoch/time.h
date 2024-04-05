
  double ymd2jd(int y, int m, int d);

  void jd2ymd(int *y, int *m, int *d, int jd);

  double ymd2mjd(int y, int m, int d);

  void mjd2ymd(int *y, int *m, int *d, double mjd);

  void mjd2ymdhms(int *y, int *m, int *d, int *hr, int *mn, double *sc, double mjd);

  void mjd2gpsweek(int *gpsweek, int *dow, double mjd);

  void ydoy2gpsweek(int *gpsweek, int *dow, int y, int doy);

  void ydoy2ymd(int *y, int *m, int *d, int year, int doy);

  void mjd2ydoy(int *y, int *doy, double mjd);

  int yyyy2yy(int yyyy);

  int yy2yyyy(int yy);

  int yds2yyyy(int yy, int d, int s);
