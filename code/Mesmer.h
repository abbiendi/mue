#ifndef Mesmer_H
#define Mesmer_H

extern "C" void init_mesmer(char* inputdatacard);

extern "C" void mesmer_setup
(int* sampletag,
 char* mesmerversion,
 char* hostname,
 char* datetime,
 int* idproc,
 long int* nev,
 int* areweighted,
 char* RCorder,
 int* includedfs,
 int* radmu,
 int* rade,
 int* iseed1,
 int* iseed2,
 double* emulab,
 double* spread,
 int* extmubeam,
 double* Qmu,
 double* mumass,
 double* elmass,
 double* invalpha,
 double* wnorm,
 double* wmax,
 double* eemin,
 double* themin,
 double* themax,
 double* themumin,
 double* themumax,
 double* ethr,
 double* ththr,
 int* iacopl,
 double* acopl,
 int* iela,
 double* ela,
 int* ivpwgts,
 int* ihadon,
 int* ivpfl,
 int* nwarmup,
 int* ndistrw,
 int* isync,
 double* k0,
 double* phmass);

extern "C" void generate_event_mesmer(double* pmu, int *nfs, int *mcids, double (*pmat)[4], double *weight,
				      int *itag, long int *ievtnr, double *wnovp, double *wnohad, double *wLO,
				      double *wNLO, double *cwvp, int *ierr);

extern "C" void finalize_mesmer(double *xsw, double *exsw, long int *foohpm, long int *fooh, 
				double *truemax, long int *nabove, long int *nlt0, 
				double *xsbias, double *exsbias, double *xsbiasn, double *exsbiasn, 
				double *sumow, double *sum2ow2, double *sumnow, double *sum2now2);

extern "C" void IncomingMuonMomentum_mesmer(double* pmu);

#endif
