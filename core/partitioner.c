#include <math.h>

#include "gslib.h"
#if defined(PARRSB)
#include "parRSB.h"
#endif

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {long long vtx[MAXNV]; long long eid; int proc;} edata;

#define HEADER_LEN 132

#define writeMapFile FORTRAN_UNPREFIXED(writemapfile,WRITEMAPFILE)
int writeMapFile(int nelt,int nv,int *pmap,int *vtx,char *fileName,
		       MPI_Comm comm){
  int rank,size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

#if !defined(NOMPIIO)
  char *version="#v001";
  int nelgt=nelt;
  MPI_Allreduce(&nelt,&nelgt,1,MPI_INT,MPI_SUM,comm);
  const int npts=nelgt*nv;
  const int depth=(int)log2(1.0*nelgt);
  const int d2=(int)(pow(2,depth)+0.5);
  int nactive=nelgt,nrnk=nelgt,noutflow=0;

  char header[HEADER_LEN];
  header[HEADER_LEN]='\0';
  sprintf(header,"%5s%12d%12d%12d%12d%12d%12d%12d",version,
		  nelgt,nactive,depth,d2,npts,nrnk,noutflow);
  memset(header+strlen(header),' ',HEADER_LEN-strlen(header));
  header[HEADER_LEN]='\0';

  MPI_Info infoIn;
  MPI_Info_create(&infoIn);
  MPI_Info_set(infoIn,"access_style","write_once,random");

  int errs=0;
  MPI_File file;
  int err=MPI_File_open(comm,fileName,
		MPI_MODE_WRONLY|MPI_MODE_CREATE,
		infoIn,&file);
  if(err){
    errs++;
    MPI_Abort(comm,911);
  }

  int writeSize=0;
  if(rank==0)
    writeSize=HEADER_LEN*sizeof(char)+sizeof(float);
  writeSize+=(nv+1)*nelt*sizeof(int);

  char *buf=(char*)malloc(writeSize),*buf0=buf;

  float test=6.54321;
  if(rank==0){
    memcpy(buf0,header,HEADER_LEN*sizeof(char)),buf0+=HEADER_LEN;
    memcpy(buf0,&test,sizeof(float)),buf0+=sizeof(float);
  }

  int i;
  for(i=0;i<nelt;i++){
    memcpy(buf0,&pmap[i],sizeof(int)),buf0+=sizeof(int);
    memcpy(buf0,&vtx[i*nv],sizeof(int)*nv),buf0+=nv*sizeof(int);
  }

  MPI_Status status;
  err=MPI_File_write_ordered(file,buf,writeSize,MPI_BYTE,&status);
  if(err) errs++;

  err=MPI_File_close(&file);
  if(err) errs++;

  MPI_Barrier(comm);
  free(buf);

  return errs;
#else
  if(rank==0) printf("MPI-IO is disabled. Enable MPI-IO to dump .ma2 file.\n");
#endif
}

#ifdef PARMETIS

#include "parmetis.h"
#include "defs.h"

int parMETIS_partMesh(int *part, long long *vl, int nel, int nv, int *opt, comm_ext ce)
{
  int i, j;
  int ierrm;
  double time, time0;

  MPI_Comm comms;
  struct comm comm;
  int color;
  int ibuf;

  struct crystal cr;
  struct array A;
  edata *row;

  long long nell;
  long long *nelarray;
  idx_t *elmdist;
  idx_t *evlptr;
  idx_t *part_;
  real_t *tpwgts;
  idx_t edgecut;
  real_t ubvec;
  idx_t *elmwgt;
  idx_t wgtflag;
  idx_t numflag;
  idx_t ncon;
  idx_t ncommonnodes;
  idx_t nparts;
  idx_t nelsm;
  idx_t options[10];

  ierrm = METIS_OK;
  nell = nel;
  edgecut = 0;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  ubvec = 1.02;
  elmwgt = NULL; /* no weights */
  ncommonnodes = 2;

  part_ = (idx_t*) malloc(nel*sizeof(idx_t));

  if (sizeof(idx_t) != sizeof(long long)){
    printf("ERROR: invalid sizeof(idx_t)!\n");
    goto err;
  }
  if (nv != 4 && nv != 8){
    printf("ERROR: nv is %d but only 4 and 8 are supported!\n", nv);
    goto err;
  }

  color = MPI_UNDEFINED;
  if (nel > 0) color = 1;
  MPI_Comm_split(ce, color, 0, &comms);
  if (color == MPI_UNDEFINED)
    goto end;

  comm_init(&comm,comms);
  if (comm.id == 0)
    printf("Running parMETIS ... "), fflush(stdout);

  nelarray = (long long*) malloc(comm.np*sizeof(long long));
  MPI_Allgather(&nell, 1, MPI_LONG_LONG_INT, nelarray, 1, MPI_LONG_LONG_INT, comm.c);
  elmdist = (idx_t*) malloc((comm.np+1)*sizeof(idx_t));
  elmdist[0] = 0;
  for (i=0; i<comm.np; ++i)
    elmdist[i+1] = elmdist[i] + (idx_t)nelarray[i];
  free(nelarray);

  evlptr = (idx_t*) malloc((nel+1)*sizeof(idx_t));
  evlptr[0] = 0;
  for (i=0; i<nel; ++i)
    evlptr[i+1] = evlptr[i] + nv;
  nelsm = elmdist[comm.id+1] - elmdist[comm.id];
  evlptr[nelsm]--;

  if (nv == 8) ncommonnodes = 4;
  nparts = comm.np;

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = 0;
  options[PMV3_OPTION_SEED]   = 0;
  if (opt[0] != 0) {
    options[PMV3_OPTION_DBGLVL] = opt[1];
    if (opt[2] != 0) {
      options[3] = PARMETIS_PSR_UNCOUPLED;
      nparts = opt[2];
    }
  }

  tpwgts = (real_t*) malloc(ncon*nparts*sizeof(real_t));
  for (i=0; i<ncon*nparts; ++i)
    tpwgts[i] = 1./(real_t)nparts;

  if (options[3] == PARMETIS_PSR_UNCOUPLED)
    for (i=0; i<nel; ++i)
      part_[i] = comm.id;

  comm_barrier(&comm);
  time0 = comm_time();
  ierrm = ParMETIS_V3_PartMeshKway(elmdist,
                                   evlptr,
                                   (idx_t*)vl,
                                   elmwgt,
                                   &wgtflag,
                                   &numflag,
                                   &ncon,
                                   &ncommonnodes,
                                   &nparts,
                                   tpwgts,
                                   &ubvec,
                                   options,
                                   &edgecut,
                                   part_,
                                   &comm.c);

  time = comm_time() - time0;
  if (comm.id == 0)
    printf("%lf sec\n", time), fflush(stdout);

  for (i=0; i<nel; ++i)
    part[i] = part_[i];

  free(elmdist);
  free(evlptr);
  free(tpwgts);
  MPI_Comm_free(&comms);
  comm_free(&comm);

end:
  comm_init(&comm,ce);
  comm_allreduce(&comm, gs_int, gs_min, &ierrm, 1, &ibuf);
  if (ierrm != METIS_OK) goto err;
  return 0;

err:
  return 1;
}
#endif

#define dumpMapFile FORTRAN_UNPREFIXED(dumpmapfile,DUMPMAPFILE)
void dumpMapFile(char *casename,int *len,int *nell,int *nve,int *part,
  long long *el,long long *vl,int *fcomm,int *retval)
{
  struct comm comm;
  *retval=1;
#if defined(MPI)
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  int e, n;
  struct array eList;
  edata *data;

  int nel=*nell;
  int nv=*nve;

  array_init(edata, &eList, nel), eList.n = nel;
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid  = el[e];
    for(n = 0; n < nv; ++n) {
      data[e].vtx[n] = vl[e*nv + n];
    }
  }

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(edata,eList.ptr,(unsigned int)nel,eid,1,&buf);

  int *vli=(int*) malloc(nv*nel*sizeof(int));
  int *parti=(int*) malloc(nel*sizeof(int));
  for(data=eList.ptr, e=0; e<nel; ++e) {
    parti[e]=data[e].proc;
    for(n=0; n<nv; ++n) {
      vli[e*nv+n] = data[e].vtx[n];
    }
  }
  buffer_free(&buf);

  char mapfile[BUFSIZ];
  memcpy(mapfile,casename,sizeof(char)*(*len));
  mapfile[*len]='\0';
  char ext[BUFSIZ];
  sprintf(ext,"_%08d.ma2",comm.np);
  strcat(mapfile,ext);
  writeMapFile(nel,nv,parti,vli,mapfile,comm.c);

  free(parti);
  free(vli);

  array_free(&eList);
  comm_free(&comm);

  *retval=0;
}

#define transferElements FORTRAN_UNPREFIXED(transferelements,TRANSFERELEMENTS)
void transferElements(int *nell,int *nve,int *part,long long *el,long long *vl,
  int *lelt,int *fcomm,int *retval)
{
  *retval=1;

  struct comm comm;
  struct crystal cr;
  int ibuf;

#if defined(MPI)
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  int e, n;
  int count;

  struct array eList;
  edata *data;

  int nel=*nell;
  int nv=*nve;

  array_init(edata, &eList, nel), eList.n = nel;
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].eid  = el[e];
    for(n = 0; n < nv; ++n) {
      data[e].vtx[n] = vl[e*nv + n];
    }
  }

  crystal_init(&cr, &comm);
  sarray_transfer(edata, &eList, proc, 0, &cr);
  crystal_free(&cr);

  nel = eList.n;

  count = 0;
  if (nel > *lelt) count = 1;
  comm_allreduce(&comm, gs_int, gs_add, &count, 1, &ibuf);
  if (count > 0) {
    if (comm.id == 0)
      printf("ERROR: resulting parition requires lelt=%d!\n", nel);
    goto err;
  }

  for(data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].eid;
    for(n = 0; n < nv; ++n) {
      vl[e*nv + n] = data[e].vtx[n];
    }
  }

  array_free(&eList);
  comm_free(&comm);

  *nell = nel;
  *retval=0;

err:
  fflush(stdout);
  *retval = 1;
}

#define fpartMesh FORTRAN_UNPREFIXED(fpartmesh,FPARTMESH)
void fpartMesh(int *part,long long *el,long long *vl,const int *nell,
  const int *nve,int *fcomm,int *rtval)
{
  struct comm comm;

  int nel, nv;
  int ierr;
  int opt[3];

  nel  = *nell;
  nv   = *nve;

#if defined(MPI)
  comm_ext cext = MPI_Comm_f2c(*fcomm);
#else
  comm_ext cext = 0;
#endif
  comm_init(&comm, cext);

  ierr = 1;
#if defined(PARRSB)
  opt[0] = 1;
  opt[1] = 2; /* verbosity */
  opt[2] = 0;
  ierr = parRSB_partMesh(part, vl, nel, nv, opt, comm.c);
#elif defined(PARMETIS)
  opt[0] = 1;
  opt[1] = 0; /* verbosity */
  opt[2] = comm.np;
  ierr = parMETIS_partMesh(part, vl, nel, nv, opt, comm.c);
#endif
  if (ierr != 0) goto err;

  comm_free(&comm);
  *rtval = 0;
  return;

err:
  fflush(stdout);
  *rtval = 1;
}

void printPartStat(long long *vtx, int nel, int nv, comm_ext ce)
{
  int i,j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;

  struct gs_data *gsh;
  int b;

  int numPoints;
  long long *data;

  comm_init(&comm,ce);
  np = comm.np;
  id = comm.id;

  if (np == 1) return;

  numPoints = nel*nv;
  data = (long long*) malloc(numPoints*sizeof(long long));
  for(i = 0; i < numPoints; i++) data[i] = vtx[i];

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  pw_data_nmsg(gsh, &Nmsg);
  Ncomm = (int *) malloc(Nmsg*sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax , 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin , 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum , 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (i=1; i<Nmsg; ++i){
    nsMax = Ncomm[i] > Ncomm[i-1] ? Ncomm[i] : Ncomm[i-1];
    nsMin = Ncomm[i] < Ncomm[i-1] ? Ncomm[i] : Ncomm[i-1];
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax , 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin , 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax , 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin , 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nssSum , 1, &b);

  nsSum = nsSum/Nmsg;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum , 1, &b);

  nelMax = nel;
  nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if (id == 0) {
    printf(
      " Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
      ncMax, ncMin, (double)ncSum/np);
   printf(
      " Max nvolume: %d | Min nvolume: %d | Avg nvolume: %lf\n",
      nsMax, nsMin, (double)nsSum/np);
    printf(
      " Max volume: %d | Min volume: %d | Avg volume: %lf\n",
      nssMax, nssMin, (double)nssSum/np);
    printf(
      " Max elements: %d | Min elements: %d | Balance: %lf\n",
      nelMax, nelMin, (double)nelMax/nelMin);
    fflush(stdout);
  }

  comm_free(&comm);
}

#define fprintPartStat FORTRAN_UNPREFIXED(printpartstat,PRINTPARTSTAT)
void fprintPartStat(long long *vtx, int *nel, int *nv, int *comm)
{

#if defined(MPI)
  comm_ext c = MPI_Comm_f2c(*comm);
#else
  comm_ext c = 0;
#endif

  printPartStat(vtx, *nel, *nv, c);
}


