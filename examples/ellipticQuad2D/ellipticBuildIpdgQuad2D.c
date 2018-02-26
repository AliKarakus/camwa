
#include "ellipticQuad2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgQuad2D(mesh2D *mesh, dfloat tau, dfloat lambda, int *BCType, 
                      nonZero_t **A, int *nnzA, int *globalStarts, const char *options){

  int size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  int *globalIds = (int *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(int));
  
  // every degree of freedom has its own global id
  /* so find number of elements on each rank */
  int *rankNelements = (int*) calloc(size, sizeof(int));
  int *rankStarts = (int*) calloc(size+1, sizeof(int));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_INT,
    rankNelements, 1, MPI_INT, MPI_COMM_WORLD);
  //find offsets
  for(int r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (int e =0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) {
      globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    int *idSendBuffer = (int *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, mesh->Np*sizeof(int), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
    free(idSendBuffer);
  }


  int nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-12.;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj=0;nj<mesh->N+1;++nj){
    for(int ni=0;ni<mesh->N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh->N+1;++j){
        for(int i=0;i<mesh->N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh->Np+node] = 1;
          if(nj==j)
            Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
          if(ni==i)
            Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
          
          ++node;
        }
      }
      
      ++mode;
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));
  
  // reset non-zero counter
  int nnz = 0;
      
  // loop over all elements
  for(int eM=0;eM<mesh->Nelements;++eM){
    
    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          int base = eM*mesh->Np*mesh->Nvgeo + i;
          dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
          dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
          dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
          dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
          dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
          
          int idn = n*mesh->Np+i;
          int idm = m*mesh->Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];
          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            int baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

            // double check vol geometric factors are in halo storage of vgeo
            int idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP    = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo
            int localEP = mesh->vmapP[idM]/mesh->Np;
            int baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
            
            // grab surface geometric factors
            int base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];
            
            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM; 
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
            dfloat lmM = B[idmM];
            
            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP;
            dfloat lmP = B[idmP];
            
            dfloat penalty = tau*hinv;     

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            int eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP < 0) {
              int qSgn, gradqSgn;
              int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
              int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              Anm += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              Anm += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              Anm += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-) 
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
          }
          if(fabs(AnmP)>tol){
            // remote info
            int eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            (*A)[nnz].row = globalIds[eM*mesh->Np + n];
            (*A)[nnz].col = globalIds[eP*mesh->Np + m];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          } 
        }
        if(fabs(Anm)>tol){
          // local block
          (*A)[nnz].row = globalIds[eM*mesh->Np+n];
          (*A)[nnz].col = globalIds[eM*mesh->Np+m];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  free(globalIds);
  free(B);  free(Br); free(Bs); 
}
