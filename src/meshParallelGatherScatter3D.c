#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh3D.h"
#include "mpi.h"

extern "C"
{
  void gsParallelGatherScatter(void *gsh, void *v, const char *type);
}

// ok to use o_v = o_gsv
void meshParallelGatherScatter3D(mesh3D *mesh,
				 ogs_t *ogs, 
				 occa::memory &o_v,
				 occa::memory &o_gsv,
				 const char *type){

  // gather on DEVICE
  mesh->gatherKernel(ogs->Ngather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, ogs->o_gatherTmp);

  // extract gathered halo node data [i.e. shared nodes ]
  if(ogs->Nhalo){
    mesh->getKernel(ogs->Nhalo, ogs->o_gatherTmp, ogs->o_haloLocalIds, ogs->o_haloTmp); // subv = v[ids]
    
    // copy partially gathered halo data from DEVICE to HOST
    ogs->o_haloTmp.copyTo(ogs->haloTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(ogs->gatherGsh, ogs->haloTmp, dfloatString); // danger on hardwired type
    
    // copy totally gather halo data back from HOST to DEVICE
    ogs->o_haloTmp.copyFrom(ogs->haloTmp);
    
    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(ogs->Nhalo, ogs->o_haloTmp,ogs->o_haloLocalIds, ogs->o_gatherTmp); // v[ids] = subv
  }
  
  // do scatter back to local nodes
  mesh->scatterKernel(ogs->Nscatter, ogs->o_scatterOffsets, ogs->o_scatterLocalIds, ogs->o_gatherTmp, o_gsv);
  
}