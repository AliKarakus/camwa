
#include <stdlib.h>
#include "mesh.h"

void meshHaloExtract(mesh_t *mesh, size_t Nbytes, void *sourceBuffer, void *haloBuffer){
  
  // copy data from outgoing elements into temporary send buffer
  for(int i=0;i<mesh->totalHaloPairs;++i){
    // outgoing element
    int e = mesh->haloElementList[i];
    memcpy(((char*)haloBuffer)+i*Nbytes, ((char*)sourceBuffer)+e*Nbytes, Nbytes);
  }
  
}
