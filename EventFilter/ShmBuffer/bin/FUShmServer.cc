////////////////////////////////////////////////////////////////////////////////
//
// FUShmServer
// -----------
//
//            17/11/2006 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "EventFilter/ShmBuffer/bin/FUShmServer.h"


#include <iostream>


using namespace std;
using namespace evf;


////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
FUShmServer::FUShmServer(FUShmBuffer* buffer)
  : buffer_(buffer)
{
  
}


//______________________________________________________________________________
FUShmServer::~FUShmServer()
{
  
}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
unsigned int FUShmServer::writeNext(unsigned char *data,
				    unsigned int   nFed,
				    unsigned int  *fedSize)
{
  FUShmRawCell* cell =buffer_->rawCellToWrite();
  cell->printState();
  
  if (!cell->isEmpty()) {
    if (cell->isProcessed())
      cout<<"DISCARD "<<cell->index()<<endl;
    else if (cell->isDead())
      cout<<"dead cell "<<cell->index()<<", HANDLE&DISCARD"<<endl;
    else{
      cout<<"ERROR: unexpected state of cell "<<cell->index()<<endl;
      cell->dump();
    }
  }
  
  // write data
  cell->clear();
  unsigned int dataSize(0);
  for (unsigned int i=0;i<nFed;i++) dataSize+=fedSize[i];
  unsigned int   iCell         =cell->index();
  unsigned char *cellBufferAddr=cell->writeData(data,dataSize);
  
  if (0!=cellBufferAddr) {
    // mark feds
    unsigned int fedOffset(0);
    for (unsigned int i=0;i<nFed;i++) {
      unsigned char* fedAddr=cellBufferAddr+fedOffset;
      cell->markFed(i,fedSize[i],fedAddr);
      fedOffset+=fedSize[i];
    }
    
    // set cell state to 'written'
    buffer_->finishWritingRawCell(cell);
    //buffer_->lock();
    //cell->setStateWritten();
    //shmdt(cell);
    //buffer_->unlock();
    //buffer_->postReaderSem();
  }
  
  return iCell;
}
