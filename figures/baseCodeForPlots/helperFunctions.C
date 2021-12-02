/*******************************************************************/
/* helperFunctions.C                                               */
/* Author: Stephanie Kwan                                          */
/*******************************************************************/

#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#ifndef HELPER_FUNCTIONS_INCL
#define HELPER_FUNCTIONS_INCL

/*******************************************************************/

// Get the number of entries in an nEvents histogram, at the file
// rootFileDirectory and the folder/hist name histPath.
double getEvents(TString rootFileDirectory, TString histPath)
{
  /* Load file */
  TFile *file = new TFile(rootFileDirectory);
  if (!file->IsOpen() || file==0 )
    {
      std::cout<<"ERROR FILE "<< rootFileDirectory <<" NOT FOUND; EXITING"<<std::endl;
      return 0;
    }
  
  TH1F* nEvents = (TH1F*) file->Get(histPath);
  return (double) nEvents->GetEntries();

}

/*******************************************************************/

#endif
