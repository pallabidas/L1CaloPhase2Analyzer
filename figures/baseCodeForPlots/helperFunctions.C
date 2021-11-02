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

double getEvents(TString rootFileDirectory, TString folderName)
{
  /* Load file */
  TFile *file = new TFile(rootFileDirectory);
  if (!file->IsOpen() || file==0 )
    {
      std::cout<<"ERROR FILE "<< rootFileDirectory <<" NOT FOUND; EXITING"<<std::endl;
      return 0;
    }
  
  TH1F* nEvents = (TH1F*) file->Get(folderName + "/nEvents");
  return (double) nEvents->GetEntries();

}

/*******************************************************************/

#endif
