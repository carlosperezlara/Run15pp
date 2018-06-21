#include <iostream>
#include "Analysis.h"
#include "AT_PiZero.h"
#include "AT_EP.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();


  AT_PiZero *tsk = new AT_PiZero();
  tsk->CentralitySelection(0,5);
  TString ssys = argv[3];
  if(ssys.Contains("ERT")) {
    ssys = "";
  } else if(ssys.Contains("D0")) {
    tsk->SetDist(7.5); ssys="D0";
  } else if(ssys.Contains("D1")) {
    tsk->SetDist(8.5); ssys="D1";
  } else if(ssys.Contains("A0")) {
    tsk->SetAlpha(0.75); ssys="A0";
  } else if(ssys.Contains("A1")) {
    tsk->SetAlpha(0.85); ssys="A1";
  } else if(ssys.Contains("T0")) {
    tsk->SetTime(4.5); ssys="T0";
  } else if(ssys.Contains("T1")) {
    tsk->SetTime(5.5); ssys="T1";
  }

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("PiZero_EP/out%s/out_%s.root",ssys.Data(),run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  ana->AddTask( tsk );

  AT_EP *tsk2 = new AT_EP();
  ana->AddTask( tsk2 );

  ana->Run();

}
