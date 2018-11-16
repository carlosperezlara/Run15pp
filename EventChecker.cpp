#include "Analysis.h"
#include "AT_EventChecker.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("eventsOnly/%s.root",run.Data()) );
  ana->OutputFileName( Form("EventChecker/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_EventChecker *tsk = new AT_EventChecker();
  tsk->SetSkipDetails();
  ana->AddTask( tsk );

  ana->Run();

}
