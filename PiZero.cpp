#include "Analysis.h"
#include "AT_PiZero.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("PiZero/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_PiZero *tsk = new AT_PiZero();
  tsk->DoQA();
  //tsk->CentralitySelection(0,5);
  ana->AddTask( tsk );

  ana->Run();

}
