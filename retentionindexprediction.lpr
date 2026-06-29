program retentionindexprediction;

{$mode objfpc}{$H+}

uses {$IFDEF UNIX} {$IFDEF UseCThreads}
  cthreads, {$ENDIF} {$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms,
  tachartlazaruspkg,
  MainFormUnit,
  ResultsFormUnit,
  UtilsUnit,
  SingleMoleculeResultFormUnit,
  HelpFormUnit,
  ChartUnit { you can add units after this };

{$R *.res}

begin
  RequireDerivedFormResource := True;
  Application.Scaled := True;
  Application.Initialize;
  Application.CreateForm(TRIPredictionForm, RIPredictionForm);
  Application.CreateForm(TResultsForm, ResultsForm);
  Application.CreateForm(THelpForm, HelpForm);
  Application.CreateForm(TFormSMILES, FormSMILES);
  Application.CreateForm(TChartForm, ChartForm);
  Application.Run;
end.
