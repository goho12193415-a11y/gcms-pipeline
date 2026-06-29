unit ResultsFormUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, ComCtrls, ExtCtrls,
  CheckLst, Grids, StdCtrls;

type

  { TResultsForm }

  TResultsForm = class(TForm)
    Memo1: TMemo;
    Memo2: TMemo;
    PageControl: TPageControl;
    Grid: TStringGrid;
    TableTab: TTabSheet;
    FullTab: TTabSheet;
    TextTab: TTabSheet;
    procedure FormResize(Sender: TObject);
  private

  public

  end;

var
  ResultsForm: TResultsForm;

implementation

{$R *.lfm}

{ TResultsForm }

procedure TResultsForm.FormResize(Sender: TObject);
begin
  Memo2.Width := ResultsForm.Width;
  Memo1.Width := ResultsForm.Width;
  Grid.Width := ResultsForm.Width;
  PageControl.Width := ResultsForm.Width;
  Memo2.Height := ResultsForm.Height - 30;
  Memo1.Height := ResultsForm.Height - 30;
  Grid.Height := ResultsForm.Height - 30;
  PageControl.Height := ResultsForm.Height;
end;

end.
