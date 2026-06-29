unit ChartUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, TAGraph,
  TASeries, UtilsUnit, Math;

type

  { TChartForm }

  TChartForm = class(TForm)
    Plot: TChart;
    PlotLineSeries1: TLineSeries;
    procedure FormResize(Sender: TObject);
  private

  public

  end;

var
  ChartForm: TChartForm;

procedure drawChart(Data: TMultiMolResult);

implementation

{$R *.lfm}
procedure drawChart(Data: TMultiMolResult);
var
  i: integer;
begin
  ChartForm.PlotLineSeries1.Clear;
  for i := 1 to Data.grid[0].Count - 1 do
  begin
    ChartForm.PlotLineSeries1.AddXY(StrToFloat(Data.grid[1][i]),
      StrToFloat(Data.grid[2][i]), '', clGreen);
  end;
end;

{ TChartForm }

procedure TChartForm.FormResize(Sender: TObject);
begin
  Plot.Width := Min(ChartForm.Width, ChartForm.Height);
  Plot.Height := Min(ChartForm.Width, ChartForm.Height);
end;

end.
