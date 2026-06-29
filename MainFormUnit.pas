unit MainFormUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, ComCtrls,
  ExtDlgs, EditBtn, Grids, Types, ResultsFormUnit, UtilsUnit,
  SingleMoleculeResultFormUnit, FileUtil, ChartUnit, HelpFormUnit;

const
  ttnnNN = 0;
  ttnnDB624 = 1;
  ttnnDB17 = 2;
  ttnn1701 = 3;
  ttnn210 = 4;
  ttnnCustom = 5;
  ttnnColumns5 = 6;

type

  { TRIPredictionForm }

  TRIPredictionForm = class(TForm)
    Label11: TLabel;
    ModelComboBox2: TComboBox;
    CustomModelCheckBox: TCheckBox;
    Label10: TLabel;
    ModelFileEdit2: TFileNameEdit;
    Label9: TLabel;
    ValidateButton: TButton;
    TrainButton: TButton;
    DatasetInfoButton: TButton;
    ModelComboBox: TComboBox;
    ModelFileEdit1: TFileNameEdit;
    FilterCheckBox: TCheckBox;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    TestSetFileEdit: TFileNameEdit;
    TrainSetFileEdit: TFileNameEdit;
    PredictButton: TButton;
    SMILESInfoButton: TButton;
    OtherArticlesButton: TButton;
    ColumnTypeComboBox: TComboBox;
    ColumnComboBox: TComboBox;
    SMILESEdit: TEdit;
    SMILESFileEdit: TFileNameEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    PageControl: TPageControl;
    RIPredictionTab: TTabSheet;
    TrainTab: TTabSheet;
    procedure ColumnTypeComboBoxChange(Sender: TObject);
    procedure CustomModelCheckBoxChange(Sender: TObject);
    procedure DatasetInfoButtonClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure OtherArticlesButtonClick(Sender: TObject);
    procedure PredictButtonClick(Sender: TObject);
    procedure SMILESEditChange(Sender: TObject);
    procedure SMILESFileEditChange(Sender: TObject);
    procedure SMILESInfoButtonClick(Sender: TObject);
    procedure TrainButtonClick(Sender: TObject);
    procedure ValidateButtonClick(Sender: TObject);
  private
    procedure DisplayMultiMoleculeResultOnResultsForm(Sender: TObject;
      results: TMultiMolResult);
    procedure Predict(Sender: TObject; n: integer);
    procedure PredcitForFile(Sender: TObject; filename: string; n: integer);
    procedure PredictForSMILES(Sender: TObject; n: integer);
    procedure fileDoesNotExist(Sender: TObject; filename: string);
    procedure fileCheck(Sender: TObject; filenames: TStringArray);
    procedure longOperation(Sender: TObject; operation: string);
  public

  end;

const
  NON_POLAR_COLUMNS: array [0..14] of string =
    ('DB-1', 'SE-30', 'OV-101', 'OV-1', 'Methyl_Silicone',
    'CP_Sil_5_CB', 'BP-1', 'HP-1', 'SPB-1', 'RTX-1', 'Ultra-1',
    'Petrocol_DH', 'Polydimethyl_siloxane',
    'OV-1,_SE-30,_Methyl_silicone,_SP-2100,_OV-101,_DB-1,_etc.', 'Other_non_polar');
  SEMI_NON_POLAR_COLUMNS: array [0..20] of string =
    ('5_%_Phenyl_methyl_siloxane', 'DB-5', 'HP-5',
    'HP-5MS', 'VF-5MS', 'Squalane', 'HP-5_MS', 'SE-54', 'DB-5MS',
    'Apiezon_L', 'BPX-5', 'SE-52', 'CP_Sil_8_CB',
    'SPB-5', 'RTX-5', 'TR5-MS', 'Ultra-2', 'DB-5_MS', 'ZB-5',
    'SLB-5_MS', 'Other_semi_non_polar');
  POLAR_COLUMNS: array [0..20] of string =
    ('DB-Wax', 'Carbowax_20M', 'Supelcowax-10', 'OV-351',
    'HP-Innowax', 'PEG-20M', 'BP-20', 'FFAP', 'CP-Wax_52CB', 'HP-Innowax_FSC',
    'Innowax', 'Innowax_FSC',
    'RTX-Wax', 'PEG_4000', 'DB-FFAP', 'Carbowax', 'ZB-Wax', 'HP-Wax',
    'AT-Wax', 'Stabilwax', 'Other_polar');
  OTHER_COLUMNS: array [0..3] of string = ('DB-624', 'OV-17', 'DB-1701', 'DB-210');

var
  RIPredictionForm: TRIPredictionForm;

implementation

{$R *.lfm}

{ TRIPredictionForm }
procedure TRIPredictionForm.fileDoesNotExist(Sender: TObject; filename: string);
begin
  Dialogs.ShowMessage('Error! File ' + filename + ' does not exist!');
end;

procedure TRIPredictionForm.fileCheck(Sender: TObject; filenames: TStringArray);
var
  i: integer;
begin
  for i := 0 to Length(filenames) - 1 do
  begin
    if (not (FileExists(filenames[i]))) then
    begin
      Dialogs.ShowMessage('Error! File ' + filenames[i] + ' does not exist!');
      Halt(1);
    end;
  end;
end;


procedure TRIPredictionForm.longOperation(Sender: TObject; operation: string);
begin
  Dialogs.ShowMessage(operation + ' can take a very long time for large data sets.' +
    ' It may seem that the program is frozen ' +
    ' and does not respond. Please be patient and wait for a while.');
end;

procedure TRIPredictionForm.ColumnTypeComboBoxChange(Sender: TObject);
begin
  ColumnComboBox.Enabled := True;
  ColumnComboBox.Items.Clear;
  if ColumnTypeComboBox.ItemIndex = 0 then
  begin
    ColumnComboBox.Items.Clear;
    ColumnComboBox.Items.AddStrings(NON_POLAR_COLUMNS);
  end;
  if ColumnTypeComboBox.ItemIndex = 1 then
  begin
    ColumnComboBox.Items.Clear;
    ColumnComboBox.Items.AddStrings(SEMI_NON_POLAR_COLUMNS);
  end;
  if ColumnTypeComboBox.ItemIndex = 2 then
  begin
    ColumnComboBox.Items.Clear;
    ColumnComboBox.Items.AddStrings(POLAR_COLUMNS);
  end;
  if ColumnTypeComboBox.ItemIndex = 3 then
  begin
    ColumnComboBox.Items.Clear;
    ColumnComboBox.Items.AddStrings(OTHER_COLUMNS);
  end;
  ColumnComboBox.ItemIndex := 0;
end;

procedure TRIPredictionForm.SMILESEditChange(Sender: TObject);
begin
  SMILESFileEdit.FileName := '';
end;

procedure TRIPredictionForm.CustomModelCheckBoxChange(Sender: TObject);
begin
  if CustomModelCheckBox.Checked then
  begin
    ColumnTypeComboBox.Enabled := False;
    ColumnComboBox.Enabled := False;
    ModelFileEdit2.Enabled := True;
  end;
  if not (CustomModelCheckBox.Checked) then
  begin
    ColumnTypeComboBox.Enabled := True;
    ColumnComboBox.Enabled := True;
    ModelFileEdit2.Enabled := False;
  end;
end;

procedure TRIPredictionForm.DatasetInfoButtonClick(Sender: TObject);
begin
  HelpForm.Show;
  HelpForm.PageControl.ActivePage := HelpForm.Tab2;
end;

procedure TRIPredictionForm.FormCreate(Sender: TObject);
var
  t: double;
  b: boolean;
begin
  DefaultFormatSettings.DecimalSeparator := '.';
  fileCheck(Sender, ['./models_polar/descriptors_info.txt',
    './models_polar/cnn.nn', './models_polar/cnnPolar.nn', './models_polar/mlp.nn',
    './models_polar/mlpPolar.nn', './models_polar/db624.svr',
    './models_polar/db17.svr', './' + UtilsUnit.JARNAME]);
  t := -100;
  b := False;
  b := tryRunJava;
  if (not (b)) then
  begin
    Dialogs.ShowMessage('Cannot run Java! Probably Java is not installed.' +
      ' Java virtual machine is required!');
  end;
  b := UtilsUnit.testSMILESAndCreateImage('CCCCCCCCCCCC', t);
  if FileExists('tmp_structure1.tmp.png') then
  begin
    DeleteFile('tmp_structure1.tmp.png');
  end
  else
  begin
    b := False;
  end;

  if (not (b)) then
  begin
    Dialogs.ShowMessage(
      'Cannot run the retention index prediction Java package.' +
      ' Check if Java is set up correctly, you have rights to run it, ' +
      'and the ' + UtilsUnit.JARNAME + ' file is available and can be launched.');
    Halt(1);
  end;
  if (t < 1150) or (t > 1250) then
  begin
    Dialogs.ShowMessage('Testing of the retention prediction model is failed!');
    Halt(1);
  end;
end;

procedure TRIPredictionForm.OtherArticlesButtonClick(Sender: TObject);
begin
  HelpForm.Show;
  HelpForm.PageControl.ActivePage := HelpForm.Tab3;
end;

procedure TRIPredictionForm.DisplayMultiMoleculeResultOnResultsForm(Sender: TObject;
  results: TMultiMolResult);
var
  i, j: integer;
begin
  ResultsFormUnit.ResultsForm.Memo1.Clear;
  ResultsFormUnit.ResultsForm.Memo2.Clear;
  ResultsFormUnit.ResultsForm.Grid.Clear;

  for i := 0 to results.debug.Count - 1 do
  begin
    ResultsFormUnit.ResultsForm.Memo2.Lines.Add(results.debug[i]);
  end;
  for i := 0 to results.results.Count - 1 do
  begin
    ResultsFormUnit.ResultsForm.Memo1.Lines.Add(results.results[i]);
  end;
  ResultsFormUnit.ResultsForm.Grid.ColCount := UtilsUnit.COLS;
  ResultsFormUnit.ResultsForm.Grid.RowCount := results.grid[0].Count + 1;
  for i := 0 to results.grid[0].Count - 1 do
  begin
    for j := 0 to UtilsUnit.COLS - 1 do
    begin
      ResultsFormUnit.ResultsForm.Grid.Cells[j, i + 1] := results.grid[j][i];
    end;
  end;
  ResultsFormUnit.ResultsForm.Grid.Cells[0, 0] := 'A';
  ResultsFormUnit.ResultsForm.Grid.Cells[1, 0] := 'B';
  ResultsFormUnit.ResultsForm.Grid.Cells[2, 0] := 'C';
  ResultsFormUnit.ResultsForm.Grid.Cells[3, 0] := 'D';
  ResultsFormUnit.ResultsForm.Grid.Cells[4, 0] := 'E';
  ResultsFormUnit.ResultsForm.Grid.Cells[5, 0] := 'F';
end;

procedure TRIPredictionForm.PredcitForFile(Sender: TObject; filename: string;
  n: integer);
var
  results: TMultiMolResult;
begin
  ResultsForm.Show;
  if (n = ttnnNN) then
  begin
    results := UtilsUnit.predictUsingNNForFile(filename,
      ColumnComboBox.ItemIndex, ColumnTypeComboBox.ItemIndex);
  end;
  if (n = ttnnDB624) then
  begin
    results := UtilsUnit.predictDB624ForFile(filename);
    DisplayMultiMoleculeResultOnResultsForm(Sender, results);
  end;
  if (n = ttnnDB17) then
  begin
    results := UtilsUnit.predictDB17ForFile(filename);
    DisplayMultiMoleculeResultOnResultsForm(Sender, results);
  end;
  if (n = ttnn1701) then
  begin
    results := predictNonCommonColumnForFile(filename, UtilsUnit.COLUMN_WEIGHTS_DB1701);
  end;
  if (n = ttnn210) then
  begin
    results := predictNonCommonColumnForFile(filename, UtilsUnit.COLUMN_WEIGHTS_DB210);
  end;
  if (n = ttnnCustom) then
  begin
    results := UtilsUnit.predictSecondLevelForFile(
      filename, Trim(ModelFileEdit2.Text), ModelComboBox2.ItemIndex);
    DisplayMultiMoleculeResultOnResultsForm(Sender, results);
  end;
  if (n = ttnnColumns5) then
  begin
    results := UtilsUnit.predict5ColumnsForFile(filename);
    DisplayMultiMoleculeResultOnResultsForm(Sender, results);
  end;
  DisplayMultiMoleculeResultOnResultsForm(Sender, results);
  FreeAndNil(results);
end;

procedure TRIPredictionForm.PredictForSMILES(Sender: TObject; n: integer);
begin
  createTMPFileWithSmiles(Trim(SMILESEdit.Text));
  PredcitForFile(Sender, 'tmp_structure1.tmp.txt', n);
  DeleteFile('tmp_structure1.tmp.txt');
  ResultsForm.Memo1.Lines.Add('_____________________________________');
  ResultsForm.Memo1.Lines.Add('Retention index was predicted for this compound.');
end;

procedure TRIPredictionForm.Predict(Sender: TObject; n: integer);
var
  filename: string;
  t: double;
begin
  t := 0;
  if (not (Trim(SMILESFileEdit.Text) = '')) and (Trim(SMILESEdit.Text) = '') then
  begin
    filename := Trim(SMILESFileEdit.Text);
    if FileExists(Trim(filename)) then
    begin
      longOperation(Sender, 'Prediction');
      PredcitForFile(Sender, Trim(filename), n);
    end
    else
    begin
      fileDoesNotExist(Sender, SMILESFileEdit.Text);
    end;
  end;
  if (not (Trim(SMILESEdit.Text) = '')) then
  begin
    if UtilsUnit.testSMILESAndCreateImage(Trim(SMILESEdit.Text), t) then
    begin
      FormSMILES.ShowModal;
      PredictForSMILES(Sender, n);
    end
    else
    begin
      Dialogs.ShowMessage(
        'Incorrect SMILES string or unsupported structure.' +
        ' Elements other than C, O, N, H, Cl, F, Br, I, P, S, Si' +
        ' are not supported. Very large molecules, structures containing ' +
        'two or more fragments that are not linked with' +
        ' covalent bonds are not supported.');
    end;
  end;
end;

procedure TRIPredictionForm.PredictButtonClick(Sender: TObject);
begin
  if not (CustomModelCheckBox.Checked) then
  begin
    if ((ColumnTypeComboBox.ItemIndex > -1) and (ColumnComboBox.ItemIndex > -1)) then
    begin
      if (ColumnTypeComboBox.ItemIndex < 3) then
      begin
        Predict(Sender, ttnnNN);
      end;
      if (ColumnTypeComboBox.ItemIndex = 3) then
      begin
        if (ColumnComboBox.ItemIndex = 0) then
        begin
          Predict(Sender, ttnnDB624);
        end;
        if (ColumnComboBox.ItemIndex = 1) then
        begin
          Predict(Sender, ttnnDB17);
        end;
        if (ColumnComboBox.ItemIndex = 2) then
        begin
          Predict(Sender, ttnn1701);
        end;
        if (ColumnComboBox.ItemIndex = 3) then
        begin
          Predict(Sender, ttnn210);
        end;
      end;
    end
    else
    begin
      Predict(Sender, ttnnColumns5);
    end;
  end
  else
  begin
    if not (FileExists(Trim(ModelFileEdit2.Text))) then
    begin
      fileDoesNotExist(Sender, ModelFileEdit2.Text);
    end
    else
    begin
      Predict(Sender, ttnnCustom);
    end;
  end;
end;

procedure TRIPredictionForm.SMILESFileEditChange(Sender: TObject);
begin
  SMILESEdit.Text := '';
end;

procedure TRIPredictionForm.SMILESInfoButtonClick(Sender: TObject);
begin
  HelpForm.Show;
  HelpForm.PageControl.ActivePage := HelpForm.Tab1;
end;

procedure TRIPredictionForm.TrainButtonClick(Sender: TObject);
var
  failed: boolean;
  i: integer;
begin
  failed := False;
  if not (FileExists(Trim(TrainSetFileEdit.Text))) then
  begin
    failed := True;
    fileDoesNotExist(Sender, TrainSetFileEdit.Text);
  end;
  if Trim(ModelFileEdit1.Text) = '' then
  begin
    failed := True;
  end;
  if (FileExists(Trim(ModelFileEdit1.Text))) then
  begin
    i := 0;
    while (FileExists(Trim(ModelFileEdit1.Text) + '.' + IntToStr(i) + '.bak')) do
    begin
      i := i + 1;
    end;
    CopyFile(Trim(ModelFileEdit1.Text), Trim(ModelFileEdit1.Text) +
      '.' + IntToStr(i) + '.bak', True);
    removeFileAndCheck(Trim(ModelFileEdit1.Text));
  end;
  if ((FilterCheckBox.Checked) and (not (FileExists(Trim(TestSetFileEdit.Text))))) then
  begin
    if (not (Trim(TestSetFileEdit.Text) = '')) then
    begin
      failed := True;
      fileDoesNotExist(Sender, TestSetFileEdit.Text);
    end;
  end;
  if (not (Failed)) then
  begin
    longOperation(Sender, 'Training');
    trainSecondLevel(Trim(TrainSetFileEdit.Text), Trim(TestSetFileEdit.Text),
      Trim(ModelFileEdit1.Text), ModelComboBox.ItemIndex, FilterCheckBox.Checked);
    if (FileExists(Trim(ModelFileEdit1.Text))) then
    begin
      Dialogs.ShowMessage(
        'Model was trained! Testing will be performed now if a test set is specified.');
      ValidateButtonClick(Sender);
    end
    else
    begin
      Dialogs.ShowMessage('Training failed. Check your data set file.');
    end;
  end;
end;

procedure TRIPredictionForm.ValidateButtonClick(Sender: TObject);
var
  failed: boolean;
  results: TMultiMolResult;
begin
  failed := False;
  if not (FileExists(Trim(TestSetFileEdit.Text))) then
  begin
    failed := True;
    fileDoesNotExist(Sender, TestSetFileEdit.Text);
  end;
  if (Trim(ModelFileEdit1.Text) = '') or (Trim(TestSetFileEdit.Text) = '') then
  begin
    failed := True;
  end;
  if not (FileExists(Trim(ModelFileEdit1.Text))) then
  begin
    failed := True;
    fileDoesNotExist(Sender, ModelFileEdit1.Text);
  end;
  if (not (Failed)) then
  begin
    ResultsForm.Show;
    longOperation(Sender, 'Validation');
    results := validateSecondLevel(Trim(TestSetFileEdit.Text),
      Trim(ModelFileEdit1.Text), ModelComboBox.ItemIndex);
    DisplayMultiMoleculeResultOnResultsForm(Sender, results);
    drawChart(results);
    ChartForm.Show;
    FreeAndNil(results);
  end;
end;


end.
