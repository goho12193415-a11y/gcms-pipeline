unit UtilsUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Process, FileUtil;

const
  JARNAME = 'retentionprediction4-0.0.6-jar-with-dependencies.jar';
  ClassName = 'ru.ac.phyche.gcmsburyak.retentionprediction4.App2';
  COLS = 6;
  NON_POLAR = 0;
  SEMI_NON_POLAR = 1;
  POLAR = 2;
  MIN_NUM_POL = 15;
  COLUMN_WEIGHTS_DB1701: array[0..4] of double = (0.0, 0.0, 0.0, 0.949, 0.0769);
  COLUMN_WEIGHTS_DB210: array[0..4] of double = (-3.3267, -0.3007, 3.0737, 1.72, -0.1071);


type
  TMultiMolResult = class
  public
    results: TStrings;
    debug: TStrings;
    grid: array[0..COLS - 1] of TStrings;
    constructor Create;
    destructor Destroy; override;
  end;

function predictUsingNNForFile(filename: string; isPolar: boolean;
  columnNum: integer): TMultiMolResult;
function predictUsingNNForFile(filename: string; columnNum: integer;
  columnType: integer): TMultiMolResult;
function parseOutput(output: TStringArray; nNumCols: integer): TMultiMolResult;
function predict5ColumnsForFile(filename: string): TMultiMolResult;
function predictSecondLevelForFile(filename: string; secondLevelModelFile: string;
  secondLevelModelType: integer): TMultiMolResult;
function predictDB624ForFile(filename: string): TMultiMolResult;
function predictDB17ForFile(filename: string): TMultiMolResult;
function predictNonCommonColumnForFile(filename: string;
  weights: array of double): TMultiMolResult;
function testSMILESAndCreateImage(SMILES: string; var riDBWAX: double): boolean;
procedure createTMPFileWithSmiles(SMILES: string);
function removeFileAndCheck(filename: string): boolean;
procedure trainSecondLevel(TrainSetFile: string; ValidationSetFile: string;
  ModelFile: string; ModelType: integer; FilterTrain: boolean);
function ValidateSecondLevel(ValidationSetFile: string; ModelFile: string;
  ModelType: integer): TMultiMolResult;
function tryRunJava: boolean;

implementation

constructor TMultiMolResult.Create();
var
  i: integer;
begin
  results := TStringList.Create;
  debug := TStringList.Create;
  for i := 0 to COLS - 1 do
  begin
    grid[i] := TStringList.Create;
  end;
end;

destructor TMultiMolResult.Destroy;
var
  i: integer;
begin
  results.Free;
  debug.Free;
  for i := 0 to COLS - 1 do
  begin
    grid[i].Free;
  end;
  inherited;
end;

function parseOutput(output: TStringArray; nNumCols: integer): TMultiMolResult;
var
  i, j: integer;
  results: TMultiMolResult;
  found: boolean;
  split: TStringArray;
begin
  results := TMultiMolResult.Create;
  found := False;
  for i := 0 to (Length(output) - 1) do
  begin
    results.debug.Add(output[i]);
    if Trim(output[i]) = 'Predictions:' then
    begin
      found := True;
    end;
    split := output[i].Split([' ']);
    if found and ((Length(split) = nNumCols + 1) or
      ((split[1] = 'ERROR') and ((Length(split) = 2)))) then
    begin
      results.results.Add(output[i]);
      for j := 0 to COLS - 1 do
      begin
        if (j <= 1) then
        begin
          results.grid[j].Add(split[j]);
        end;
        if ((j > 1) and (j <= nNumCols)) then
        begin
          if (split[1] = 'ERROR') then
          begin
            results.grid[j].Add('');
          end
          else
          begin
            results.grid[j].Add(split[j]);
          end;
        end;
        if (j > nNumCols) then
        begin
          results.grid[j].Add('');
        end;
      end;
    end;
  end;
  Result := results;
end;

function predictUsingNNForFile(filename: string; isPolar: boolean;
  columnNum: integer): TMultiMolResult;
var
  results: TMultiMolResult;
  mlp: string;
  cnn: string;
  StringArray: TStringArray;
  s: string;
  i: integer;
begin
  mlp := './models_polar/mlp.nn';
  cnn := './models_polar/cnn.nn';
  if isPolar then
  begin
    mlp := './models_polar/mlpPolar.nn';
    cnn := './models_polar/cnnPolar.nn';
  end;
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'PredictNN',
    mlp, cnn, './models_polar/descriptors_info.txt', IntToStr(columnNum),
    filename], s, i);
  stringArray := s.Split([LineEnding]);
  results := parseOutput(StringArray, 3);
  Result := results;
end;

function predictUsingNNForFile(filename: string; columnNum: integer;
  columnType: integer): TMultiMolResult;
var
  isPolar: boolean;
  columnNumForNN: integer;
begin
  if (columnType = NON_POLAR) then
  begin
    columnNumForNN := columnNum;
    isPolar := False;
  end;
  if (columnType = SEMI_NON_POLAR) then
  begin
    columnNumForNN := columnNum + MIN_NUM_POL;
    isPolar := False;
  end;
  if (columnType = POLAR) then
  begin
    columnNumForNN := columnNum + MIN_NUM_POL;
    isPolar := True;
  end;
  Result := predictUsingNNForFile(filename, isPolar, columnNumForNN);
end;

function predict5ColumnsForFile(filename: string): TMultiMolResult;
var
  StringArray: TStringArray;
  s: string;
  i: integer;
begin
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'Predict5',
    './models_polar/', './models_polar/descriptors_info.txt',
    './models_polar/db624.svr', filename], s, i);
  stringArray := s.Split([LineEnding]);
  Result := parseOutput(StringArray, 5);
end;

function predictSecondLevelForFile(filename: string; secondLevelModelFile: string;
  secondLevelModelType: integer): TMultiMolResult;
var
  StringArray: TStringArray;
  s: string;
  i: integer;
begin
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'PredictSL',
    './models_polar/', './models_polar/descriptors_info.txt',
    secondLevelModelFile, IntToStr(secondLevelModelType), filename], s, i);
  stringArray := s.Split([LineEnding]);
  Result := parseOutput(StringArray, 1);
end;

function predictDB624ForFile(filename: string): TMultiMolResult;
begin
  Result := predictSecondLevelForFile(filename, './models_polar/db624.svr', 0);
end;

function predictDB17ForFile(filename: string): TMultiMolResult;
begin
  Result := predictSecondLevelForFile(filename, './models_polar/db17.svr', 0);
end;

function predictNonCommonColumnForFile(filename: string;
  weights: array of double): TMultiMolResult;
var
  stringArray: TStringArray;
  s: string;
  i, j: integer;
  results: TMultiMolResult;
  found: boolean;
  split: TStringArray;
  ri: double;
begin
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'Predict5',
    './models_polar/', './models_polar/descriptors_info.txt',
    './models_polar/db624.svr', filename], s, i);
  stringArray := s.Split([LineEnding]);
  results := TMultiMolResult.Create;
  found := False;
  for i := 0 to (Length(stringArray) - 1) do
  begin
    results.debug.Add(stringArray[i]);
    if Trim(stringArray[i]) = 'Predictions:' then
    begin
      found := True;
    end;
    split := stringArray[i].Split([' ']);
    if found and ((Length(split) = 6) and (not (split[0] = 'SMILES'))) then
    begin
      ri := 0;
      for j := 0 to 4 do
      begin
        ri := ri + weights[j] * StrToFloat(split[j + 1]);
      end;
      results.results.Add(split[0] + ' ' + FloatToStr(ri));
      for j := 0 to COLS - 1 do
      begin
        if (j = 0) then
        begin
          results.grid[j].Add(split[j]);
        end;
        if (j = 1) then
        begin
          results.grid[j].Add(FloatToStr(ri));
        end;
        if (j > 1) then
        begin
          results.grid[j].Add('');
        end;
      end;
    end;
    if found and ((split[1] = 'ERROR') and ((Length(split) = 2))) then
    begin
      results.results.Add(split[0] + ' ' + split[1]);
      for j := 0 to COLS - 1 do
      begin
        if (j <= 1) then
        begin
          results.grid[j].Add(split[j]);
        end;
        if (j > 1) then
        begin
          results.grid[j].Add('');
        end;
      end;
    end;
  end;
  Result := results;
end;

function removeFileAndCheck(filename: string): boolean;
var
  noFile: boolean;
begin
  noFile := True;
  if FileExists(filename) then
  begin
    noFile := False;
    if DeleteFile(filename) then
    begin
      if not (FileExists(filename)) then
      begin
        noFile := True;
      end;
    end;
  end;
  if not (noFile) then
  begin
    raise Exception.Create('Failed to create and remove temporary file.');
  end;
  Result := noFile;
end;

procedure createTMPFileWithSmiles(SMILES: string);
var
  f: System.Text;
begin
  removeFileAndCheck('tmp_structure1.tmp.txt');

  System.Assign(f, 'tmp_structure1.tmp.txt');
  System.Rewrite(f);
  Write(f, Trim(SMILES));
  System.Close(f);

end;

function testSMILESAndCreateImage(SMILES: string; var riDBWAX: double): boolean;
var
  found: boolean;
  split: TStringArray;
  s: string;
  i: integer;
begin
  riDBWAX := -999.0;
  removeFileAndCheck('tmp_structure1.tmp.png');
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'SMILESToDepiction',
    SMILES, 'tmp_structure1.tmp.png'], s, i);
  if FileExists('tmp_structure1.tmp.png') then
  begin
    RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'PredictNN',
      './models_polar/mlpPolar.nn', './models_polar/cnnPolar.nn',
      './models_polar/descriptors_info.txt', IntToStr(MIN_NUM_POL), SMILES], s, i);
    found := False;
    split := s.Split([LineEnding]);

    for i := 0 to Length(split) do
    begin
      if (split[i].Contains('Prediction Average:')) then
      begin
        found := True;
        riDBWAX := StrToFloat(split[i].Split(' ')[2]);
      end;
    end;
    if found then
    begin
      Result := True;
    end
    else
    begin
      Result := False;
    end;
  end
  else
  begin
    Result := False;
  end;
end;

procedure trainSecondLevel(TrainSetFile: string; ValidationSetFile: string;
  ModelFile: string; ModelType: integer; FilterTrain: boolean);
var
  TrainSetFileF: string;
  ValidationSetFileF: string;
  oldClassName: string;
  s: string;
  i: integer;
begin
  TrainSetFileF := 'tmp_train_set31.ri';
  ValidationSetFileF := 'tmp_validation_set31.ri';
  removeFileAndCheck(TrainSetFileF);
  removeFileAndCheck(ValidationSetFileF);
  CopyFile(TrainSetFile, TrainSetFileF);
  CopyFile(ValidationSetFile, ValidationSetFileF);

  oldClassName := 'ru.ac.phyche.gcmsburyak.retentionprediction4.App';
  RunCommandIndir('.', 'java', ['-cp', JARNAME, oldClassName,
    '-RemoveUnsupportedCompounds', TrainSetFile, TrainSetFileF, 'tmp_10.ri'], s, i);
  RunCommandIndir('.', 'java', ['-cp', JARNAME, oldClassName,
    '-RemoveUnsupportedCompounds', ValidationSetFile, ValidationSetFileF,
    'tmp_10.ri'], s, i);

  if ((not (Trim(ValidationSetFile) = '')) and FilterTrain) then
  begin
    RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName,
      'RemoveOverlap', TrainSetFileF, ValidationSetFileF, TrainSetFileF], s, i);
  end;
  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'Train',
    './models_polar', './models_polar/descriptors_info.txt', ModelFile,
    IntToStr(ModelType), TrainSetFileF], s, i);
  removeFileAndCheck(TrainSetFileF);
  removeFileAndCheck(ValidationSetFileF);
  removeFileAndCheck('tmp_10.ri');
end;

function tryRunJava: boolean;
var
  s: string;
  i: integer;
begin
  i := 0;
  s := '';
  Result := (RunCommandIndir('.', 'java', ['-version'], s, i) = 0);
end;

function ValidateSecondLevel(ValidationSetFile: string; ModelFile: string;
  ModelType: integer): TMultiMolResult;
var
  ValidationSetFileF: string;
  oldClassName: string;
  s: string;
  i: integer;
  results: TMultiMolResult;
  stringArray: TStringArray;
  f: System.Text;
begin
  ValidationSetFileF := 'tmp_validation_set31.ri';
  removeFileAndCheck(ValidationSetFileF);
  CopyFile(ValidationSetFile, ValidationSetFileF);

  oldClassName := 'ru.ac.phyche.gcmsburyak.retentionprediction4.App';
  RunCommandIndir('.', 'java', ['-cp', JARNAME, oldClassName,
    '-RemoveUnsupportedCompounds', ValidationSetFile, ValidationSetFileF,
    'tmp_10.ri'], s, i);

  RunCommandIndir('.', 'java', ['-cp', JARNAME, ClassName, 'Validate',
    './models_polar', './models_polar/descriptors_info.txt', ModelFile,
    IntToStr(ModelType), ValidationSetFileF, 'tmp_10.ri'], s, i);
  removeFileAndCheck(ValidationSetFileF);
  stringArray := s.Split([LineEnding]);
  results := TMultiMolResult.Create;
  for i := 0 to (Length(stringArray) - 1) do
  begin
    results.debug.Add(stringArray[i]);
    if (Trim(stringArray[i])).Contains('RMSE') and
      (Trim(stringArray[i])).Contains('MdAE') then
    begin
      results.results.Add('Validation results:');
      results.results.Add(stringArray[i]);
      results.results.Add('');
      results.results.Add('');
      results.results.Add('');
      results.results.Add('');
      results.results.Add('');
    end;
  end;
  System.Assign(f, 'tmp_10.ri');
  System.Reset(f);
  s := 'SMILES Reference Predicted';
  while (not System.EOF(f)) do
  begin
    stringArray := s.Split(' ');
    results.results.Add(s);
    results.grid[0].Add(stringArray[0]);
    results.grid[1].Add(stringArray[1]);
    results.grid[2].Add(stringArray[2]);
    for i := 3 to COLS - 1 do
    begin
      results.grid[i].Add('');
    end;
    System.readln(f, s);
  end;
  System.Close(f);
  removeFileAndCheck('tmp_10.ri');
  Result := results;
end;



end.
