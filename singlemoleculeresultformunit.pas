unit SingleMoleculeResultFormUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, ExtCtrls, StdCtrls;

type

  { TFormSMILES }

  TFormSMILES = class(TForm)
    ButtonOK: TButton;
    Image: TImage;
    Label1: TLabel;
    procedure ButtonOKClick(Sender: TObject);
    procedure FormShow(Sender: TObject);
  private

  public

  end;

var
  FormSMILES: TFormSMILES;

implementation

{$R *.lfm}

{ TFormSMILES }

procedure TFormSMILES.FormShow(Sender: TObject);
begin
  Image.Picture.LoadFromFile('tmp_structure1.tmp.png');
  DeleteFile('tmp_structure1.tmp.png');
end;

procedure TFormSMILES.ButtonOKClick(Sender: TObject);
begin
  FormSMILES.Close;
end;

end.
