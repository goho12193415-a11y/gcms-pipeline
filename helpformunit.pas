unit HelpFormUnit;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, ComCtrls, StdCtrls, Types;

type

  { THelpForm }

  THelpForm = class(TForm)
    Memo1: TMemo;
    Memo2: TMemo;
    Memo3: TMemo;
    PageControl: TPageControl;
    Tab1: TTabSheet;
    Tab2: TTabSheet;
    Tab3: TTabSheet;
    procedure FormResize(Sender: TObject);
  private

  public

  end;

var
  HelpForm: THelpForm;

implementation

{$R *.lfm}

{ THelpForm }


procedure THelpForm.FormResize(Sender: TObject);
begin
  Memo3.Width := HelpForm.Width;
  Memo2.Width := HelpForm.Width;
  Memo1.Width := HelpForm.Width;
  PageControl.Width := HelpForm.Width;
  Memo2.Height := HelpForm.Height - 30;
  Memo1.Height := HelpForm.Height - 30;
  Memo3.Height := HelpForm.Height - 30;
  PageControl.Height := HelpForm.Height;
end;


end.
