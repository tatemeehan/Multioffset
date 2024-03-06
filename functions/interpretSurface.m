function [GPR] = interpretSurface(GPR,gainType,pickerType)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    pickerType = 'polarPicker';
else
    % Manual Picking Not Yet Available
    pickerType = 'polarPicker';
end

if nargin < 2
    gainType = 'AGC';
else
    % t-squared Gain Not Yet Available
    gainType = 'AGC';
end
for ii = 1:GPR.MD.nFiles
% Extract Variables
Rad = GPR.D.Radar{ii};
nChan = GPR.Geometry.nChan{ii};
offsets = GPR.Geometry.offset{ii};
[nearOffset,nearIx] = min(offsets);
% [nearOffset,nearIx] = min(abs(offsets-0.5));
if strcmp(gainType,'AGC')
%     for jj = 1:nChan
%         Rad{jj} = AGCgain(Rad{jj},100,2);
%     end
Rad{nearIx} = AGCgain(Rad{nearIx},100,2);
end
[pickX,pickT] = polarPickerT8({Rad{nearIx}});
% Layman's Deconvolution
prompt = ['*',blanks(2),'Perform Layman`s Deconvolution?',blanks(2),'*'];
titleTxt = ['Surface Reflection Picker'];
laymanQuest = questdlg(prompt, titleTxt,'Yes','No','Yes');
if strcmp(laymanQuest,'Yes')
        laymanDecon = inputdlg('  Number of Samples to Subtract from Pick  :',...
        'Layman`s Deconvolution', 1,{'10'});
        laymanDecon = str2num(cell2mat(laymanDecon));
else
    laymanDecon = 0;
end
% Store Picks
surfaceTWT = (pickT{1}'-laymanDecon).*GPR.D.dt{ii};
% NMO Correction
va = 0.3;
surfaceT0 = sqrt(surfaceTWT.^2-(nearOffset.^2)./va.^2);
% Synthesize Surface Picks at all offsets
surfaceTWT = round(sqrt(surfaceT0.^2+(offsets'./va).^2),1);
GPR.D.surfaceTWT{ii} = surfaceTWT;
GPR.D.surfaceT0{ii} = surfaceT0;
GPR.D.surfaceH{ii} = surfaceT0.*va./2;
end
end
