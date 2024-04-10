function [GPR] = interpretGround(GPR,gainType,pickerType)
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
    gainType = 'SEC';
end
for ii = 1:GPR.MD.nFiles
% Extract Variables
Rad = GPR.D.Radar{ii};
nChan = GPR.Geometry.nChan{ii};
offsets = GPR.Geometry.offset{ii};
[nearOffset,nearIx] = min(offsets);
if strcmp(gainType,'AGC')
    for jj = 1:nChan
        Rad{jj} = AGCgain(Rad{jj},500,2);
    end
Rad{nearIx} = AGCgain(Rad{nearIx},500,2);
elseif strcmp(gainType,'SEC')
    for jj = 1:nChan
        Rad{jj} = SECgain(Rad{jj},.0005,1);
    end
end
[pickX,pickT] = polarPickerT8(Rad);
% Layman's Deconvolution
prompt = ['*',blanks(2),'Perform Layman`s Deconvolution?',blanks(2),'*'];
titleTxt = ['Layman`s Deconvolution'];
laymanQuest = questdlg(prompt, titleTxt,'Yes','No','Yes');
if strcmp(laymanQuest,'Yes')
        laymanDecon = inputdlg('  Number of Samples to Subtract from Pick  :',...
        'Layman`s Deconvolution', 1,{'10'});
        laymanDecon = str2num(cell2mat(laymanDecon));
else
    laymanDecon = 0;
end
% Store Picks
picks = horzcat(pickT{:});
nHorizon = size(pickT,2);
groundTWT = (picks'-laymanDecon).*GPR.D.dt{ii};
% LQ Filter
% LocalQuantile Filter
X = GPR.Geolocation.X{ii};
Y = GPR.Geolocation.Y{ii};
Distance = GPR.Geolocation.Distance{ii};
r = 5;%.5; % Bin Radius (m)
q = [0.05,0.95]; % Quantile Threshold
for kk = 1:nChan.*nHorizon
[filtPicks, ~] = lqfilter(groundTWT(kk,:),X,Y,r,q);
% Smooth Picks
smoothr = 2.5;%0.5;%2;
[ smoothFiltPicks ] = nonParametricSmooth( Distance,filtPicks,Distance,smoothr );
groundTWT(kk,:) = smoothFiltPicks;
end
GPR.D.groundTWT{ii} = groundTWT;
GPR.D.nHorizon{ii} = nHorizon;
end
end
