function [GPR] = pickEditorMx(GPR,isLoadTmpPicks)
%% editPicks.m 
% Scrolls through the data file 1km at a time. After first displaying the 
% radar channels and coherency, the user is prompted to edit picks and
% asked which image to edit on. The manual picker then opens a new window.
% The user defines an editable segment of picks beginning and ending each
% segment with a right mouse click. Picks, including end markers, can be
% placed in any order and will be sorted. The detele and backspace keys
% will delete the previous pick. Pressing the z key enables 10 seconds of
% the zoom tool.
%
% If either radargram channel is selected the picks will snap to the
% neareast zero crossing within a 5 sample window. If coherency is used,
% picks will snap to the maximum within a 7 pixel window, and then laymans
% deconvolution will adjust the peak by -1 ns. The Pick Snapping
% happens during a step of interpolation between the picks -- which is
% displayed in brefiely in real time upon commiting the picks.
%
% Commands: Right Click - Begin Editing Segment
%           Left  Click - Manual Pick
%           z key       - 10 seconds of zoom
%        bckspc/del key -  Delete Previous Pick
%           Right Click - Close Editing Segment
%           Enter Key   - Commit Picks and Advance
%
% In the event of a crash or ungraceful exit, picks upto the current edited
% segment are stored in tmpPicks.mat. To continue from the checkpoint, use
% the logical flag for isLoadTmpPicks = 1 as [GPR] = pickEditor(GPR,1).
% 
% Tate Meehan ~ CRREL 2022
% Notes to self:
% Polarization is not defined, assumes HH = Chan1; HV = Chan2;

% Check inputs
if nargin == 1
    isLoadTmpPicks = 0;
end
for ii = 1:GPR.MD.nFiles
    nChan = GPR.Geometry.nChan{ii};
% ii = 1;
% km get segment Ix
tmpDist = (GPR.Geolocation.Distance{ii} - GPR.Geolocation.Distance{ii}(1))./1000; % [km]
maxDist = max(tmpDist);
rmSeg = 0;
if maxDist > 1.5 % Combine First Two Segments if Seg1 < 1.5 km
    floorDist = floor(maxDist);
    kmArray = 1:floorDist;
    segmentIx = zeros(length(kmArray),1);
    for kk=kmArray
        if kk == kmArray(end)
            if maxDist-kmArray(kk) < 0.5
                segmentIx(kk) = find(tmpDist>=maxDist,1);
            else
                segmentIx(kk)=find(tmpDist>=kmArray(kk),1);
                % Needs Recursion here for  500 > Additional Piece < 1000 m
                segmentIx = [segmentIx;find(tmpDist>=maxDist,1)];
                kmArray = [kmArray,kmArray(kk)+1];
            end
        else
                segmentIx(kk)=find(tmpDist>=kmArray(kk),1);
        end
    end
else
    segmentIx = length(GPR.Geolocation.Distance{ii});
end
% Load Picks from Tempory Save Point
if isLoadTmpPicks
    load([GPR.MD.dataDir,'\tmpPicks.mat'])
    currentMT = find(cellfun(@isempty,picks'),1);
    [row, col] = ind2sub([size(picks,1), size(picks,2)],currentMT);
    currentsegment = col;
    currentChannel = row;
    if isempty(currentsegment)
        currentsegment = length(picks); % Catch for unedited segments
    end
    display(['The current Channel is ', num2str(currentChannel), ' of ', num2str(size(picks,1))])
    display(['The current Segment is ', num2str(currentsegment), ' of ', num2str(size(picks,2))])
    for dum = 1:nChan
        chans{dum} = num2str(dum);
    end
    ChanChoice = menu('Choose a Channel to Resume Pick Editor',chans);
    chanLoopIx = ChanChoice:nChan;
    % Edit a Completed Segment?
    editAns = questdlg('Will you make edits to a previous segment?','Pick Editor','Yes','No','No');
if strcmp(editAns,'Yes')
    prevSegmentIx = 1:currentsegment;
    for ss = 1:currentsegment
        segments{ss} = num2str(prevSegmentIx(ss));
    end
    choice = menu('Choose a Segment Index',segments);
    display(['Pick Editing will continue from the recent checkpoint on segment ',num2str(choice),'.'])
    loopIx = choice:length(segmentIx);
else
    display(['Pick Editing will continue from the recent checkpoint on segment ',num2str(currentsegment),'.'])
    loopIx = find(cellfun(@isempty,picks),1):length(segmentIx);
    choice = 1;
end
else
% Allocate Picks
picks = cell(nChan,length(segmentIx)); pickIx = picks;
loopIx = 1:length(segmentIx);
chanLoopIx = 1:nChan;
end
% Only Works for one Horizon...
for cc = chanLoopIx
for kk = loopIx
    % Extract km Chunks of the Data
    if kk == 1
        xix = 1:segmentIx(kk);
    else
        xix = segmentIx(kk-1)+1:segmentIx(kk);
    end
    % Set the Image TWT window
    twtlim = max(GPR.D.groundTWT{ii}(xix))+.25.*max(GPR.D.groundTWT{ii}(xix));
    if twtlim > max(GPR.D.TimeAxis{ii})
        twtlim = max(GPR.D.TimeAxis{ii});
    end
    
% Multi-offset
isPicks = 1;
if isLoadTmpPicks & kk == choice
    oldPicks = picks{cc,choice};
else
oldPicks = GPR.D.groundTWT{ii}(cc,xix);
end
isResampled = abs((GPR.Geolocation.Distance{ii}(2)-GPR.Geolocation.Distance{ii}(1))-...
        (GPR.Geolocation.Distance{ii}(end)-GPR.Geolocation.Distance{ii}(end-1))) <= 0.1;
figure(100);clf;
% set(gcf,'position',[4 140 1361 472])
set(gcf,'position',[4 223 1361 461])

if isResampled
    p1 = imagesc(GPR.Geolocation.Distance{ii}(xix),GPR.D.TimeAxis{ii},GPR.D.Radar{ii}{cc}(:,xix));colormap(bone);hold on;caxis([quantile(GPR.D.Radar{ii}{cc}(:),[0.005,0.995])]);
else
p1 = pcolor(GPR.Geolocation.Distance{ii}(xix),GPR.D.TimeAxis{ii},GPR.D.Radar{ii}{cc}(:,xix));
shading interp; colormap(bone);axis ij;hold on;caxis([quantile(GPR.D.Radar{ii}{cc}(:),[0.005,0.995])]);
end
if isPicks
 plot(GPR.Geolocation.Distance{ii}(xix),oldPicks,'.','color',[.95,0,1]);
end
% c.FontSize = 12; c.Label.FontSize = 12;
set(gca,'fontsize',12,'fontweight','bold','layer','top')
         title(['Offset ',num2str(round(GPR.Geometry.offset{ii}(cc),2)),' (m)'])
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
ylim([0 twtlim]);
ax = ancestor(p1,'axes');
ax.XAxis.Exponent = 0;

figure(100)
pause(5)
% Begin Picker
editAns = questdlg('Will you make edits here?','Pick Editor','Yes','No','No');
if strcmp(editAns,'Yes')
    % Set Editor Parameters
    signal = GPR.D.Radar{ii}{cc}(:,xix);
    method = 'zero';
    snapWin = 2;

% 1km Image Segments are Brought into Picking Environment
[picks{cc,kk},pickIx{cc,kk}] = ManualPickEditor(signal,oldPicks,xix,method,snapWin,GPR.D.dt{ii},GPR.Geolocation.Distance{ii}(xix),GPR.D.TimeAxis{ii});
% Save Temporary Checkpoint 
save([GPR.MD.dataDir,'\tmpPicks.mat'],'picks','pickIx')
clear('signal','method','snapWin')
close(998);close(100);
else
    % No Edits to Picks
    picks{cc,kk} = GPR.D.groundTWT{ii}(cc,xix);
    pickIx{cc,kk} = xix(:);
end
end
end
uniqueIx = zeros(nChan,max(segmentIx));
newPicks = uniqueIx;
for cc = 1:nChan
    tmp1 = [];
    tmp2 = [];
for kk = 1:length(segmentIx)
    tmp1 = [tmp1;pickIx{cc,kk}(:)];
    tmp2 = [tmp2;picks{cc,kk}(:)];
end

uniqueIx(cc,:) = unique(tmp1);
newPicks(cc,:) = tmp2;
newPicks(cc,:) = newPicks(uniqueIx(cc,:));
end
% Ensure no Picks are Negative
for cc = 1:nChan
for kk = 1:length(newPicks)
zeroPick(kk) = (newPicks(cc,kk)<GPR.D.surfaceTWT{ii}(1,kk));
end
zeroPickIx = find(zeroPick);
newPicks(cc,zeroPickIx) = GPR.D.surfaceTWT{ii}(1,zeroPickIx);
% smooth picks
newPicks(cc,:) = nonParametricSmooth(GPR.Geolocation.Distance{ii},newPicks(cc,:),GPR.Geolocation.Distance{ii},0.5);
end
% Save Picks to Structure
GPR.D.groundTWT{ii} = newPicks';
% D.groundIx{ii} = round(newPicks./GPR.D.dt{ii});
clear('picks','pickIx','uniqueIx','newPicks','oldPicks','ax','c','editAns','p1','twtlim')
end
end

