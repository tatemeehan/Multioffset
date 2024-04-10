function [picks,pickIx] = ManualPickEditor(varargin)
% 
% INPUT:     varargin{1} signal = matrix of input data
%            varargin{2} picks  = horizon of picks needing editing
%            varargin{3} xix    = x indicies of matrix (for segmenting)
%            varargin{4} method = 'min', 'max', or 'zero' for peak snap
%            varargin{5} window = [samples] number of samples to look 
%                                   +- the raw pick for a peak to snap to
%            varargin{6} dt     = [ns] sample rate
%            varargin{7} x      = [m] x axis (Distance or Offset)
%            varargin{8} y      = [ns] y axis (Travel-time)
%            varargin{9} isAGC  = binary option (1 = AGC, 0 = t^2)
% 
% OUTPUT:     picks = [index] interpolated, rounded pick values

signal      = varargin{1};
ntr         = length(signal(1,:)); % number of traces

if length(varargin) == 5
    oldPicks = varargin{2};
    xix     = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;
elseif length(varargin) == 4
    oldPicks = varargin{2};
    xix     = varargin{3};
    method  = varargin{4};
    window  = 3;
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;    
elseif length(varargin) == 2
    oldPicks = varargin{2};
    xix = 1:length(oldPicks);
    window  = 3;
    method  = 'zero';
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;    
elseif length(varargin) == 3
    oldPicks = varargin{2};
    xix = varargin{3};
    window  = 3;
    method  = 'zero';
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;  
elseif length(varargin) == 6
    oldPicks = varargin{2};
    xix     = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = varargin{6};
    x       = 1:ntr;
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;    
elseif length(varargin) == 7
    oldPicks = varargin{2};
        xix = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = varargin{6};
    x       = varargin{7};
    y       = 1:size(signal,1);
    isAGC   = 0;
    L       = 80;    
elseif length(varargin) == 8
    oldPicks = varargin{2};
    xix     = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = varargin{6};
    x       = varargin{7};
    y       = varargin{8};
    isAGC   = 0;
    L       = 80;    
elseif length(varargin) == 9
    oldPicks = varargin{2};
        xix = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = varargin{6};
    x       = varargin{7};
    y       = varargin{8};
    isAGC   = varargin{9};
    L       = 80;   
elseif length(varargin) == 10
    oldPicks = varargin{2};
        xix = varargin{3};
    window  = varargin{5};
    method  = varargin{4};
    dt      = varargin{6};
    x       = varargin{7};
    y       = varargin{8};
    isAGC   = varargin{9};
    L       = varargin{10};
elseif length(varargin) == 1
    error('Not enough input arguments')
end
% Check if pcolor is needed
isResampled = abs((x(2)-x(1))-(x(end)-x(end-1))) <= 0.1;
% Plot the data
fig = figure(998); clf
dcm = datacursormode(fig);
dcm.Enable = 'off';
% removeToolbarExplorationButtons(fig)
c = get(gca,'Children');%clf
set(gcf,'position',[3 365 1361 301])
% Hack for Colormap Determination
% Coherence
if strcmp(method,'max')
% cm = load('yetBlack.txt');
    cm = 'bone';
else
    % Amplitude
    cm = 'bone';
end
try c(length(c)).CData;
catch
    if isAGC % Apply AGC
        signal = AGCgain(signal,L,2);
        if isResampled
            imagesc(x,y,signal);colormap(cm);caxis([quantile(signal(:),[0.005,0.995])]);hold on
        else
        pcolor(x,y,signal);colormap(cm);shading interp; axis ij;caxis([quantile(signal(:),[0.005,0.995])]); hold on
        end
        if length(varargin) >= 8
                ylim([0,max(y)])
        end
    else % t-squared scaling
        if isResampled
            imagesc(x,y,signal);colormap(cm); caxis([quantile(signal(:),[0.005,0.995])]); hold on
        else
        pcolor(x,y,signal);colormap(cm); shading interp; axis ij; caxis([quantile(signal(:),[0.005,0.995])]); hold on
        end
        if length(varargin) >= 8
                ylim([0,max(y)])
        end
    end
    dx = mean(diff(x));
    global caxe 
    global tmpcaxe
    caxe = caxis;
    tmpcaxe = abs(caxe(2));
    plot(x,oldPicks,'.','color',[.95,0,1]);
    [tb,btns] = axtoolbar({'zoomin','zoomout','restoreview'});
end
% Pick the surface
% ginput rules: left mouse click = pick, delete or backspace key = delete
% last pick, z = 10 second zoom , any other key = quit
% Add Rule: Right Click Mark Editable Segment End (Must be Paired)
% Chang rpx to rpx(xix)
% Loop over interpolation, iteration per edit segment
% Addrule 10-19-23: Mark Bad Segment and Remove
button = 1;
kk = 1;
ee = 1;
rpx = []; rpy = [];
epx = []; epy = [];
ax = gca;
disableDefaultInteractivity(ax)
btn1 = uicontrol('style','slider','String','Gain',...
    'FontWeight','bold','FontSize',12,'Position', [15 225 115 25],...
    'Callback',{@displayGain});
set(btn1,'Value',tmpcaxe,'min',tmpcaxe./5,'max',3.*tmpcaxe);
annotation(gcf,"textbox",[.035 .8 .1, .1],'String','Gain','FontWeight','bold','FontSize',12,'FitBoxToText','on',LineStyle='none')
while button == 1 || button == 8 || button == 127 || button == 122 || button == 3 || button == 0 || button == 113
[tmpx, tmpy, button] = ginput(1);
xLimValues = get(gca,'XLim'); 
yLimValues = get(gca,'YLim');
% datacursormode(fig,'off')
if isempty(button)
    if mod(length(epx),2) ~= 0
        warning('You must Pick an additional segment endpoint!')
        button = 0;
        keyboard
    else
        break
    end
elseif (tmpx < xLimValues(1)-dx) | (tmpx > xLimValues(end)+dx) | (tmpy < yLimValues(1)-dt) | (tmpy > yLimValues(end)+dt)
    continue
elseif button == 1
    rpx(kk) = tmpx; rpy(kk) = tmpy;
    plot(rpx(kk),rpy(kk),'rx')
    drawnow
    kk = kk+1;
elseif button == 3
    rpx(kk) = tmpx; rpy(kk) = tmpy;
    epx(ee) = tmpx; epy(ee) = tmpy;
    plot(rpx(kk),rpy(kk),'bx')
    drawnow
    kk = kk+1;
    ee = ee+1;
elseif button == 8 || button == 127 && kk > 1
    c = get(gca,'Children');
    if length(c) > 2
    kk = kk - 1;
    if isempty(epx)
        % Delete Pick and Edit point
        rpx(kk) = []; rpy(kk) = [];
    elseif rpx(kk)==epx(ee-1)
        ee = ee-1;
        % Delete Pick and Edit point
        rpx(kk) = []; rpy(kk) = [];
        epx(ee) = []; epy(ee) = [];
    else
        % Delete Pick
        rpx(kk) = []; rpy(kk) = [];
    end
    % Remove Pick from Plot
    delete(c(1))
    end
elseif button == 122
    figure(998)
    zoom on
    pause(10)
    figure(998)
    zoom off
    continue
elseif button == 113 && kk > 1
    % Surgical Deletion
    c = get(gca,'Children');
    if length(c) > 2
    tmpDist = sqrt((tmpx-rpx).^2 + (tmpy-rpy).^2);
    [~,tmpIx] = min(abs(tmpDist));
    tmpx = rpx(tmpIx);
    tmpy = rpy(tmpIx);
    [isEP,epIx] = ismember(tmpx,epx);
    % Delete Points
    rpx(tmpIx) = [];
    rpy(tmpIx) = [];
    kk = kk-1;
    if isEP && ee > 1
        epx(epIx) = [];
        epy(epIx) = [];
        ee = ee-1;
    end
    % Remove Children From Figure
    for cc = 1:length(c)-2
        cx(cc) = c(cc).XData;
        cy(cc) = c(cc).YData;
    end
    tmpDist = sqrt((tmpx-cx).^2 + (tmpy-cy).^2);
    [~,tmpIx] = min(abs(tmpDist));
    delete(c(tmpIx))
    else
        continue
    end
else
    disp('Not a recognized command. Press the Enter key to commit picks.')
    button = 0;
    continue
end
end
if isempty(rpx)
    warning('No Edits were Made!')
    picks = oldPicks(:);
    pickIx = xix(:);
else
rpx=rpx(:);rpy = rpy(:);
epx = epx(:); epy = epy(:);
nrp = length(rpx);
rp = [rpx,rpy];
% Sort as ascending x values
[sorted,IxS]    = sort(rp(:,1),'ascend');
RawPicks(:,1)   = sorted;
RawPicks(:,2)   = rp(IxS,2);
% Now get only unique values 
[RawX,IxU]      = unique(RawPicks(:,1));
RawY            = RawPicks(IxU,2);
% Sort Editable Ends
if isempty(epx)
    epx(1) = RawX(1);epx(2) = RawX(end);
    epy(1) = RawY(1);epy(2) = RawY(end);
    warning('Segments Were Not Chosen: Defaulting to Picked Region!')
end
ep = [epx,epy];
[sorted,IxS]    = sort(ep(:,1),'ascend');
RawEdits(:,1)   = sorted;
RawEdits(:,2)   = ep(IxS,2);
nEdits = ceil(length(RawEdits(:,1))./2);
for jj = 1:nEdits
    % Bin Segment Edges
%     editx = find(x>=RawEdits(2.*jj-1,1)& x<=RawEdits(2.*jj,1));
%     editIx = find(ismember(editx,x));
    editIx = find(x>=RawEdits(2.*jj-1,1)& x<=RawEdits(2.*jj,1));
            editx = x(editIx);
% 
%     if editIx(1) == 1 || editIx(end) == length(x)
% %         editx = x(editIx);
%     else
%         editx = x(editIx(1)-5:editIx(end)+5);
%     end
%     editIx = find(ismember(x,editx));
    % Bin Picks within Segments
%     binIx = find(RawX>=min(editx) & RawX<=max(editx));
    % Bin Picks including end of Segments
    binIx = 1:length(RawX);
% Fill in
RawPicksAll     = interp1(RawX(binIx),RawY(binIx),editx,'pchip','extrap');
IxN             = find(isnan(RawPicksAll));
IxT             = diff(IxN)-1;
mid             = find(IxT>0);
RawPicksAll(1:IxN(mid))     = RawPicksAll(mid+1);
RawPicksAll(IxN(mid+1):end) = RawPicksAll(IxN(mid+1)-1);

% Get the actual peaks within the window
% Change picks to number of traces in segment, and concatenate segments
newpicks = zeros(length(RawPicksAll),1);
for n = 1:length(RawPicksAll)
    if nargin >= 8
        ixRange        = int64((round(RawPicksAll(n),1)-window.*dt)./dt):...
            int64((round(RawPicksAll(n),1)+window.*dt)./dt);
    else
        ixRange        = int64(round(RawPicksAll(n))-window):...
            int64(round(RawPicksAll(n))+window);
    end
    if any(ixRange <= 0)
        error('Pick Extrapolation went Awry.. Pick More Carefully.')
    end
    % Trave-Time Picks Snapped to Method
    [newpicks(n)]	= snapPicks(signal(:,editIx(n)),RawPicksAll(n),ixRange,dt,window,method);
end
% Overwrite Picks
if strcmp(method,'max')
    % Perform Layman's Deconvolution
    newpicks = newpicks - 10.*dt;
end
% oldPicks(editIx) = RawPicksAll; % This is from Debugging
oldPicks(editIx) = newpicks; % Use the New Picks!
% Interpolate New Picks Smoothly with Existing Neighboring Picks
% joinIx = editIx(1)-3:editIx(end)+3;
% joinPicks    = interp1(RawX(binIx),RawY(binIx),editx,'pchip','extrap');
% Store Picks for Export
picks = oldPicks(:);
pickIx = xix(:);
figure(998)
% Refresh Axis Children
c = get(gca,'Children');
% Remove Raw Picks from figure(998)
delete(c(1:end-1))
% Refresh Axis Children
c = get(gca,'Children');
% Plot Interpolated Picks
plot(x,picks,'.m')
pause(5)
end

end
end
function displayGain(src,evnt)
global tmpcaxe
global caxe
tmpcaxe = round(evnt.Source.Value);
caxe = [-tmpcaxe,tmpcaxe];
figure(998);caxis(caxe);
end