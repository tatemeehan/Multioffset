function [GPR] = autoPickSurface2(GPR,isQuickLook)
% autoPickSurface.m applies a gain function to the data to illuminate the
% snow surface reflection, while muting out the later reflection arrivals.
% This is accomplished by assuming that the antenna height is roughly
% 1 m above the snow then the gain taper is selected via the NMO equation.
% The picks are made automatically using the Modified Energy Ratio (MER)
% (Wong,2008). Then the zero offset travel-time and antenna height are 
% estimated by least-squares following the x^2-t^2 method of Green (1938).
% The algorithm can be optionally set to apply bootstraping to the
% least-squares estimate. However, due to smoothing applied to the raw
% travel-time picks, little benefit is gained by the Monte Carlo process.
% Lastly, the zero-offset travel-times are corrected to ensure that the 
% EM free-space wave speed of 0.3 m/ns is true. In experiments we found
% that travel-time uncertainties are on the order of +- 2 samples, leading
% to velocity innacuracies of +- 0.01 m/ns

% Input GPR structure array
%       isQuickLook - Logical Option to Plot Surface Picks
% Output GPR struture array with additional variables:
%       GPR.D.surfaceT0 - the estimated zero-offset snow surface traveltime
%       GPR.D.surfaceTWT - the surface reflection travel-time picks 
%       GPR.D.surfaceH - the estimated thickness of the air layer
%   

% Tate Meehan, 1-22-21, SnowEx 2021

% Check inputs
if nargin == 1
    isQuickLook = 0;
end

for ii = 1:GPR.MD.nFiles
    % Allocation
    offsets = GPR.Geometry.offset{ii};
    surfaceTWT = zeros(numel(offsets),size(GPR.D.Radar{ii}{1},2));
    nChan = GPR.Geometry.nChan{ii};
    for jj = 1:nChan
        % Extract Variables
        t0 = GPR.D.t0{ii}(jj);
        dt = GPR.D.dt{ii};
        Rad = GPR.D.Radar{ii}{jj};
        TWT = GPR.D.TimeAxis{ii};
        offset = GPR.Geometry.offset{ii}(jj);
        [nsamp, ntrcs] = size(Rad);
        
        %------------------------------------------------------------------
        % Exponential Time Dependant Gain
        
        % Parameters
        tpow = 2; % Filter Order, 1 is Linear
        % Rolloff assumes antennas h meters off of surface ~ roughly
        h = .5;%1;
        va = 0.3;
        tmpt0 = (2.*h./va); tmpt = sqrt(tmpt0.^2+(offset./va).^2);
        rollOff = round((tmpt./dt));
        
        % Function
%         Rad = gain(Rad, tpow, rollOff );
        % Top  and Bottom Mute
        taper = rollOff+25;
        tukeywindow = tukeywin(taper,1);
        window = [tukeywindow;zeros(round(nsamp-length(tukeywindow)),1)];
        mute = window*ones(1,ntrcs);
        Rad = Rad.*mute;
        
        %------------------------------------------------------------------
        % First Break Picking
        
        % Parameters
        merR = 0.6;  %Rank of MER window [ns]
        powMER = 3; % Power of MER operation
        ix1 = 1;%round((t0./dt)+1); % Because Padding - removed this 11/6/23
        R = 1;%51; % Smoothing Operator length;
        CommonFBPick = zeros(1,ntrcs);
        for kk = 1:length(Rad(1,:))
            % First Break Picking
            tmp = Rad(ix1:end,kk);
            %         tmp = AGCgain(Rad2(:,kk),50,2);
            [CommonFBix, ~, ~] = wong_mer(tmp, merR, dt, powMER);
            CommonFBix = CommonFBix+ix1;
            CommonFBPick(kk) = TWT(CommonFBix);
        end
        % Filter Picks
%         q = quantile(CommonFBPick,q);
%         badix = find(CommonFBPick<q(1) | CommonFBPick>q(2));
%         goodix = find(~ismember(1:ntrcs,badix));
%         filtPicks = interp1(goodix,CommonFBPick(goodix),1:ntrcs,'linear','extrap');
%         filtPicks = movmean(filtPicks,R);
        % LocalQuantile Filter
        X = GPR.Geolocation.X{ii};
        Y = GPR.Geolocation.Y{ii};
        Distance = GPR.Geolocation.Distance{ii};
        r = 1;%2.5;%;%1;%50;%.5; % Bin Radius (m)
        q = [0.05,0.95]; % Quantile Threshold
        [filtPicks, badix] = lqfilter(CommonFBPick,X,Y,r,q);
        % Smooth Picks
        smoothr = 1.5;%%5;%7.5;%1;;%25;%0.5;%2;
        [ filtPicks ] = nonParametricSmooth( Distance,filtPicks,Distance,smoothr );
        % Quik Look
        if isQuickLook
%             figure();imagesc(GPR.Geolocation.Distance{ii},TWT,Rad);colormap(bone);
            figure();pcolor(GPR.Geolocation.Distance{ii},TWT,Rad); axis ij;
            colormap(bone); shading interp;
            hold on;
            plot(GPR.Geolocation.Distance{ii},filtPicks,'.m');
            plot(GPR.Geolocation.Distance{ii}(badix),CommonFBPick(badix),'.b')
%             plot(filtPicks,'.m');
%             plot(badix,CommonFBPick(badix),'.b')
            title(['Offset ',num2str(offset)])
            ylabel('Travel-time (ns)')
            xlabel('Distance (m)')
        end
        % Store Surface Picks
        surfaceTWT(jj,:) = filtPicks;
    end
    GPR.D.surfaceTWT{ii} = surfaceTWT;
%     GPR.D.surfaceH{ii}




end

