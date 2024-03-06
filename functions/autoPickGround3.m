function [GPR] = autoPickGround3(GPR,isQuickLook)
% Check inputs
if nargin == 1
    isQuickLook = 0;
end

for ii = 1:GPR.MD.nFiles
    % Allocation
    offsets = (GPR.Geometry.offset{ii});
    groundTWT = zeros(numel(offsets),size(GPR.D.Radar{ii}{1},2));
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
        % Spreading and Exponential Compensation Gain
        
        % Parameters
        exppow = 0.0001;
        tpow = 1;%10; % Filter Order, 1 is Linear
        
        % Function
%         Rad = SECgain(Rad, exppow, tpow);
        
        %------------------------------------------------------------------
        % First Break Picking
        
        % Parameters
        merR = 1.5;  %Rank of MER window [ns]
        powMER = 3; % Power of MER operation
        ix1 = round((t0./dt)+1); % Because Padding
        CommonFBPick = zeros(1,ntrcs);
%         for kk = 1:length(Rad(1,:))
        parfor kk = 1:length(Rad(1,:))
            % First Break Picking
            tmp = Rad(ix1:end,kk);
            %         tmp = AGCgain(Rad2(:,kk),50,2);
            [CommonFBix, ~, ~] = wong_mer(tmp, merR, dt, powMER);
            CommonFBix = CommonFBix+(ix1-1);
            CommonFBPick(kk) = TWT(CommonFBix);
        end
        
        % LocalQuantile Filter
        X = GPR.Geolocation.X{ii};
        Y = GPR.Geolocation.Y{ii};
        Distance = GPR.Geolocation.Distance{ii};
        r = 5;%.5; % Bin Radius (m)
        q = [0.05,0.95]; % Quantile Threshold
        [filtPicks, badix] = lqfilter(CommonFBPick,X,Y,r,q);
        % Smooth Picks
        smoothr = 2.5;%0.5;%2;
        [ filtPicks ] = nonParametricSmooth( Distance,filtPicks,Distance,smoothr );
        
        % Quik Look
        if isQuickLook
%             if jj == 1
%             figure();imagesc(GPR.Geolocation.Distance{ii},TWT,Rad);
            figure();pcolor(GPR.Geolocation.Distance{ii},TWT,Rad); axis ij;
            colormap(bone); shading interp;
            hold on;
            plot(GPR.Geolocation.Distance{ii},filtPicks,'.m');
            plot(GPR.Geolocation.Distance{ii}(badix),CommonFBPick(badix),'.b')
            title(['Offset ',num2str(offset)])
            ylabel('Travel-time (ns)')
            xlabel('Distance (m)')
%             end
        end
        % Store Surface Picks
        groundTWT(jj,:) = filtPicks;
    end
        GPR.D.groundTWT{ii} = groundTWT;
end
end

