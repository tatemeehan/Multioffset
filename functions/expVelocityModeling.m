function [GPR] = expVelocityModeling(GPR)
% expVelocityModeling.m fits a stacking velocity model of the form
% v = b1 + b2 .* exp(-b3.*t) to the velocity spectra, by minimizing for the
% model parameters (b1, b2, & b3) via coherence weighted nonlinear
% least-squares.

% Exponetial Modeling using nonlinear least-squares
warning('off','all');
parfevalOnAll(@warning,0,'off','all');
vOut = [];vOutIx = [];
for ii = 1 : GPR.MD.nFiles
    % Extract Variables
    tau = GPR.D.stackingTimeAxis{ii}; % Stacking Time Axis
    v = GPR.D.velocityAxis{ii};       % Velocity Axis
    C = GPR.D.velocityCoherence{ii};  % Velocity Spectra
    Cdistance = GPR.Geolocation.cmpDistance{ii}; % MidPoint Distance
    Z = GPR.Geolocation.cmpZ{ii};
    dt = GPR.D.dt{ii};
    timeAxis = GPR.D.TimeAxis{ii};
    R = round(length(timeAxis)./length(tau)); % Coherence Comp. TWT Interval
    % Allocate Model Space
    vTrue = zeros(length(tau),length(C));
    %     vStack = vTrue;
    %     Zed = vTrue;
    %     vRMSo = vTrue;
    flagIx = zeros(length(C),1); badModIx = flagIx; badVRMSix = flagIx;
    count = 0;
    modCount = 0;
    vrmscount = 0;
%     initBeta = zeros(numel(Cdistance),1);
%     parfor jj = 1:length(C)
% %     for jj = 1:length(C)
%         warning('') % Clear last warning message
%         % Curve Fitting
%         tmpC = C{jj};
%         % Taper Mute
%         muteIx = round((185./dt)./R); % 185 NanoSecond Mute
%         muteWin = tukeywin(20); muteWin = muteWin(11:end);
%         % Apply Bottom Mute to Coherence
%         tmpC(muteIx:end,:) = 0;
%         tmpC(muteIx-numel(muteWin):muteIx-1,:) = tmpC(muteIx-numel(muteWin):muteIx-1,:).*muteWin;
%         Cix = find(tmpC>0.4);
%         %         if isempty(Cix)
%         %             count = count +1;
%         %             flagIx(jj) = 1;
%         %         else
%         [tauIx,vIx] = ind2sub(size(tmpC),Cix);
%         % Find Maximum Coherence
%         [~,maxCix] = max(tmpC(Cix));
%         % rmvTauIx = find(tauIx>50);
%         rmvTauIx = find(tauIx > tauIx(maxCix)-10 | tauIx > tauIx(maxCix)+10 );
%         Cix(rmvTauIx) = [];
%         tauIx(rmvTauIx) = [];
%         vIx(rmvTauIx) = [];
%         % Coherence Weights
%         cWeight = tmpC(Cix);
%         % Convert tau and v into a table, which is the form fitnlm() likes the input data to be in.
%         dFit = table(tau(tauIx)',v(vIx)');
%         modelfun = @(b,x) b(1) + b(2) .* exp(-b(3).*x(:, 1));
%         tmp = 0.135:0.001:0.195;
%         %             meanDistW = zeros(length(tmp),1);
%         meanC= zeros(length(tmp),1);
% 
%         for ff = 1:length(tmp)
%             % Initial Model Guess
%             %             beta0 = [0.180,0.05,0.005];
%             beta0 = [tmp(ff),0.05,0.005];
%             mod0 = modelfun(beta0,tau(:));
%             %             % Calculate Distance Weights
%             %             distW = zeros(size(cWeight));
%             %             for kk = 1:length(dFit.Var1)
%             %                 tmpdist = sqrt((dFit.Var1(kk)-tau(:)).^2);
%             %                 [~,minIx] = min(tmpdist);
%             %                 tmpdist = sqrt((dFit.Var2(kk)-mod0(minIx)).^2)./cWeight(kk).^2;
%             %                 [minDist,~] = min(tmpdist);
%             %                 distW(kk) = minDist;
%             %             end
%             %             meanDistW(ff) = mean(distW);
%             %             end
%             % Calculate Coherence Weights
% %             distW = zeros(size(cWeight));
%             Cbin = [];
%             for kk = 1:length(dFit.Var1)
%                 tmpdist = sqrt((dFit.Var1(kk)-tau(:)).^2);
%                 [~,minIx] = min(tmpdist);
%                 tmpdist = sqrt((dFit.Var2(kk)-mod0(minIx)).^2);
%                 biffIx = find(tmpdist<0.001);
%                 if isempty(biffIx)
%                 else
%                     Cbin = [Cbin,cWeight(kk)];
%                 end
%                 %                 [minDist,~] = min(tmpdist);
%                 %                 distW(kk) = minDist;
%             end
%             if isempty(Cbin)
%                 meanC(ff) = 0;
%             else
%                 meanC(ff) = mean(Cbin);
%             end
%         end
%         [~,initIx] = max(meanC);
%         initBeta(jj) = tmp(initIx);
%     end
%     % Convolution Smoothing
%     nc = length(initBeta);
%     R=125;
%     kernel = hamming(2.*R +1);
%     % Pad Data with End Traces prior to Convolution
%     padl = ones(1,R).*mean(initBeta(1:R));
%     padr = ones(1,R).*mean(initBeta(nc-2.*R+1:nc));
%     convdata = [padl,initBeta',padr];
%     initBeta = conv(convdata,kernel'./sum(kernel)','valid')';
%     %     initBeta = movmean(initBeta,250);
    tmpFun = ((0.18-.15)./(1850-1100)).*Z;
    tmpFun0 = 0.15-0.0440;%min(tmpFun);
    initBeta = tmpFun+tmpFun0;
    parfor (jj = 1:length(C),GPR.MD.nWorkers)
%     for jj = 1:length(C)
        modelfun = @(b,x) b(1) + b(2) .* exp(-b(3).*x(:, 1));
        tmpC = C{jj};
        Cix = find(tmpC>0.4);
        [tauIx,vIx] = ind2sub(size(tmpC),Cix);
        % Coherence Weights
        cWeight = tmpC(Cix);
        % Convert tau and v into a table, which is the form fitnlm() likes the input data to be in.
        dFit = table(tau(tauIx)',v(vIx)');
        beta0 = [initBeta(jj),0.05,0.005];
        mod0 = modelfun(beta0,tau(:));
        % Calculate Distance Weights
        distW = zeros(size(cWeight));
        for kk = 1:length(dFit.Var1)
            %                 tmpdist = sqrt((dFit.Var1(kk)-tau(:)).^2);
            %                 [~,minIx] = min(tmpdist);
            tmpdist = sqrt((dFit.Var1(kk)-tau(:)).^2+(dFit.Var2(kk)-mod0(:)).^2);
            [minDist,~] = min(tmpdist);
            distW(kk) = minDist;
        end
        %             meanDistW(ff) = mean(distW);
        %             end
        radius = 0.01;
        rmvIx = find(distW > radius);
        % Remove Outlier Coherence Bubbles
        Cix(rmvIx) = [];
        numC = 5;
        if numel(Cix)<numC
            iter = 0;
            cthresh = 0.4;
            while numel(Cix) < numC
                radius = radius + 0.0025;
                rmvIx = find(distW > radius);
                Cix = find(tmpC>cthresh);
                Cix(rmvIx) = [];
                iter = iter +1;
                if iter > 3
                    radius = 0.1;
                    cthresh = cthresh-0.025;
                    iter = 0;
                    if cthresh < 0.1
                        keyboard
                    end
                end

            end
        end
        [tauIx,vIx] = ind2sub(size(C{jj}),Cix);
        cWeight = tmpC(Cix);
        dFit = table(tau(tauIx)',v(vIx)');
        % Fit Exponential Model
        %             modelfun = @(b,x) tmp(initIx) + b(2) .* exp(-b(3).*x(:, 1));
        modelfun = @(b,x) b(1) + b(2) .* exp(-b(3).*x(:, 1));
        beta0 = [initBeta(jj),0.05,0.005];
        mdl = fitnlm(dFit, modelfun, beta0,'Weights',cWeight);
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            modCount = modCount+1;
            badModIx(jj) = 1;
        end
        % Forward Model
        vRMS = zeros(length(tau),1);
        for kk = 1:length(tau)
            vRMS(kk,1) = modelfun(mdl.Coefficients.Estimate,tau(kk));
        end
        % if VRMS is not generally decreasing with depth mark as bad
        if vRMS(25) < vRMS(end-25)
            vrmscount = vrmscount+1;
            badVRMSix(jj) = 1;
        end
        %             vRMSo(:,jj) = vRMS;
        %             Zed(:,jj) = (vRMS(:).*tau(:))./2;
        % Interval Velocity
        % Depth
        % [vInt] = intervalVelocity(vRMS,Zed(:,jj));
        % Time
        [vInt] = intervalVelocity(vRMS,tau(:));
        % Fib this inverval such that it is bounded
        % vInterval = rescale(vInterval,0.17,vRMS(1));
        fibIx = find(vInt<0.169,1);
        vInt(fibIx:end)=0.169;
        vTrue(:,jj) = vInt;
        % Re-stack Stacking Velocity
        % Time
        %             vStack(:,jj) = stackingVelocity(vInt,tau(:));
        % Depth
        % vStack(:,jj) = stackingVelocity(vInt,Zed(:,jj));
        %         end
    end
    %% Interpolate to CMP Gather Domain
    % Filter Interval Velocity Model
    rmvIx = find(badVRMSix);
    if ~isempty(rmvIx)
        vTrueF = vTrue;
        CdistanceF = Cdistance;
        rmvIx = unique(rmvIx);
        vTrueF(:,rmvIx) = [];
        CdistanceF(:,rmvIx) = [];
        Vtrue = medfilt2(vTrueF,[25,250],'symmetric');
        % Interpolate b/w removed traces
        VtrueOut = zeros(size(vTrue));
        for kk = 1:numel(tau)
            VtrueOut(kk,:) = interp1(CdistanceF,Vtrue(kk,:),Cdistance,'linear','extrap');
        end
        VtrueOut = medfilt2(VtrueOut,[25,250],'symmetric');
    else
        % Median Filter Interval Velocity
        vTrue = medfilt2(vTrue,[25,250],'symmetric');
        VtrueOut = vTrue;
    end
    % Stitch all Models together
    vOut = [vOut,VtrueOut];
    vOutIx = [vOutIx,size(VtrueOut,2)];
end
vOutIx = cumsum(vOutIx);
%     nc = numel(Cdistance);
%     nr = numel(tau);
[nr,nc] = size(vOut);
% Convolution Smoothing
R=100;
kernel = hamming(2.*R +1);
% Pad Data with End Traces prior to Convolution
padl = ones(nr,R).*vOut(:,1);
padr = ones(nr,R).*mean(vOut(:,nc-2.*R+1:nc),2);
convdata = [padl,vOut,padr];
% Row Wise Convolution
padu = ones(R,nc+2.*R).*(convdata(1,:));
padd = ones(R,nc+2.*R).*convdata(nr,:);
convdata = [padu;convdata;padd];
Vsmooth = conv2(kernel./sum(kernel),kernel./sum(kernel),convdata(:,:),'valid');
clear('convdata');
%     % Convolution Smoothing
%     R=100;
%     kernel = hamming(2.*R +1);
%     % Pad Data with End Traces prior to Convolution
%     padl = ones(nr,R).*VtrueOut(:,1);
%     padr = ones(nr,R).*mean(VtrueOut(:,nc-2.*R+1:nc),2);
%     convdata = [padl,VtrueOut,padr];
%     % Row Wise Convolution
%     padu = ones(R,nc+2.*R).*(convdata(1,:));
%     padd = ones(R,nc+2.*R).*convdata(nr,:);
%     convdata = [padu;convdata;padd];
%     Vtrue = conv2(kernel./sum(kernel),kernel./sum(kernel),convdata(:,:),'valid');
%     clear('convdata');
zzTop = [];zzTopIx = [];
for ii = 1 : GPR.MD.nFiles
    tau = GPR.D.stackingTimeAxis{ii}; % Stacking Time Axis
    Cdistance = GPR.Geolocation.cmpDistance{ii}; % MidPoint Distance
    % Interval Velocity
    % Re-Segment Velocity Model
    if ii == 1
        Vtrue = Vsmooth(:,1:vOutIx(ii));
    else
        Vtrue = Vsmooth(:,vOutIx(ii-1)+1:vOutIx(ii));
    end
    % Interpolate in Distance
    nt = numel(tau); nx = numel(GPR.Geolocation.Distance{ii});
    tmpVtrue = Vtrue; Vtrue = zeros(nt,nx);
    % Final Stacking Velocity Model vs tau and GPR Distance
    for kk = 1:nt
        Vtrue(kk,:) = interp1(Cdistance,tmpVtrue(kk,:),GPR.Geolocation.Distance{ii},'linear','extrap');
    end
    % Interpolate to GPR CMP Gather Travel Time Axis
    Tau = GPR.D.TimeAxis{ii};
    nT = numel(Tau);
    tmpVtrue = Vtrue; Vtrue = zeros(nT,numel(GPR.Geolocation.Distance{ii}));
    for kk = 1:nx
        Vtrue(:,kk) = interp1(tau(:),tmpVtrue(:,kk),Tau(:),'linear','extrap');
    end
    % % Now Convert to Stacking Velocity and Depth and Density
    Vstack = zeros(size(Vtrue));
    for kk = 1:nx
        Vstack(:,kk) = stackingVelocity(Vtrue(:,kk),Tau(:));
    end
    % Calculate Depth
    Zstack = (Vstack.*Tau)./2;
    zzTop = [zzTop,Zstack];
    zzTopIx = [zzTopIx,size(Zstack,2)];
    % Axis for Depth Image
    DepthAxis = (linspace(min(Zstack(:)),max(Zstack(:)),nT))';

    %     % Calculate Interval Velocity in Depth
    %     Vz = zeros(size(Vtrue));
    %     for kk = 1:nx
    %         Vz(:,kk) = intervalVelocity(Vstack(:,kk),Zstack(:,kk));
    %     end

    % Calculate Density
    FirnDensity = DryCrim(Vtrue);
    LWC = 0.005;
    wetFirnDensity = WetCrim(Vtrue,LWC);

    % Store Output in GPR Structure
    GPR.D.stackingVelocity{ii} = Vstack;
    GPR.D.intervalVelocity{ii} = Vtrue;
    GPR.D.Depth{ii} = Zstack;
    GPR.D.DepthAxis{ii} = DepthAxis;
    GPR.D.Density{ii} = FirnDensity;
    GPR.D.wDensity{ii} = wetFirnDensity;
    GPR.D.LWC{ii} = LWC;
end
zzTopIx = cumsum(zzTopIx);
% Axis for Depth Image & Models (Same for all Lines of Batch)
DepthAxis = (linspace(min(zzTop(:)),max(zzTop(:)),nT))';
% Velocity Models are in Time; Density Models are in Depth
for ii = 1 : GPR.MD.nFiles
    nx = numel(GPR.Geolocation.Distance{ii});
    %     % Re-Segment Depth Model
    %     if ii == 1
    %         Zstack = zzTop(:,1:zzTopIx(ii));
    %     else
    %         Zstack = zzTop(:,zzTopIx(ii-1)+1:zzTopIx(ii));
    %     end
    % Re-Segment Depth Model
    Zstack = GPR.D.Depth{ii};

    GPR.D.DepthAxis{ii} = DepthAxis;

    % Convert Firn Density to Depth
    FirnDensity = GPR.D.Density{ii};
    wetFirnDensity = GPR.D.wDensity{ii};
    FirnDensityZ = zeros(size(FirnDensity));
    FirnDensityWZ = zeros(size(FirnDensity));

    for kk = 1:nx
        FirnDensityZ(:,kk) = interp1(Zstack(:,kk),FirnDensity(:,kk),DepthAxis(:),'linear','extrap');
        FirnDensityWZ(:,kk) = interp1(Zstack(:,kk),wetFirnDensity(:,kk),DepthAxis(:),'linear','extrap');
    end
    GPR.D.Density{ii} = FirnDensityZ;
    GPR.D.wDensity{ii} = FirnDensityWZ;
end
% End of Function
end