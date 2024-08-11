function [GPR] = NormalMoveoutCorrection(GPR,isNMOstacking)
%% Normal Moveout Correction
% A step of 2D interpolation of the NMO velocity model is computed during
% an initial iteration over the first channel. Subsequent channels are
% processed in parallel.
%
% Output Variables where mxn is #samples x #traces
% RadarNMO ~ the NMO corrected common-offset gathers

if nargin < 2
    isNMOstacking = 1;
end

% Maximum Allowable Wavelet Stretch [%]
Stretch = 100;
for ii = 1:GPR.MD.nFiles
    % Allocation for Parfor Overhead
    dumX = GPR.Geolocation.Distance{ii};
    dumT = GPR.D.TimeAxis{ii};
    dumV = GPR.D.stackingVelocity{ii};
    offsetArray = GPR.Geometry.offset{ii};
    dt = GPR.D.dt{ii};
    f0 = GPR.D.f0{ii};
    isInterpolate = 0;
    Radar = GPR.D.Radar{ii};
    parfor (jj = 1:GPR.Geometry.nChan{ii}, GPR.MD.nWorkers)
        %         for jj = 1:GPR.Geometry.nChan{ii}
        % Perform NMO Correction
        [RadarNMO{jj},~,~,~,~] = ...
            commonOffsetNMO(Radar{jj},dt,f0,offsetArray(jj),dumX,dumT,dumV,Stretch,isInterpolate);
    end
    GPR.D.RadarNMO{ii} = RadarNMO;
    if isNMOstacking
        nearChans = [1:3];
        nearStack = mean(cat(3,RadarNMO{nearChans}),3);
        farChans = [13:16];
        farStack = mean(cat(3,RadarNMO{farChans}),3);
        % Weighted Sum of Near and Far Stack
        twtOverlap = 225; % [ns]
        %             % Correct Mean Bias
        % %             nearMean = mean(mean(nearStack(twtOverlap./dt - 50:twtOverlap./dt+50,:)));
        % %             farMean = mean(mean(farStack(twtOverlap./dt - 50:twtOverlap./dt+50,:)));
        %             nearMean = (mean(nearStack(twtOverlap./dt - 50:twtOverlap./dt+50,:)));
        %             farMean = (mean(farStack(twtOverlap./dt - 50:twtOverlap./dt+50,:)));
        % %             nearMean = mean(nearStack(:));
        % %             farMean = mean(farStack(:));
        %             cormean = movmean(farMean-nearMean,251);
        %             farStack = farStack-cormean;
        % Design Blend Window
        win = triang(200);
        nearW = zeros(length(dumT),1);
        nearWin = [ones(twtOverlap./dt,1);win(101:end)];
        nearW(1:numel(nearWin)) = nearWin;
        farW = ones(length(dumT),1);
        farWin = [zeros(twtOverlap./dt,1);win(1:100)];
        farW(1:numel(farWin)) = farWin;
        % Stack Near and Far Offsets
        nearStack = nearStack.*nearW;
        farStack = farStack.*farW;
        stack = zeros(size(nearStack,1),size(nearStack,2),2);
        stack(:,:,1) = nearStack; stack(:,:,2) = farStack;

        % Offset Stacking
        GPR.D.RadarStack{ii} = sum(stack,3); % Windowed
        %            GPR.D.RadarStack{ii} = mean(cat(3,RadarNMO{stackChans}),3);

        % Post-Stack Processing

        % Top Mute
        ns = 40./dt; % 40 ns mute tapered window
        win = hamming(2.*ns+1);
        win = [win(1:ns);ones(length(dumT)-(ns),1)];
        GPR.D.RadarStack{ii} = GPR.D.RadarStack{ii}.*win;

        % RMS Amplitude
        % RMS Win
        ns = 2./dt; % 2 ns RMS Window;
        GPR.D.RadarStack{ii} = rmsAmplitude(GPR.D.RadarStack{ii},ns);

        % Wiener Filter
        %             GPR.D.RadarStack{ii} = wiener2(GPR.D.RadarStack{ii},[3,3]);
        % Trace Stacking
        R = 25; % Number of Traces in Weighted Stack 2R+1
        GPR.D.RadarStack{ii} = StakR(GPR.D.RadarStack{ii},R);

        % Anisotropic Diffusion Filter
        opts.T = 15; % Diffusion Time Steps
        opts.dt = 2; % Time Step Interval
        opts.Scheme = 'R'; % Rotational Invariant Implicit Solver
        opts.eigenmode = 2; % Edge-Enhancment Diffustion Tensor
        opts.lambda_e = 0.25; % EED Planar Structure Contrast
        opts.sigma = 3; % Image Kernal Smoothing
        opts.rho = 3; % Hessian Kernal Smoothing
        GPR.D.RadarStack{ii} = CoherenceFilter(GPR.D.RadarStack{ii},opts);

        % AGC Gain
        %             GPR.D.RadarStack{ii} = AGCgain(GPR.D.RadarStack{ii},100,2);

        % Ampltudes are Lost at the Beginnings and Ends of Files Due to
        % Data Fold of CMP Stacking.. Let's Recover them!
        % Recover Amplitudes at Start of Image
        if ii == 1
            nTrcs = 50; % Window of Traces to Gain
            % Extract Data Gather 1
            rad1 = GPR.D.RadarStack{ii}(:,1:nTrcs);
            % Automatically Determine Bottom of Image (lol)
            ix = wong_mer(flipud(mean(rad1,2)),10,dt,3);
            ns = size(rad1,1)-ix;% TimeWindow Threshold of Image Bottom;
            % Amplitude Recovery Algorithm
            rad = zeros(size(rad1)); % Allocation
            % Average Amplitude of Traces Away from Boundary
            meanTrc = mean(GPR.D.RadarStack{ii}(:,nTrcs+1:2.*nTrcs),2);
            % Amplitude Bias of Each Trace
            rad(1:ns,:) = (meanTrc(1:ns) - (rad1(1:ns,:)));
            % RMS Amplitude
            if any(sign(meanTrc(:)) == -1)
                rad = rmsAmplitude(rad,10);
            end
            % Smooth Amplitude Biases
            amps = conv2(rad,hamming(250)./sum(hamming(250)),'same');
            amps = conv2(amps,hamming(25)'./sum(hamming(25)),'same');
            % Preserve Sign of Signal
            amps = amps.*(sign(meanTrc(:))*ones(1,nTrcs));
            amps(ns:end,:) = meanTrc(ns:end)*ones(1,nTrcs);
            % Apply Amplitude Gain
            radout = rad1+amps;
            % Pad Image Boundary For Smoothing
            radout2 = [radout,GPR.D.RadarStack{ii}(:,nTrcs+1:nTrcs+51)];
            % Smooth Across Boundary
            radout2 = StakR(radout2,25,'kernel');
            % Remove Padded Traces
            radout2(:,nTrcs+1:nTrcs+51) = [];
            % Stitch Recovered Amplitudes into Stacked Gathers
            GPR.D.RadarStack{ii}(:,1:nTrcs) = radout2;
        end

        % Recover Amplitudes Across File Boundaries
        if ii > 1
            nTrcs = 50; % Window of Traces to Gain
            % Extract Data Gather 1
            rad1 = GPR.D.RadarStack{ii-1}(:,end-nTrcs+1:end);
            % Extract Data Gather 2
            rad2 = GPR.D.RadarStack{ii}(:,1:nTrcs);
            % Concatenate Image Files
            rad3 = [rad1,rad2];
            % Automatically Determine Bottom of Image (lol)
            ix = wong_mer(flipud(mean(rad3,2)),10,dt,3);
            ns = size(rad3,1)-ix;% TimeWindow Threshold of Image Bottom;
            % Amplitude Recovery Algorithm
            rad = zeros(size(rad3)); % Allocation
            % Average Amplitude of Traces Away from Boundary
            meanTrc1 = mean(GPR.D.RadarStack{ii-1}(:,end-2.*nTrcs+1:end-nTrcs),2);
            meanTrc2 = mean(GPR.D.RadarStack{ii}(:,nTrcs+1:2.*nTrcs),2);
            meanTrc = mean([meanTrc1(1:ns),meanTrc2(1:ns)],2);
            % Amplitude Bias of Each Trace
            rad(1:ns,:) = (meanTrc - (rad3(1:ns,:)));
            % RMS Amplitude
            if any(sign(meanTrc(:)) == -1)
                rad = rmsAmplitude(rad,10);
            end
            % Smooth Amplitude Biases
            amps = conv2(rad,hamming(250)./sum(hamming(250)),'same');
            amps = conv2(amps,hamming(25)'./sum(hamming(25)),'same');
            % Preserve Sign of Signal
            amps = amps.*(sign(mean([meanTrc1(:),meanTrc2(:)],2))*ones(1,2.*nTrcs));
            amps(ns:end,:) = mean([meanTrc1(ns:end),meanTrc2(ns:end)],2)*ones(1,2.*nTrcs);
            % Apply Amplitude Gain
            radout = rad3+amps;
            % Pad Image Boundary For Smoothing
            radout2 = [GPR.D.RadarStack{ii-1}(:,end-nTrcs+1-50:end-nTrcs+1),...
                radout,GPR.D.RadarStack{ii}(:,nTrcs+1:nTrcs+51)];
            % Smooth Across Boundary
            radout2 = StakR(radout2,25,'kernel');
            % Remove Padded Traces
            radout2(:,1:51) = []; radout2(:,2.*nTrcs+1:end) = [];
            % Stitch Recovered Amplitudes into Stacked Gathers
            GPR.D.RadarStack{ii}(:,1:nTrcs) = radout2(:,nTrcs+1:end);
            GPR.D.RadarStack{ii-1}(:,end-nTrcs+1:end) = radout2(:,1:nTrcs);
        end

        % Recover Amplitudes at End of Image
        if ii == GPR.MD.nFiles
            nTrcs = 50; % Window of Traces to Gain
            % Extract Data Gather 1
            rad1 = GPR.D.RadarStack{ii}(:,end-nTrcs+1:end);
            % Automatically Determine Bottom of Image (lol)
            ix = wong_mer(flipud(mean(rad1,2)),10,dt,3);
            ns = size(rad1,1)-ix;% TimeWindow Threshold of Image Bottom;
            % Amplitude Recovery Algorithm
            rad = zeros(size(rad1)); % Allocation
            % Average Amplitude of Traces Away from Boundary
            meanTrc = mean(GPR.D.RadarStack{ii}(:,end-2.*nTrcs+1:end-nTrcs),2);
            % Amplitude Bias of Each Trace
            rad(1:ns,:) = (meanTrc(1:ns) - (rad1(1:ns,:)));
            % RMS Amplitude
            if any(sign(meanTrc(:)) == -1)
                rad = rmsAmplitude(rad,10);
            end
            % Smooth Amplitude Biases
            amps = conv2(rad,hamming(250)./sum(hamming(250)),'same');
            amps = conv2(amps,hamming(25)'./sum(hamming(25)),'same');
            % Preserve Sign of Signal
            amps = amps.*(sign(meanTrc(:))*ones(1,nTrcs));
            amps(ns:end,:) = meanTrc(ns:end)*ones(1,nTrcs);
            % Apply Amplitude Gain
            radout = rad1+amps;
            % Pad Image Boundary For Smoothing
            radout2 = [GPR.D.RadarStack{ii}(:,end-nTrcs+1-50:end-nTrcs+1),radout];
            % Smooth Across Boundary
            radout2 = StakR(radout2,25,'kernel');
            % Remove Padded Traces
            radout2(:,1:51) = [];
            % Stitch Recovered Amplitudes into Stacked Gathers
            GPR.D.RadarStack{ii}(:,end-nTrcs+1:end) = radout2(:,1:nTrcs);
        end
    end
end

if isNMOstacking
    % Concatenate Profile
    radarCat = []; radarCatIx = [];
    % RMS Amplitude Sign Recovery
    for ii = 1:GPR.MD.nFiles
        % Detrend Amplitudes
        if any(sign(meanTrc(:)) == -1)
        else
            nTrcs = size(GPR.D.RadarStack{ii},2);
            time = GPR.D.TimeAxis{ii};
            tmpRadar = GPR.D.RadarStack{ii};
            bpns = 55; % Break Point [ns]
            %             parfor (jj = 1:nTrcs, GPR.MD.nWorkers) % parfor is slower
            for jj = 1:nTrcs
                % Fit Piece Wise Model with Breakpoint at 55 ns;
                tmp = tmpRadar(:,jj);
                dum = zeros(length(tmp),1);
                breakIx = fix(bpns./dt);
                % Piece 1
                G = [ones(breakIx,1),time(1:breakIx)];
                d = tmp(1:breakIx);
                m = G\d;
                dcal1 = G*m;
                % Piece 2
                G = [ones(ns-breakIx,1),time((breakIx+1):ns)];
                d = tmp((breakIx+1):ns);
                m = G\d;
                dcal2 = G*m;
                dcal = [dcal1;dcal2];
                dcal = movmean(dcal,25);
                % Detrend Amplitudes
                dum(1:ns) = tmp(1:ns) - dcal;
                dum(ns-1:end) = movmean(dum(ns-1:end),25);
                tmpRadar(:,jj) = dum;
                % Fit Linear Model to RMS Amplitudes
                %                 G = [ones(ns,1),time(1:ns)];
                %                 d = tmpRadar(1:ns,jj);
                %                 m = G\d;
                %                 dcal = G*m;
                %                 tmpRadar(:,jj) = tmp - dcal;
            end
            GPR.D.RadarStack{ii} = tmpRadar;
        end
        radarCat = [radarCat,GPR.D.RadarStack{ii}];
        radarCatIx = [radarCatIx,size(GPR.D.RadarStack{ii},2)];
    end
    radarCatIx = cumsum(radarCatIx);

    % Processing on Entire Image
            % Spiking Deconvolution
%             tmpRadar = GPR.D.RadarStack{ii};
            NF = 250;     % Operator Length
            mu = 100;     % Pre-whitening
            [~,radarCat] = spiking(radarCat,NF,mu);
%             GPR.D.RadarStack{ii} = tmpRadar;
    
            % Anisotropic Diffusion Filter
            opts.T = 15; % Diffusion Time Steps
            opts.dt = 2; % Time Step Interval
            opts.Scheme = 'R'; % Rotational Invariant Implicit Solver
            opts.eigenmode = 2; % Edge-Enhancment Diffustion Tensor
            opts.lambda_e = 0.25; % EED Planar Structure Contrast
            opts.sigma = 1; % Image Kernal Smoothing
            opts.rho = 1; % Hessian Kernal Smoothing
            GPR.D.RadarStack{ii} = CoherenceFilter(GPR.D.RadarStack{ii},opts);
            radarCat = CoherenceFilter(radarCat,opts);
    for ii = 1:GPR.MD.nFiles
        % Re-Segment Imagery
        if ii == 1
            GPR.D.RadarStack{ii} = radarCat(:,1:radarCatIx(ii));
        else
            GPR.D.RadarStack{ii} = radarCat(:,radarCatIx(ii-1)+1:radarCatIx(ii));
        end
    end
end
% End of Function
end