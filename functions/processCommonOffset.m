function [ Rad, TWT, t0 ] = processCommonOffset( Rad, f0, dt, TWT, offset, x, dx, chanNo )
% processCommonOffset is a data processing Subroutine
%   
%   Written by Tate Meehan, Boise State University, GreenTrACS 2016
%
% Process Data
% Filter Parameters are Established within Function
isDisplay = 0;       % Control Text Display to Screen
isNormalize = 0;     % Control Flag to Normalize Data
isMedFilt = 0;       % Control Flag for Median Subtraction Filter
isBandPass = 1;      % Control Flag to Band-Pass Filter Data
isFK = 0;            % Control Flag for FK Filter
isSpiking = 0;       % Control Flag for Spiking Deconvolution
isDemultiple = 0;    % Control Flag for Predictive Deconvolution
isCalibration = 1;   % Control Flag for Air-wave Calibration using Priors
isTimeZero = 0;      % Control Flag for Time-Zero Correction
isShearletFilt = 0;  % Control Flag for Shearlet Background Subtraction
isSVDfilt = 0;       % Control Flag for SVD Component Subtraction Filter
isSpatialMedFilt = 1;% Control Flag for Spatial Median Subtraction Filter
isExpGain = 0;       % Control Flag for Ramped Gain of Data
isAGCgain = 1;       % Control Flag for AGC Gain of Data
isSECgain = 0;       % Control Flag for SEC Gain of Data
isStak = 0;          % Control Flag to Stack Data
isKuwahara = 0;      % Control Flag for Kuwahara Filter
isWiener = 1;        % Control Flag for Wiener Filter
isRMSamplitude = 0;  % Control Flag for RMS Amplitdues
isMute = 0;          % Control Flag for Taper Mute

    %----------------------------------------------------------------------      
    % Convert Units
    f0Hz = f0 * 1e6;        % [Hz]
    dtSec = dt * 1e-9;      % [s]
    [nsamp, ntrcs] = size(Rad);

    %----------------------------------------------------------------------
    % Normalize Data
    if isNormalize
        if isDisplay
        display( 'Begin Normalize')
        tic
        end
        
        Rad = normalize( Rad );
        
        if isDisplay
        display( 'Normalize Done')
        toc
        display(' ')
        end
    end    
    
        
%     %----------------------------------------------------------------------
%     % Median Subtraction Filter
%     if isMedFilt
%         if isDisplay
%         display( 'Begin Median Subtraction Filter')
%         tic
%         end
%         
%         % Parameters
%         % Rank of Median Subtraction Filter
% %                 MedFiltR = 10;
%         % Nominal Frequency Pass
%         MedFiltR = 2.*(ceil(1/((f0Hz)*dtSec))-1)+1;
%         % Low Pass
% %          MedFiltR = 2.*(ceil(1/((f0Hz/2.5)*dtSec))-1)+1;
%         % High Pass
% %           MedFiltR = 2.*(ceil(1/((f0Hz*2.5)*dtSec))-1)+1;
% 
%         % Function
%         Rad = medfilt1( Rad, MedFiltR, [], 2,'omitnan','truncate' );
%         
%         if isDisplay
%         display( 'Median Subtraction Filter Done')
%         toc
%         display(' ')
%         end
%     end
    
    %----------------------------------------------------------------------    
    % Band-Pass Filter
    if isBandPass
        if isDisplay
        display( 'Begin Band-Pass Filter')
        tic
        end
        
        % Parameters
        % Build Two Octave Filter About Nominal Frequency
        fMin = f0Hz/2; % [Hz]   
        fMax = f0Hz*2; % [Hz]
%         fMin = f0Hz/2; % [Hz]   
%         fMax = f0Hz; % [Hz]
        
        % Function
        Rad = bpfilter( Rad, dtSec, fMin, fMax, 8 );
        
        if isDisplay
        display( 'Band-Pass Filter Done')
        toc
        display(' ')
        end
    end
    %----------------------------------------------------------------------
    % FK Filter
    if isFK
    [ Rad ] = fk_t8( Rad, dt, dx );
    end
    
    %----------------------------------------------------------------------
    % Spiking Deconvolution
    if isSpiking
        NF = 200;    % Operator Length
        mu = 1e6;     % Pre-whitening
        [~,Rad] = spiking(Rad,NF,mu);

    end
    %----------------------------------------------------------------------
    % Predictive Deconvolution
    if isDemultiple
        NF = 1;    % Operator Length
        L = 10;    % Prediction Length
        mu = 1;    % Pre-whitening
        % Demultiple
        Rad = demultiple(Rad,NF,L,mu);
    end    
    %----------------------------------------------------------------------
    % Apriori Air-wave Calibration
    if isCalibration
        % Airwave Arrival Times are Picked using pickAirwaveCalibration.m
        % This could be better automated
        path = 'E:\JIF24\JIF_MO';
        file = 'AirwaveCalibration.mat';
        [Rad, TWT, t0] = airwaveCalibrationJIF( Rad, TWT, dt, chanNo, offset, path, file);
%         if chanNo == 1
%             disp('Choose an Airwave Calibration File')
%             [file,path] = uigetfile;
%             currentDir = pwd;
%             cd('C:\Users\snowfield\Desktop\');
%             save('calibrationLocation.mat','file','path')
%             cd(currentDir);
%             calibration = load(fullfile(path,file));
%             tmpfieldname = fieldnames(calibration);
%             calibrationIx = getfield(calibration,tmpfieldname{1});
%         else
%             load('C:\Users\snowfield\Desktop\calibrationLocation.mat')
%             calibration = load(fullfile(path,file));
%             tmpfieldname = fieldnames(calibration);
%             calibrationIx = getfield(calibration,tmpfieldname{1});
%         end
%          % Parameters      
%         merR = 1;   %Rank of MER window [ns]
%         powMER = 3; % Power of MER operation
%        
%         % Function
% %         [Rad, TWT, t0] = airwaveCalibration( Rad, TWT, dt, merR, powMER, offset, chanNo, calibrationIx, path );
%         [Rad, TWT, t0] = airwaveCalibrationDux( Rad, TWT, dt, merR, powMER, offset, chanNo, calibrationIx, path );

    end
    %----------------------------------------------------------------------
    % Time-Zero Correction
    if isTimeZero
        if isDisplay
        display( 'Begin Time Zero Correction')
        tic
        end
        
        
        % Parameters      
        merR = 1;   %Rank of MER window [ns]
        powMER = 3; % Power of MER operation
       
        % Function
        [Rad, TWT, t0] = timeZero( Rad, TWT, dt, merR, powMER, offset );        
        [nsamp, ~] = size(Rad);

        if isDisplay
        display( 'Time Zero Correction Done')
        toc
        display(' ')
        end
    elseif ~isCalibration
        t0 = 0;
    end
    %----------------------------------------------------------------------
    % Shearlet Transform Background Subtraction Filter    
    if isShearletFilt
        numScales = 5;
        pow = 5;
        [Rad] = shearletFilter(Rad,t0,dt,numScales,pow);
    end
    %----------------------------------------------------------------------
    % SVD Principle Component Subtraction Filter
    if isSVDfilt
        if isDisplay
        display( 'Begin Singular Value Decomposition Filter')
        tic
        end
        PCA = 0.9; % PCA percentage threshold
        Rad = SVDSfilter( Rad, PCA );
        
        if isDisplay
        display( 'Singular Value Decomposition Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Spatial Median Subtraction
    %     SpatialMedianSubtract
    if isSpatialMedFilt
        if isDisplay
            display('Begin Coherent Noise Subtraction')
            tic
        end
        R = 500;
                [Rad] = backSubtractBigMess(Rad,x,R);
%         [Rad] = backSubtract(Rad,x,R);
%         [ Rad ] = movingMedianSubtraction( Rad, round(size(Rad,2).*.1) );
        
        if isDisplay
            display('Coherent Noise Removal Done')
            toc
            display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Median Subtraction Filter
    if isMedFilt
        if isDisplay
        display( 'Begin Median Subtraction Filter')
        tic
        end
        
        % Parameters
        % Rank of Median Subtraction Filter
%                 MedFiltR = 10;
        % Nominal Frequency Pass
        MedFiltR = 2.*(ceil(1/((f0Hz)*dtSec))-1)+1;
        % Low Pass
%          MedFiltR = 2.*(ceil(1/((f0Hz/2.5)*dtSec))-1)+1;
        % High Pass
%           MedFiltR = 2.*(ceil(1/((f0Hz*2.5)*dtSec))-1)+1;

        % Function
        Rad = medfilt1( Rad, MedFiltR, [], 2,'omitnan','truncate' );
        
        if isDisplay
        display( 'Median Subtraction Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Exponential Time Dependant Gain
    if isExpGain
        if isDisplay
        display( 'Begin Power Gain')
        tic
        end
        
        % Parameters
        tpow = 2; % Filter Order, 1 is Linear
        
        % Function
        Rad = gain(Rad, tpow, rollOff ); 
        
        if isDisplay
        display( 'Power Gain Done')
        toc
        display(' ')
        end
    end
    
    % ---------------------------------------------------------------------
    % Automatic Gain Control
    if isAGCgain
        if isDisplay
            tic
            display( 'Begin AGC')
        end
        
        % Parameters
        R = 100;%ceil(nsamp/2);
        type = 2; % Trace Normalize: 0 = None, 1 = amplitude, 2 = RMSnorm

        Rad = AGCgain(Rad,R,type);
        
        if isDisplay
            display( 'AGC Done')
            toc
            display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Spreading and Exponential Compensation Gain
    if isSECgain
        if isDisplay
            display( 'SEC Gain')
            tic
        end
        
        %         % Parameters
        %         exppow = 0.001; % Exponential Order
        %         tpow = 2; % Power Order, 1 is Linear
        %
        % Parameters
        exppow = -0.005; % Exponential Order
        tpow = 2; % Power Order, 1 is Linear
        
        % Function
        Rad = SECgain(Rad, exppow, tpow);
    end

    %----------------------------------------------------------------------
    % RMS Amplitude
    if isRMSamplitude
        if isDisplay
            display( 'Begin RMS Amplitude')
            tic
        end
        % Parameters
        L = 25;
        % Function
        [Rad] = rmsAmplitude(Rad,L);
        if isDisplay
            display( 'RMS Amplitude Done')
            toc
            display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Edge Preserving Horizon Filter    
    if isKuwahara
        Rad = Kuwahara(Rad,5);
    end

    %----------------------------------------------------------------------
    % Wiener Filter
    if isWiener
        % 5 ns by 50 m Filter Window
        nr = 5./dt; nc = 50./dx;
        Rad = wiener2(Rad,[nr nc]);
    end
    %----------------------------------------------------------------------
    % Stacking Filter
    if isStak
        if isDisplay
        display( 'Begin Stack')
        tic
        end
        
        % Parameters
        StakFiltR = 10; % Filter Rank
        
        % Function
        Rad = StakR(Rad,StakFiltR);
        
        if isDisplay
        display( 'Stacking Done')
        toc
        display(' ')
        end
    end

    %----------------------------------------------------------------------
    % Top  and Bottom Mute
    if isMute
        % Parameters
        rollOff = m - dt.*5;
        R = rollOff./nsamp;%0.25;
        tukeywindow = tukeywin(rollOff.*2,1);
        window = [tukeywindow;zeros(nsamp-2.*rollOff,1)];
        mute = window*ones(1,ntrcs);
        Rad = Rad.*mute;
        % Top  and Bottom Mute
        R = rollOff./nsamp;%0.25;
        L = 101;
        window = tukeywin(L,1); zed = zeros(nsamp,1);
        zed(rollOff-ceil(L/2)+1:rollOff-ceil(L/2)+L) = window;
        window = zed;
        mute = window*ones(1,ntrcs);
        Rad = Rad.*mute;
    end


end

