function [GPR] = predictiveDeconvolution(GPR)
isDepthFXdecon = 1;
isSpatialFXdecon = 0;

if isDepthFXdecon
    mu = 1e7; % Weighting
    nfR = 50; % Length of Deconvolution Filter
    fxProRadar = cell(GPR.MD.nFiles,1); % Allocation

    for ii = 1:GPR.MD.nFiles
        dt = GPR.D.dt{ii};
        f0 = GPR.D.f0{ii};
        % Spatial Predictive Filtering for Random Noise Suression
        Window = 20/dt; % Number of Temporal Samples in Fx-Decon Window Per Computation
        Overlap = round(Window/10); % Number of Overlapping Computations
        RadarStack = GPR.D.RadarStack{ii};
        %         parfor (jj = chan, nWorkers)
        %         for jj = chan
        tmp = 1:length(RadarStack(:,1));
        filterLap = Overlap:Overlap:length(tmp);
        % Design Data Window
        for kk = filterLap
            if kk == filterLap(1)
                run = 1;
                traceIx = 1:kk+(Overlap);
                tmpData = RadarStack(traceIx,:);
            elseif kk == filterLap(end)
                run = length(filterLap);
                traceIx = kk-(Overlap):length(tmp);
                tmpData = RadarStack(traceIx,:);
            else
                run = find(kk==filterLap);
                traceIx = kk-(Overlap):kk+(Overlap);
                tmpData = RadarStack(traceIx,:);
            end
            % Perform Fx-Decon within Windowed Data
            fxProRadar{ii}{run,1} = fx_deconT8(tmpData,dt,nfR,mu,f0);
            fxProRadar{ii}{run,2} = traceIx;
            %             end
        end

        fxRadar = cell(GPR.MD.nFiles,1);   % Allocate Temporary Arrays
        fxCatRadar = cell(GPR.MD.nFiles,1);% Temporary Array

        %         parfor (jj = chan, nWorkers)
        % Concatenate Deconvolved Matricies
        fxCatRadar{ii}{1} = cat(1,fxProRadar{ii}{:,1});
        % Sort Trace Indicies
        [fxCatRadar{ii}{2},fxCatRadar{ii}{3}] = sort([fxProRadar{ii}{:,2}]);
        fxRadar{ii} = fxCatRadar{ii}{1}(fxCatRadar{ii}{3},:);

        % Fold Replicated Traces by Mean Stacking
        for kk = 1:length(RadarStack(:,1));
            RadarStack(kk,:) = mean(fxRadar{ii}(find(fxCatRadar{ii}{2}==kk),:),1);
        end
        %         end
        % Replace NaN with Zero (Bottom of Image has zero frequency Content)
        nanIx = find(isnan(RadarStack));
        RadarStack(nanIx) = 0;
        % Normalize Image
        RadarStack = trcNormalize(RadarStack);
        GPR.D.RadarStack{ii} = RadarStack;
    end
    clear('fxCatRadar','fxProRadar','fxRadar','stackIx')
end
if isSpatialFXdecon
    Window = 100; % Number of Traces in Fx-Decon Window Per Computation
    Wn = 2;      % Number of Overlapping Computations
    Overlap = round(Window/Wn);
    mu = 1e3;   % Weighting
    nfR = 50; % Length of Deconvolution Filter
    fxProRadar = cell(GPR.MD.nFiles,1); % Allocation

    for ii = 1:GPR.MD.nFiles
        dt = GPR.D.dt{ii};
        f0 = GPR.D.f0{ii};
        %             for jj = 1
        %             parfor (jj = chan, nWorkers)
        RadarStack = GPR.D.RadarStack{ii};
        tmpfx = 1:length(RadarStack);
        filterLap = Overlap:Overlap:length(tmpfx);
        % Design Data Window
        for kk = filterLap
            if kk == filterLap(1)
                run = 1;
                traceIx = 1:kk+(Overlap);
                tmpData = RadarStack(:,traceIx);
            elseif kk == filterLap(end)
                run = length(filterLap);
                traceIx = kk-(Overlap):length(tmpfx);
                tmpData = RadarStack(:,traceIx);
            else
                run = find(kk==filterLap);
                traceIx = kk-(Overlap):kk+(Overlap);
                tmpData = RadarStack(:,traceIx);
            end
            % Perform Fx-Decon within Windowed Data
            fxProRadar{ii}{run,1} = fx_deconT8(tmpData,dt,nfR,mu,f0);
            fxProRadar{ii}{run,2} = traceIx;
        end
        %             end

        fxRadar = cell(GPR.MD.nFiles,1);   % Allocate Temporary Arrays
        fxCatRadar = cell(GPR.MD.nFiles,1);% Temporary Array

        %             parfor (jj = chan, nWorkers)
        %             for jj = 1;
        % Concatenate Deconvolved Matricies
        fxCatRadar{ii}{1} = cat(2,fxProRadar{ii}{:,1});
        % Sort Trace Indicies
        [fxCatRadar{ii}{2},fxCatRadar{ii}{3}] = sort([fxProRadar{ii}{:,2}]);
        fxRadar{ii} = fxCatRadar{ii}{1}(:,fxCatRadar{ii}{3});

        % Fold Replicated Traces by Mean Stacking
        for kk = 1:length(RadarStack(1,:));
            RadarStack(:,kk) = mean(fxRadar{ii}(:,find(fxCatRadar{ii}{2}==kk)),2);
        end
        %             end
        % Replace NaN with Zero (Bottom of Image has zero frequency Content)
        nanIx = find(isnan(RadarStack));
        RadarStack(nanIx) = 0;
        % Normalize Image
        RadarStack = trcNormalize(RadarStack);
        GPR.D.RadarStack{ii} = RadarStack;
    end
end
clear('fxCatRadar','fxProRadar','fxRadar','out','stackIx')
% End of Function
end