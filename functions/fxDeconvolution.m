function [GPR] = fxDeconvolution(GPR)
for ii = 1:GPR.MD.nFiles
% Extract Variables
dt = GPR.D.dt{ii};
f0 = GPR.D.f0{ii};
nChan = GPR.Geometry.nChan{ii};
Radar = GPR.D.Radar{ii};
    % Spatial Predictive Filtering for Random Noise Suression
    Window = round(2/dt); % Number of Temporal Samples in Fx-Decon Window Per Computation
    Overlap = round(Window/10); % Number of Overlapping Computations
    mu = 1e8; % Weighting
    nfR = 15; % Length of Deconvolution Filter
    fxProRadar = cell(nChan,1); % Allocation
    % Begin the Processing
    %         parfor (jj = chan, nWorkers)
        for jj = 1:nChan
            [m,n] = size(Radar{jj});
            filterLap = Overlap:Overlap:m;
            % Design Data Window
            for kk = filterLap
                if kk == filterLap(1)
                    run = 1;
                    traceIx = 1:kk+(Overlap);
                    tmpData = Radar{jj}(traceIx,:);
                elseif kk == filterLap(end)
                    run = length(filterLap);
                    traceIx = kk-(Overlap):m;
                    tmpData = Radar{jj}(traceIx,:);
                else
                    run = find(kk==filterLap);
                    traceIx = kk-(Overlap):kk+(Overlap);
                    tmpData = Radar{jj}(traceIx,:);
                end
                % Perform Fx-Decon within Windowed Data
                fxProRadar{jj}{run,1} = fx_deconT8(tmpData,dt,nfR,mu,f0);
                fxProRadar{jj}{run,2} = traceIx;
            end
%         end
        
        fxRadar = cell(nChan,1);   % Allocate Temporary Arrays
        fxCatRadar = cell(nChan,1);% Temporary Array
        
%         parfor (jj = chan, nWorkers)
            % Concatenate Deconvolved Matricies
            fxCatRadar{jj}{1} = cat(1,fxProRadar{jj}{:,1});
            % Sort Trace Indicies
            [fxCatRadar{jj}{2},fxCatRadar{jj}{3}] = sort([fxProRadar{jj}{:,2}]);
            fxRadar{jj} = fxCatRadar{jj}{1}(fxCatRadar{jj}{3},:);

            % Fold Replicated Traces by Mean Stacking
            for kk = 1:length(Radar{jj}(:,1))
                Radar{jj}(kk,:) = mean(fxRadar{jj}(find(fxCatRadar{jj}{2}==kk),:),1);
            end
%         end
            % Replace NaN with Zero (Bottom of Image has zero frequency Content)
            nanIx = find(isnan(Radar{jj}));
            Radar{jj}(nanIx) = 0;
            % Normalize Image
            Radar{jj} = trcNormalize(Radar{jj});
        end
        GPR.D.Radar{ii} = Radar;
    clear('fxCatRadar','fxProRadar','fxRadar','stackIx')
end
end

