function [GPR] = processSWEPR(GPR,isTrimTWT,isReduceData,rmNtrc,padding,isKillChan,killArray)
for ii = 1 : GPR.MD.nFiles
    % Check Argin
    if nargin == 1
        % Controls
        isTrimTWT = 0;
        isReduceData = 0;
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        padding = 0;
        isKillChan = 0;
%         endsamp = 350-28
    elseif nargin == 2
        % Controls
        isReduceData = 0;
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        padding = 0;
        isKillChan = 0;
    elseif nargin == 3
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        padding = 0;
        isKillChan = 0;
    elseif nargin == 4
        padding = 0;
        isKillChan = 0;
    elseif nargin == 5
        isKillChan = 0;
    elseif nargin == 6
        error('Must Define Channels to Remove!')
    end

    % Extract Basic GPR paramters
    f0 = GPR.D.f0{ii}; dt = GPR.D.dt{ii};dx = GPR.D.dx{ii} ;
    if isKillChan
        rmvChan = find(ismember(GPR.D.trhd{ii}(3,:),killArray));
        GPR.D.trhd{ii}(:,rmvChan) = [];
        GPR.D.MxRadar{ii}(:,rmvChan) = [];
        GPR.Geometry.Chan{ii} = unique(GPR.D.trhd{ii}(3,:));
        GPR.Geometry.nChan{ii} = numel(GPR.Geometry.Chan{ii});
        GPR.Geometry.offset{ii} = GPR.Geometry.offset{ii}(GPR.Geometry.Chan{ii});
    end
    offset = GPR.Geometry.offset{ii};
    % Pad Data with Instrument Zero
    instrumentPad = zeros(padding,size(GPR.D.MxRadar{ii},2));
    if padding ~= 0
        for jj = 1:size(GPR.D.MxRadar{ii},2)
            instrumentZero = GPR.D.MxRadar{ii}(1,jj);
            instrumentPad(:,jj) = ones(padding,1).*instrumentZero;
        end
        GPR.D.MxRadar{ii} = [instrumentPad;GPR.D.MxRadar{ii}];
    end
    
    % Trim Time Window
    if isTrimTWT
        endSamp = 350-27; %200HHHV The new end of the columns
        reSample = endSamp + padding;   % Number of Wanted Samples
        GPR.D.MxRadar{ii} = GPR.D.MxRadar{ii}(1:reSample,:);
        TWT = GPR.D.TimeAxis{ii}; TWT = TWT(1:reSample);
    else
        % Store Travel-Time Axis
        TWT = GPR.D.TimeAxis{ii};
    end
    
    % Allocation Here
    if ii == 1
        Radar = cell(GPR.Geometry.nChan{ii},1,GPR.MD.nFiles); traceIx = cell(GPR.Geometry.nChan{ii},GPR.MD.nFiles);
    end
    TimeAxis = cell(GPR.Geometry.nChan{ii},1,1); traceIx = cell(GPR.Geometry.nChan{ii},GPR.MD.nFiles);
    
    for jj = 1:GPR.Geometry.nChan{ii}
        % DeMux Sequential Data
        % GPS DeadReckoning completed in preProcessing
        [Radar{jj,ii},~,traceIx{jj,ii},~] = deMuxNSIDC(GPR.D.MxRadar{ii},GPR.D.trhd{ii},GPR.Geometry.Chan{ii}(jj));
    end
        
    % parfor (jj =  1:GPR.Geometry.nChan{ii}, nWorkers)
    for jj = 1:GPR.Geometry.nChan{ii}
        % Reduce Data Volume
        if isReduceData
            Radar{jj,ii} = Radar{jj,ii}(:,1:rmNtrc:end);
            traceIx{jj,ii} = traceIx{jj,ii}(1:rmNtrc:end);
        end
        
        % Process Common Offset Channels

%         if ~ismember(jj,GPR.Geometry.badChan) % Hack for Bad Chans
        [Radar{jj,ii}, TimeAxis{jj},t0(jj)] = processCommonOffset(Radar{jj,ii}, ...
            f0, dt, TWT, offset(jj), GPR.Geolocation.Distance{ii}, dx, GPR.Geometry.Chan{ii}(jj) );
%         else
%             TimeAxis{jj} = TWT; t0(jj) = 0;
%         end
    end
    if isReduceData
        tmpIx = sort(cat(2,traceIx{:,ii}));
        GPR.D.trhd{ii} = GPR.D.trhd{ii}(:,tmpIx);
    end
    
    % Trim all time windows to same length
    for jj = 1:GPR.Geometry.nChan{ii}
        if jj == 1
            minIx = length(TimeAxis{jj});
            minChan = jj;
        elseif minIx > length(TimeAxis{jj})
            minIx = length(TimeAxis{jj});
            minChan = jj;
        end
    end
    clear('TimeAxis')
    % Store Travel-Time Axis
    TimeAxis = [0:dt:(dt.*(minIx-1))]';
    % Trim Channels to Consistent Travel Time Axis
    for jj = 1:GPR.Geometry.nChan{ii}
        Radar{jj,ii} = Radar{jj,ii}(1:minIx,:);
    end
    % Write to Out Structure
    GPR.D.TimeAxis{ii} = TimeAxis;
    GPR.D.t0{ii} = t0;
    GPR.D.Radar{ii} = {Radar{:,ii}}; 
end
end

