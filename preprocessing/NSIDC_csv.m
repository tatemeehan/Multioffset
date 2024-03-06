%% Save .csv for NSIDC
dataDir = 'D:\GPRdata\031721\IDBRBO031721\processed';
nsidcDir = 'D:\GPRdata\031721\IDBRBO031721\NSIDC\';
[status, msg, msgID] = mkdir(nsidcDir);
files = dir([dataDir,'\*.csv']);
nFiles = length(files);
for ii = 2:nFiles
    T = readtable([files(ii).folder,'\',files(ii).name]);
    mm = month(T.DateTime);
    dd = day(T.DateTime);
    yy = year(T.DateTime);
    if mm < 10
        MM = [repmat('0',numel(mm),1),num2str(mm)];
    else
        MM = num2str(mm);
    end
    if dd < 10
        DD = [repmat('0',numel(dd),1),num2str(dd)];
    else
        DD = num2str(dd);
    end
    YY = num2str(yy);
    MMDDYY = [MM,DD,YY];
    for kk = 1:length(MMDDYY)
    MMDDYY(kk,:) = regexprep(MMDDYY(kk,:), ' ', '0');
    end
    Date = num2str([mm,dd,yy]);
    [H,M,S] = hms(T.DateTime);
    S = round(S);
    if H < 10
        HH = [repmat('0',numel(H),1),num2str(H)];
    else
        HH = num2str(H);
    end
    if M < 10
        MN = [repmat('0',numel(M),1),num2str(M)];
    else
        MN = num2str(M);
    end
    if S < 10
        SS = [repmat('0',numel(S),1),num2str(S)];
    else
        SS = num2str(S);
    end
    HHMMSS = [HH,MN,SS];
    for kk = 1:length(HHMMSS)
    HHMMSS(kk,:) = regexprep(HHMMSS(kk,:), ' ', '0');
    end
    existCol = strcmp('iceTWT',T.Properties.VariableNames);
    val = existCol(existCol==1) ;
    if val

        % Density
        rho = 310;
        % Depth
        v = DryCrimVRMS(rho);
        Depth = 100.*(T.snowTWT.*(v)./2); % cm
        % SWE
        SWE = rho.*Depth./100;
        % Comments
        Comments = strings(size(T.snowTWT));

        rhoIce = 917;
        vIce = DryCrimVRMS(rhoIce);
        iceThickness = 100.*T.iceTWT.*(vIce)./2;
        iceWE = iceThickness.*rhoIce./100;
        % Output Table
        Tnisdc = table(MMDDYY,HHMMSS,T.Longitude,T.Latitude,T.ElevationEMG96,T.Northing,T.Easting,T.UTMzone,T.snowTWT,Depth,SWE,T.iceTWT,iceThickness,iceWE,Comments,'VariableNames',["Date [mmddyy]","Time [HHMMSS]","Longitude [DD]","Latitude [DD]","Elevation [masl]","Northing","Easting","UTM_Zone","TWT [ns]","Depth [cm]","SWE [mm]","Ice TWT [ns]","Ice Thickness [cm]","IceWE [mm]","comments"]);
        writetable(Tnisdc,[nsidcDir,files(ii).name])
    else
        % Density
        rho = 310;
        % Depth
        v = DryCrimVRMS(rho);
        Depth = 100.*(T.TWT.*(v)./2); % cm
        % SWE
        SWE = rho.*Depth./100;
        % Comments
        Comments = strings(size(T.TWT));
        % Output Table
        Tnisdc = table(MMDDYY,HHMMSS,T.Longitude,T.Latitude,T.ElevationEMG96,T.Northing,T.Easting,T.UTMzone,T.TWT,Depth,SWE,Comments,'VariableNames',["Date [mmddyy]","Time [HHMMSS]","Longitude [DD]","Latitude [DD]","Elevation [masl]","Northing","Easting","UTM_Zone","TWT [ns]","Depth [cm]","SWE [mm]","comments"]);
        writetable(Tnisdc,[nsidcDir,files(ii).name])
    end
end

pitDensity = [301,290;248,268;296,292;309,332;324,357];
mean(pitDensity(:))