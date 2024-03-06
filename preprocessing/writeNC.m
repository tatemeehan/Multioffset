function writeNC(trhd,Rad,filename,writeDir,workingDirectory)
%% Write to NetCDF
% Attach trace header to data for NetCDF Export
DATA = [trhd;Rad];

cd(writeDir)
ncdfName = filename;

if length(strsplit(ncdfName,'.')) == 1
    ncdfName = [ncdfName,'.nc'];

fprintf(['Writing ',ncdfName, '\n'])

numrow = size(DATA,1);

numcol = size(DATA,2);

% Create the NetCDF-4 file with Dimensions of DATA
netcdf.setDefaultFormat('FORMAT_NETCDF4') ;

ncid = netcdf.create(ncdfName,'CLOBBER');

dimidrow = netcdf.defDim(ncid,'rows',numrow);

dimidcol = netcdf.defDim(ncid,'length',numcol);

% Define and put Variable DATA
varid = netcdf.defVar(ncid,'DATA','NC_DOUBLE',[dimidrow dimidcol]);

% Increase NetCDF Compression
netcdf.defVarDeflate(ncid,varid,true,true,5);

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid,DATA);

% Close the .nc file after writing
netcdf.close(ncid);

% Return to the working directory for next iteration
cd(workingDirectory)
fprintf(['Wrote ',ncdfName,'\n'])

display(' ')
end

