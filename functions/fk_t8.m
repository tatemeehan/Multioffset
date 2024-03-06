function [ fkRad ] = fk_t8( rad, dt, dx )
% fk_T8 Performs the Half-Nyquist WaveNumber Filter on NMO Corrected CMPs
%   Operation of this filter follows the outline provieded within Oz Yilmaz 
%   Seismic Data Analysis, Ch. 1.2.

%   Inputs - rad: The NMO Corrected CMP Gather
%             dt: Time Sample Interval [ns]
%             dx: CMP Trace Interval [m]
%
%   Output - fkRad: The f-k Filtered CMP Gather
%
%   Written by Tate Meehan, Boise State University, GreenTrACS 2017

% Pad the input data array with zeros to smoothly interpolate fft  
[ns, ntr] = size(rad); 
ns2 = nextpow2(ns); ntr2 = nextpow2(ntr);
padtop = nan(floor((2.^ns2-ns)./2),2.^ntr2);
padbot = nan(ceil((2.^ns2-ns)./2),2.^ntr2);
padl = nan(ns,floor((2.^ntr2-ntr)./2));
padr = nan(ns,ceil((2.^ntr2-ntr)./2));
Rad = rad; % Unpadded CMP
% Pad Rad
Rad = [padl,Rad,padr]; Rad = [padtop;Rad;padbot];
radIx = find(~isnan(Rad));
Rad(~ismember(1:numel(Rad),radIx)) = 0;
% FFT2
Rad = fft2(Rad); 
Rad = fftshift(Rad); 
% imresize
% Rad = imresize(Rad,[size(Rad,1),size(Rad,2)./4]);

%----------------------------------------------------------------------- 
%  Useful parameters 
[nf,nk] = size(Rad); 
nfhalf  = nf/2; 
fn      = 1./(2*dt);                 % Nyquist frequency 
df      = 2*fn/(nf - 1);             % frequency interval 
ff      = -fn:df:fn;                 % frequency spectrum 
nkhalf  = nk/2; 
if dx < 0 
    dx = -dx; 
end 
kn      = 1./(2*dx);                 % Nyquist wavenumber
knhalf  = kn./2;                     % Half Nyquist Wavenumber
dk      = 2*kn/(nk - 1);             % wavenumber interval 
kk      = -kn:dk:kn;                 % wavenumber spectrum   
% Build Diamond Window
[~,fkLix] = min(abs(kk+kn./8.5));
[~,fkRix] = min(abs(kk-kn./8.5));
[~,ffTix] = min(abs(ff-fn./1.75));
[~,ffBix] = min(abs(ff+fn./1.75));
ffCix = length(ff)./2;
fkCix = length(kk)./2;
% Extract Indicies
fix = length(ff)./2:-1:1;
fix = fix(fix>ffBix & fix <ffTix);
fix = [fix,fliplr(fix)];
kix = 1:length(kk);
kix = kix(kix>fkLix & kix<fkRix);
fix = round(interp1(linspace(1,length(kix),length(fix)),fix,1:length(kix)));
% Create Mask
mask = ones(size(Rad,1)./2,size(Rad,2));
    tmpfix = fix(1:length(kix)./2);
    tmpkix = kix(1:length(kix)./2);
for ii = 1:length(kix)
    tmp = fix(ii)-1:-1:1;
    try    ix = sub2ind(size(mask),tmp,kix(ii).*ones(size(tmp)));
    catch
        keyboard
    end

    mask(ix) = 0;
end
mask(:,~ismember(1:size(Rad,2),kix)) = 0;
mask = [mask;flipud(mask)];
win = hamming(round(size(Rad,1).*7.5))';
mask = conv2(mask,win,'same')./sum(win);
iskeyhole = 0;
isdoublediamond = 1;
if iskeyhole
    % Resize 1/2.^3
smallmask = imresize(mask,.125);
% The Steps of Padding 2.^3
smallmask = (padarray(smallmask,[size(smallmask)./2],'both'));
smallmask = (padarray(smallmask,[size(smallmask)./2],'both'));
smallmask = (padarray(smallmask,[size(smallmask)./2],'both'));
% Invert Weights
smallmask = abs(smallmask - 1);
% Average the Masks
theMask = (mask+smallmask)./2;
elseif isdoublediamond
    % Resize 1/2
    smallmask = imresize(mask,[size(mask,1)./2,size(mask,2)]);
    % Double Diamond
    theMask = [smallmask;smallmask];
    % Shift Padding from Center of Diamonds to Top and Bottom
    tmp = theMask;
    tmp(size(theMask,1)./2+1 - size(theMask,1)./8:size(theMask,1)./2 + size(theMask,1)./8,:) = [];
    theMask = [zeros(size(theMask,1)./8,size(theMask,2));tmp;zeros(size(theMask,1)./8,size(theMask,2))];
else
    % Diamond Mask
    theMask = mask;
end
fkRad = Rad.*theMask;
% fkRad = Rad.*mask;
% % Build Half-Nyquist WaveNumber Filter
% fkX = find(kk <= knhalf & kk >= -knhalf);
% fkH = zeros(size(Rad)); fkH(:,fkX) = 1;
% fkH = conv2(fkH,hamming(ntr)','same'); fkH = fkH./max(fkH(:));
% % Apply WaveNumber Filter
% fkRad = Rad.*fkH;

% imresize
% fkRad = imresize(fkRad,[size(fkRad,1),size(fkRad,2).*4]);

% Perform the Inverse FFT
 fkRad = real(ifft2(ifftshift(fkRad)));
 % Grab Original Data Gather without Padding
 fkRad = reshape(fkRad(radIx),size(rad));
%  fkRad = fkRad(1:ns,1:ntr); % Grab Original CMP Gather without Padding
 
% Plot f-k Plane
%     ffk = figure('Numbertitle','off','tag','fkfiltfig',... 
%         'Name','F-K spectrum','menubar','none'); 
% % Display every second element of the F-K spectrum for faster rendering 
%     pcolor(kk(1:1:nk),ff(nfhalf+1:1:nf),log10(abs(Rad(nfhalf+1:1:nf,1:1:nk))));  
%     shading('flat'); 
%     set(gca,'ydir','reverse'); 
%     ylabel('Frequency ( GHz )'); 
%     xlabel('Wavenumber (m^-^1)') 


end

