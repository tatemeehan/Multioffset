function [t,p,L,raycoord]=traceray_pp(vp,zp,zs,zr,zd,x,xcap,pfan,itermax,optflag,pflag,dflag,kol)
% TRACERAY_PP: traces a P-P (or S-S) reflection for v(z)
%
% [t,p,L,raycoord]=traceray_pp(vp,zp,zs,zr,zd,x,xcap,pfan,itermax,optflag,pflag,dflag,kol)
%
% TRACERAY_PP uses the modified bisection algorithm 
% to trace a p-p reflection (or an s-s reflection)
% in a stratified medium.
%
% vp,zp   = velocity and depth vectors giving the p wave model. Depths are 
% 			considered to be the tops of homogeneous layers. 
% 	    Can also be used for an s-s reflection
% zs  	  = depth of the source (scalar)
% zr		  = depth of receivers (scalar)
% zd      = depth of the reflection (scalar)
% x       = vector of desired source-receiver offsets (horizontal)
%           MUST be non-negative numbers
% xcap 	  = (scalar) capture radius
% pfan    = vector of ray parameters to use for the initial fan of rays
%	    By default, the routine generates the fan using straight rays. 
%	    Setting pfan to -1 or omitting it activates default behavior.
%           Setting pfan to -2 causes the routine to use the pfan used in the
%           last call to this program. (pfan of -1 and -2 are identical on the
%	    first call.) 
% itermax = maximum number of iterations allowed for ray capture
% =========================== Default = 1 ================================
% optflag = if nonzero then refine captured rays by linear interpolation
% =========================== Default = 1 ================================
% pflag   = if nonzero, then print information about all failed rays
% =========================== Default = 0 ================================
% dflag   = 0 no action
%			 = 1 then a new figure window will be opened and the raypaths plotted
% 			 = 2 raypaths will be plotted in the current window
% =========================== Default = 0 ================================
% kol = color to draw rays
% =========================== Default = 'k' (black) =================
% NOTE: All z variables must be specified relative to a common datum
%
% t 	  = vector of traveltimes for each x
% p 	  = vector of ray parameters for each x
% L     = vector of geometrical spreading factors for each x
%     (if nargout<3, L will not be calculated)
% raycoord = cell array of same size as p. raycoord{k} is a n-by-2 matrix
%       whose first column is the ray depths at layer interfaces and the second
%       column is the ray x coordinates. Thus,
%       >> rc=raycoord{k};
%       >> plot(rc(:,1),rc(:,2));flipy
%       will plot the kth ray in the current figure.
%
% NOTE: failed rays are flagged with inf (infinity) for traveltime and 
%      nan for ray parameter.
%
% G.F. Margrave, CREWES Project, June 1995
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

if(nargin<13)
	kol='k';
end
if(nargin<12)
	dflag=0;
end
if(nargin<11)
	pflag=0;
end
if(nargin<10)
	optflag=1;
end

if(nargin<9)
	itermax=4;
end
if( nargin<8)
	pfan=-1;
end

if(~prod(double(size(vp)==size(zp))))
	error('vp and zp must be the same size');
end

if(length(zs)~=1 | length(zr)~=1 | length(zd)~=1 | length(xcap)~=1 )
	error(' zs,zr,zd and xcap must be scalars ')
end

ind=find(x<0);
if(~isempty(ind))
	error('offsets must be nonnegative');
end

%make sure zs < Zd and zr<zd
if(zs > zd )
	error(' zs must be less than zd');
end
if(zr > zd )
	error(' zr must be less than zd');
end

%adjust zs,zr,and zd and z2 so that they won't be exactly on layer boundaries
 	zs=zs+100000*eps;
 	zr=zr+100000*eps;
 	zd=zd-100000*eps;

 if( zs < zp(1) )
 	error(' source depth outside model range');
 elseif(zr < zp(1))
	error('receiver depth outside model range');
 end
 
 %determine the layers we propagate through
 %down leg
 ind=find(zp>zs);
 if(isempty(ind))
 	ibegd=length(zp);
 else
 	ibegd=ind(1)-1;
 end
 ind=find(zp>zd);
 if(isempty(ind))
 	iendd=length(zp);
 else
 	iendd=ind(1)-1;
 end
 
 %test for sensible layer indices
 if(iendd<1 | iendd > length(zp) | ibegd <1 | iendd > length(zp))
     %somethings wrong. return nans
     t=inf*ones(size(x));
     p=nan*ones(size(x));
     return;
 end
 
 %create v and z arrays for down leg
 vp=vp(:);zp=zp(:);
 vpd=[vp(ibegd:iendd);vp(iendd)];%last v is irrelevent.
 zpd=[zs;zp(ibegd+1:iendd);zd];

%up leg
 ind=find(zp>zr);
 if(isempty(ind))
 	ibegu=length(zp);
 else
 	ibegu=ind(1)-1;
 end
 ind=find(zp>zd);
 if(isempty(ind))
 	iendu=length(zp);
 else
 	iendu=ind(1)-1;
 end
 
 %test for sensible layer indices
 if(iendu<1 | iendu > length(zp) | ibegu <1 | iendu > length(zp))
     %somethings wrong. return nans
     t=inf*ones(size(x));
     p=nan*ones(size(x));
     return;
 end
 
 %create v and z arrays for up leg
 vpu=[vp(ibegu:iendu);vp(iendu)];%last v is irrelevent.
 zpu=[zr;zp(ibegu+1:iendu);zd];
 
 % combine into one model. We shoot oneway rays through an
 % equivalent model consiting of the down-leg model followed by
 % the upleg model (upside down).
 vp1=[vpd(1:end-1);flipud(vpu(1:end-1))];
 %vp1=[vpd;flipud(vpu(1:end-1))];
 tmp=flipud(zpu);
 zp1=[zpd;cumsum(abs(diff(tmp)))+zpd(end)];
 %zp1=[zpd;cumsum(abs(diff(tmp(1:end-1))))+zpd(end)];
 z1=zs;%start depth
 z2=max(zp1)+100000*eps;%end depth

 %determine pmax
 vmax=max([vp1]);
 vmin=min([vp1]);
 nearlyone=1-10*eps;
 %pmax=nearlyone/vmax;
 pmax=nearlyone/vmin;
 
 %the meaning of pmax is that it is the maximum ray parameter which will not
 %be critically refracted anywhere in vp or vs.

if(pfan==-2)
	global PFAN
	if(isempty(PFAN))
	  pfan=-1;
	else
	  pfan=PFAN;
	  PFAN=-2;
	end
else
	PFAN=0;
end

if(pfan==-1)
	%shoot a fan of rays Iterate once if the fan does not span the xrange
	%guess angle range
    nrays=max([3*length(x),10]);
    pfan=linspace(0,1/max(vp1),nrays);
else
	%make sure its a row vector
	ind=find(isnan(pfan));
	if(~isempty(ind))
		pfan(ind)=[];
	end
	pfan=pfan(:)';
	pfan=sort(pfan);
	ind=find(pfan==0);
	if(isempty(ind))
		pfan=[0 pfan];
	end
	ind=find(pfan>=pmax);
	if(~isempty(ind))
		pfan(ind)=[];
	end
	pfan=[pfan pmax];
	nrays=length(pfan);
end

%shoot first fan
% 
[xx,tt]=shootray(vp1,zp1,pfan);

xmax=max(xx);
xmin=min(xx);

%see if we have spanned the x range and revise if needed
imin=surround(xx,min(x));
imax=surround(xx,max(x));
newfan=0;
if(isempty(imax))
	if(max(pfan)<pmax)
		pm= .5*(max(pfan)+pmax);
		pfan=[pfan pm];
		imax=length(pfan)-1;
		newfan=1;
	end
end
if(~isempty(imax))
	pfan=pfan(imin:imax+1);
end

%shoot new fan (if p changed)
if(length(pfan)~=nrays | newfan)
	% first p
	[xx,tt]=shootray(vp1,zp1,pfan);
	% invoke symmetry to double these
	xmax=max(xx);
	xmin=min(xx);
end

%
% loop over x
%	-1 find xx which brackets x(k)
%	-2 shoot a finer fan of rays
%	-3 find the rays which bracket x(k)
%	-4 repeat 2 and 3 until a ray is captured at x(k)
%

t=inf*ones(size(x));
p=nan*ones(size(x));

for k=1:length(x)

	ind=surround(xx,x(k));
	if(isempty(ind))
		if(x(k)>=xmax) %test for off high end
			i1=length(xx);
			i2=[];
		else
			i1=[]; %off low end
			i2=1;
		end
	else
		if(length(ind)>1) %use first arrival
			it=find(tt(ind)==min(tt(ind)));
			ind=ind(it);
		end
		i1=ind;
		i2=ind+1;
	end
	
	%i1 and i2 are indices into pfan of the rays that bracket x(k)
	
	captured=0;
	hopeless=0;
	iter=0;
	x1=xx(i1); x2=xx(i2);
	t1=tt(i1); t2=tt(i2);
	p1=pfan(i1);  p2=pfan(i2);
	if( isempty(i2) & max(pfan)==pmax)
		hopeless=1;
	end
	
	while(~captured & iter<itermax & ~hopeless)
		
		iter=iter+1;
	
		%test for capture, x1 and x2 are the offsets for rays pfan(i1) and pfan(i2)
		xtst1=abs(x(k)-x1);
		xtst2=abs(x(k)-x2);
		if( xtst1< xcap & xtst1<=xtst2)
			captured=1;
			p(k)=p1;
			t(k)=t1;
			%final adjustment
			if(optflag)
				%linear interpolation
				t(k)= t1 + (x(k)-x1)*(t2-t1)/(x2-x1);
				p(k)= p1 + (x(k)-x1)*(p2-p1)/(x2-x1);
			end
					
			%disp([' capture on iter= ' int2str(iter)])
		elseif( xtst2 < xcap & xtst2<=xtst1)
			captured=1;
			p(k)=p2;
			t(k)=t2;
			%final adjustment
			if(optflag)
				%linear interpolation
				t(k)= t2 + (x(k)-x2)*(t1-t2)/(x1-x2);
				p(k)= p2 + (x(k)-x2)*(p1-p2)/(x1-x2);
			end
			
			%disp([' capture on iter= ' int2str(iter)])
		end
		
		%shoot a new fan
		if(~captured)
			%end cases first
			if(isempty(i1))
			    p2=p1;
				p1=0;
			elseif(isempty(i2))
				%p2=.5*(p2+pmax);
				p1=p2;
				p2=pmax;
			end
			dp= (x(k)-x1)*(p2-p1)/(x2-x1);
            if(dp<.1*(p2-p1)) dp=.1*(p2-p1); end
			p0= p1+ dp;
			%p2= min([p1+2*dp, pmax]);
			%p2=min([p1+2*dp p2]);
			%p1= p1+.1*dp;
			
			%p0=.5*(p1+p2);
			
			%pnew now contains only 1 new ray. The code below expects this
			pnew=[p1 p0 p2];
% 			nnew=3;
% 			pnew=linspace(min([p1 p2]) ,max([p1 p2]),nnew);
			
			if(isempty(pnew))
				hopeless=1;
			else
				%shoot and check for nans
                if length(pnew) == 1
                    pnew = [pnew,pnew,pnew];
                    [xtmp,ttmp]=shootray(vp1,zp1,pnew(2));
                else
				%[xxnew,ttnew]=shootray(vp1,zp1,pnew);
				[xtmp,ttmp]=shootray(vp1,zp1,pnew(2));
                end
				xxnew=[x1 xtmp x2];
				ttnew=[t1 ttmp t2];
	
				xmx=max(xxnew);
				xmn=min(xxnew);
				
				%analyze the fan. see if we have bracketed our target
				ind=surround(xxnew,x(k));
				if(isempty(ind))
					if(~isempty(xmx))
						if(x(k)>=xmx)
							n2=2*(x(k)-xxnew(2))/(xxnew(3)-xxnew(2));
							p2=min([p1+n2*dp, pmax]);
							pnew(3)=p2;
						end
						[xxnew,ttnew]=shootray(vp1,zp1,pnew);
	
						xmx=max(xxnew);
						xmn=min(xxnew);
						ind=surround(xxnew,x(k));
						if(isempty(ind));
							i1=3;
							i2=[];
						else
							i1=ind;
							i2=ind+1;
						end
					else
						i1=[];
						i2=1;
					end

				else
					if(length(ind)>1) %use first arrival
						it=find(ttnew(ind)==min(ttnew(ind)));
						ind=ind(it);
					end
					i1=ind;
					i2=ind+1;
				end
				if( isempty(i2) & max(pnew)==pmax)
					hopeless=1;
				else
					x1=xxnew(i1); x2=xxnew(i2);
					t1=ttnew(i1); t2=ttnew(i2);
					p1=pnew(i1);  p2=pnew(i2);
				end
			end
		end
		
	end
	if(hopeless)
	  disp(['hopeless after ' int2str(iter) ' iterations'])
	end
	%end the while loop
	if(~captured & pflag)
		disp([' after ' int2str(iter) ' iterations, capture failed for z1,z2= ' num2str(z1) ',' num2str(z2) ...
		' and x= ' num2str(x(k))]);
	end
	
	%end the for loop
end

if(PFAN)
	PFAN=pfan;
end

if(dflag)
	if(dflag==1)
		figure;
        	[hray,xray]=drawray(vpd,zpd,zs,zd,0,p,kol);
            [hray2,xray2]=drawray(vpu,zpu,zd,zr,xray,p,kol);
    else
        	[hray,xray]=drawray(vpd,zpd,zs,zd,0,p,'r');
            [hray2,xray2]=drawray(vpu,zpu,zd,zr,xray,p,'r');
	end

	if(dflag==1)
	    %#function flipy
		flipy;xlabel('offset');ylabel('depth');
		if(length(x)>1)
			dx=x(2)-x(1);
		else
			dx=100;
		end	
		axis([min([x(:);0])-dx max([x(:);0])+dx min(zp) max([zp;zs;zr;zd])])
	end
end

if(nargout>2)
   L=sphdiv(vp1,zp1,p);
end
if(nargout>3)
    %create ray coordinates
    raycoord=cell(1,length(p));
    for k=1:length(p)
        %coordinates for down leg
        sn=vpd(1:end-1)*p(k);%sin(theta)
        cs=sqrt(1-sn.^2);%cos(theta)
        xrayd=[0;cumsum(abs(diff(zpd)).*sn./cs)];
        %coordinates for up leg
        sn=flipud(vpu(1:end-1))*p(k);
        cs=sqrt(1-sn.^2);%cos(theta)
        xrayu=xrayd(end)+[0;cumsum(abs(diff(flipud(zpu))).*sn./cs)];
        xray=[xrayd;xrayu];
        zray=[zpd;flipud(zpu)];
        raycoord{k}=[xray zray]; 
    end
end