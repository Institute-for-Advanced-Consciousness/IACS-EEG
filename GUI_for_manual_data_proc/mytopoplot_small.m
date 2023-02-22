function []=mytopoplot(dat,pos,interpmethod,interpfunction,npoints,isolines)

% GENERATES TOPOPLOT FROM THE INTENSITIES MEASURED AT A GIVEN SET OF
% CHANNELS

% INPUT:
% -  dat: Nchannels x 1 vector with the intensities to plot
% -  pos: Nchannels x 2 matrix with the channel positions (x: ...
%       left to right, y: occipital to frontal) 
% - (optional) interpmethod: Interpolation method: possible values:
%       If using griddata:
%                   'nearest'   - Nearest neighbor interpolation
%                   'linear'    - Linear interpolation (default)
%                   'natural'   - Natural neighbor interpolation
%                   'cubic'     - Cubic interpolation 
%                   'v4'        - MATLAB 4 griddata method 
%           (Note: I disrecomment linear, natural,cubic because area 
%           outside channels is not extrapolated - ugly plot)
%       If using TriScatteredInterp only 'natural','linear' or 'nearest'
% - (optional) interpfunction : matlab function to use to interpolate: 
%       'TriScatteredInterp' or 'griddata' (default = griddata)
%  - (optional) npoints: number of points (or pixels) of the interpolated 
%       image (default = 100)
%  - (optional) isolines: number of isolines to plot (default = 0 : no
%       isolines are plotted) 
            
% OUTPUT 
% No output. The topoplot is plotted in current figure (if existent)

% Last change 2016.08.19 P.Garces
            
%% PARAMETERS

if nargin < 6 || isempty(isolines)
    isolines=0; 
elseif ~isnumeric(isolines) || isolines<0
    fprintf('Isolines must be an integer >= 0 - using 0 isolines \n');
    isolines=0; 
else
    isolines = round(isolines);
end

if nargin < 5 || isempty(npoints)
    npoints=100;
elseif ~isnumeric(npoints) || npoints<=4
    fprintf('npoints must be an integer >= 5 - using 100 npoints \n');
    npoints=100;
else
    npoints = round(npoints);
end


if (nargin<4) || isempty(interpfunction)
    interpfunction='griddata';
elseif ~ismember(interpfunction,{'TriScatteredInterp' ,'griddata'})
    fprintf('Desired interpolation function not possible - using griddata \n')
    interpfunction='griddata';
end


if nargin < 3 || isempty(interpmethod)
    if strcmpi(interpfunction,'griddata')
        interpmethod = 'v4';
    else
        interpmethod = 'nearest';
    end
elseif strcmpi(interpmethod,'griddata') && ~ismember(interpmethod,{'nearest' ,'linear','natural','cubic','v4' })
    interpmethod = 'natural';
    fprintf('Desired interpolation method not possible - using natural \n')
elseif strcmpi(interpmethod,'griddata') && ~ismember(interpmethod,{'nearest' ,'linear','natural'})
    interpmethod = 'natural';
    fprintf('Desired interpolation method not possible - using natural \n')
end

if nargin<2
    error('minimum dat and pos inputs are required')
end

dat = dat(:);
if length(dat)~=size(pos,1)
   error('dimensions do not match:  length(dat)~=size(pos,1)')
end


% Channel positions
chanX=pos(:,1);
chanY=pos(:,2);


%% PREPARE MASK AND TOPOPLOT OUTLINE

center=[mean(chanX),mean(chanY)];
Nptsmask=npoints*15;
phi=linspace(0,2*pi,Nptsmask)';

% eliptica
growing_factor=[1.1 1.1];
paramhor=growing_factor(1)*max(abs(chanX-center(1))); 
paramver=growing_factor(1)*max(abs(chanY-center(2)));
mask=[center(1)+paramhor*cos(phi), center(2)+paramver*sin(phi)];

% Creates outline: elipse + nose + ears
outline=cell(1,1);
outline{1}=mask;
% Nose
phi=(pi/2-1/8:1/8:pi/2+1/8)';
outline{2}=[center(1)+paramhor*cos(phi), center(2)+paramver*sin(phi)];
outline{2}(2,2)=outline{2}(2,2)+paramver*1/10;
% Ear
phi2=linspace(pi/2,3*pi/2,20)'; earradius=paramver/10;
oreja=[earradius*cos(phi2), earradius*sin(phi2)]+kron([center(1)-paramhor,center(2)],ones(length(phi2),1));
outline{3}=oreja;
oreja=[-earradius*cos(phi2), earradius*sin(phi2)]+kron([center(1)+paramhor,center(2)],ones(length(phi2),1));
outline{4}=oreja;


%% PREPARES GRID FOR INTERPOLATION

% Extremum coordinates
vlim=[min(mask(:,2)) max(mask(:,2))];
hlim=[min(mask(:,1)) max(mask(:,1))];

xi        = linspace(hlim(1), hlim(2), npoints);   % eje x (Npuntos x 1)
yi        = linspace(vlim(1), vlim(2), npoints);   % eje y (Npuntos x 1)
[Xi,Yi]   = meshgrid(xi', yi);

% Mask: Determines which points are inside the topoplot outline
maskimage= false(npoints);
for ipt=1:npoints*npoints
    myangle=atan2(mask(:,2)-Yi(ipt),mask(:,1)-Xi(ipt));
    % unwrapping of myangle
    d=diff(myangle);
    indx = find(abs(d)>pi);
    for ijump=indx(:)'
        if d(ijump)>0
            myangle((ijump+1):end) = myangle((ijump+1):end) - 2*pi;
        else
            myangle((ijump+1):end) = myangle((ijump+1):end) + 2*pi;
        end
    end
    maskimage(ipt)=max(myangle)-min(myangle)>3*pi/2;
end


%% ESTIMATES DATA INTENSITIES AT INTERPOLATION POINTS (Zi)

if strcmp(interpfunction,'TriScatteredInterp')
    F = TriScatteredInterp(chanX,chanY,dat,interpmethod);
    [Xi,Yi] = meshgrid(xi,yi);
    Zi = F(Xi,Yi);
else
    [Xi,Yi,Zi] = griddata(chanX', chanY, dat, xi', yi, interpmethod);
end
% Sets points ouside the mask to Nan
if ~isempty(maskimage)
    Zi(~maskimage) = NaN; 
end

%% PLOTS

% uses current figure
 hold on

 
% - plots color
% OLD OPTIONS - Using surf
% shading='flat';
% deltax = xi(2)-xi(1); % length of grid entry
% deltay = yi(2)-yi(1); % length of grid entry
% surf(Xi-deltax/2,Yi-deltay/2,zeros(size(Zi)), Zi, 'EdgeColor', 'none', 'FaceColor', shading);
% NEW OPTION - Using imagesc
myim = imagesc(xi,yi,Zi);
set(myim,'AlphaData',~isnan(Zi));

% Plots outline
for i=1:length(outline)
    xval = outline{i}(:,1);
    yval = outline{i}(:,2);
    plot(xval, yval, 'k', 'LineWidth',0.5)
end

% % Crea isolineas
% if isolines>0
%     contour(Xi,Yi,Zi,isolines,'k');
% end
% 
% 
% % Plots sensor positions
% % plot(chanX,chanY,'k.','MarkerSize',1);
% plot(chanX,chanY,'ko','MarkerSize',1,'MarkerFaceColor',[0 0 0]);
% 
% axis square; axis off
% colorbar
% 
% 
% xlim([min(outline{3}(:,1)) max(outline{4}(:,1))]);
% ylim([min(outline{1}(:,2)) max(outline{2}(:,2))]);


