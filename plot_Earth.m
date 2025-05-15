function [gh] = plot_Earth(npts,GST)
% Plots the Earth in the orientation determined by GST, the Greenwich Sidereal
% Time.
% 
% 
% Note: This function requires an image of the Earth (a texturemap) to plot on
% the ellipsoidal surface. You can get one you like from:
% https://visibleearth.nasa.gov/collection/1484/blue-marble
% 
% The one included with this function is the map from July:
% https://visibleearth.nasa.gov/images/73751/july-blue-marble-next-generation-w-
% topography-and-bathymetry
% 
% Inputs:
%   npts: [](1x1)<integer> number of equatorial points to use to generate the
%         surface of the Earth
%   GST: [deg](1x1)<double> Greenwich Sidereal Time, the counter-clockwise angle
%        between the x-axis (J2000 vernal equinox) and Greenwich, England
% 
% Outputs:
%   gh: [](1x1)<graphics handle> A graphics handle that contains data about the
%        plot if needed
% Equatorial and polar radii of the Earth:
RE_eq = 6378.137; %[km] Earth equatorial radius (WGS 84)
RE_pol = 6356.752314245; %[km] Earth polar radius (WGS 84)

% Get the points on the ellipsoid defined by the two radii:
[Earthx,Earthy,Earthz] = ellipsoid(0,0,0,RE_eq,RE_eq,RE_pol,npts);

% Rotate the points to orient the Earth properly:
C = [cosd(GST),-sind(GST);
     sind(GST),cosd(GST)];
Rotated = C*[Earthx(:).';Earthy(:).'];

% Reshape the x and y matrices to be the same size as the z matrices
Earthx = reshape(Rotated(1,:),npts+1,npts+1);
Earthy = reshape(Rotated(2,:),npts+1,npts+1);   

% Read in the texturemap for the Earth
cdata = imread('Earth_texturemap.jpg');

% Plot the Earth's surface and use the texturemap as the color:
gh = surf(Earthx,Earthy,-Earthz,'FaceColor','texturemap','CData',cdata, ...
    'EdgeColor','none', 'DisplayName','Earth','FaceAlpha',1);
axis equal