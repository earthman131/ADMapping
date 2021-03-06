%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A wavenumber-domain iterative approach for apparent density mapping of an undulant layer
% Author: Lianghui Guo (guolh@cugb.edu.cn)
% Organization: China University of Geosciences (Beijing), School of Geophysics and Information Technology
% Compiled version: MATLAB R2017b
% Reference:
%       Guo L H, Shi L, Li S L. A wavenumber-domain iterative approach for apparent density
%       mapping of an undulant layer and its application in central south China. Geophysics, 2019,
%       84(1): G1-G11.
% Description of the input parameters: 
%       infile_ano??observed anomaly data file, unit: mGal
%       infile_h1??top depth data file, unit: m
%       infile_h2??bottom depth data file, unit: m
%		pb??2.67??background density, unit: g/cm^3
%       h0??height of observation surface, unit: m
%       nn??series
%       ni??number of internal iterations
%       iter??number of external iterations
% Description of the output parameters: 
%       outfile_p??apparent density, unit: g/cm^3
%       outfile_cal??calculated anomaly, unit: mGal
%       outfile_res??residual anomaly, unit: mGal
%       outfile_rms??root mean square (rms), unit: mGal
% Description of primary identifiers??
%       x, y: x, y verctor
%       nx, ny: number of points in x and y directions
%       dx, dy: spacing in x and y directions
%       obs: observed anomaly data
%       h1: top depth
%       h2: bottom depth
%       npts: extension points
%       pcal: apparent density
%       cal: calculated anomaly
%       res: residual anomaly
%       error??root mean square (rms)
% Description of subroutine function: 
%       readgrd.m: read surfer text grd file
%       fdeninv.m: inversion    
%       fdenfor.m: forward
%       taper2d.m: cosine attenuation edge extension
%       savebln.m: save rms file
%       savegrd.m: save surfer text grd file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%% I/O parameters %%%%%%%%%%%%%
infile_ano = 'anomaly.grd';
infile_h1 = 'h1.grd';
infile_h2 = 'h2.grd';
pb = 0;
h0 = 0; 
nn = 5;
ni = 3; 
iter = 10; 
outfile_p = 'pn_g_inv10.grd'; 
outfile_cal = 'pn_g_inv10_cal.grd'; 
outfile_res = 'pn_g_inv10_g_res.grd'; 
outfile_rms = 'pn_g_inv10_rms.bln';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[obs,x,y,nx,ny,dx,dy] = readgrd( infile_ano );obs=obs*1e-5;g1=obs;
[h1,~,~,~,~,~,~] = readgrd( infile_h1 );h1=h1-h0;hm1=mean(mean(h1));
[h2,~,~,~,~,~,~] = readgrd( infile_h2 );h2=h2-h0;hm2=mean(mean(h2));
[X,Y] = meshgrid(x,y);
nmax = max([nx ny]);
npts = 2^nextpow2(nmax);
pcal = zeros(ny,nx);
error = zeros(iter,1);
% inversion
for i = 1:iter
    p = fdeninv(g1,h1,h2,hm1,hm2,npts,nx,ny,dx,nn,ni); 
    g2 = fdenfor(p,h1,h2,hm1,hm2,npts,nx,ny,dx,nn); 
    g1 = g1-g2;
    pcal = pcal+p;
    cal = fdenfor(pcal,h1,h2,hm1,hm2,npts,nx,ny,dx,nn); 
    res = obs-cal;
    error(i) = rms(res(:));
end
error = error*1e5;
cal = cal*1e5;
res = res*1e5;
pcal = pcal/1000;
pcal = pcal+pb;
% save
savegrd(pcal,x,y,nx,ny,outfile_p)
savegrd(cal,x,y,nx,ny,outfile_cal)
savegrd(res,x,y,nx,ny,outfile_res)
savebln(error,outfile_rms) 