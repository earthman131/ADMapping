function p=fdeninv(g,h1,h2,hm1,hm2,npts,nx,ny,dx,nn,ni)
G=6.67*1e-11;
xdiff=floor((npts-nx)/2); ydiff=floor((npts-ny)/2);
g1=taper2d(g,npts,ny,nx,ydiff,xdiff);
fg=fft2(g1); fg=fftshift(fg);
wn=2.0*pi/(dx*(npts-1));
cx=npts/2+1; cy=cx;
f0=1.0/dx;
for i=1:npts
    ky=(i-cy)*wn;
    for j=1:npts
        kx=(j-cx)*wn;
        k=sqrt(kx*kx+ky*ky);
        H=exp(-(k^2)/(2*(f0^2)));
        if k==0
            f1(i,j)=0;
        else
            f1(i,j)=H*fg(i,j)*k/(2*pi*G*(exp(-k*hm1)-exp(-k*hm2)));
        end
    end
end
f1=fftshift(f1); fp1=ifft2(f1); fp1=real(fp1);
fp0=fp1;
hh1=taper2d(h1,npts,ny,nx,ydiff,xdiff);
hh2=taper2d(h2,npts,ny,nx,ydiff,xdiff);
for J=1:ni
    f2=zeros(npts);
    for I=1:nn
        fh1=fft2(fp1.*(hh1.^I)); fh1=fftshift(fh1);
        fh2=fft2(fp1.*(hh2.^I)); fh2=fftshift(fh2);
        for i=1:npts
            ky=(i-cy)*wn;
            for j=1:npts
                kx=(j-cx)*wn;
                k=sqrt(kx*kx+ky*ky);
                H=exp(-(k^2)/(2*(f0^2)));
                if k==0
                    f2(i,j)=f2(i,j)+0;
                else
                    f2(i,j)=f2(i,j)+H*(k/(exp(-k*hm1)-exp(-k*hm2)))*((-k)^(I-1)...
                        /factorial(I))*(exp(-k*hm2)*fh2(i,j)-exp(-k*hm1)*fh1(i,j)); % ÂË²¨Ëã×Ó1
                end
            end
        end
    end
    f2=fftshift(f2); fp2=ifft2(f2); fp2=real(fp2);
    fp1=fp0-fp2;
end
p=fp1(1+ydiff:ny+ydiff,1+xdiff:nx+xdiff);
