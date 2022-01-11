function g=fdenfor(p,h1,h2,hm1,hm2,npts,nx,ny,dx,nn)
G=6.67*1e-11;
xdiff=floor((npts-nx)/2); ydiff=floor((npts-ny)/2);
p1=taper2d(p,npts,ny,nx,ydiff,xdiff);
fp=fft2(p1); fp=fftshift(fp);
wn=2.0*pi/(dx*(npts-1));
cx=npts/2+1; cy=cx;
for i=1:npts
    ky=(i-cy)*wn;
    for j=1:npts
        kx=(j-cx)*wn;
        k=sqrt(kx*kx+ky*ky);
        if k==0
            f1(i,j)=0;
        else
            f1(i,j)=(2*pi*G/k)*(exp(-k*hm1)-exp(-k*hm2))*fp(i,j);
        end
    end
end
hh1=taper2d(h1,npts,ny,nx,ydiff,xdiff);
hh2=taper2d(h2,npts,ny,nx,ydiff,xdiff);
f2=zeros(npts);
for I=1:nn
    fh1=fft2(p1.*(hh1.^I)); fh1=fftshift(fh1);
    fh2=fft2(p1.*(hh2.^I)); fh2=fftshift(fh2);
    for i=1:npts
        ky=(i-cy)*wn;
        for j=1:npts
            kx=(j-cx)*wn;
            k=sqrt(kx*kx+ky*ky);
            f2(i,j)=f2(i,j)+((-k)^(I-1)/factorial(I))*(exp(-k*hm2)*fh2(i,j)-exp(-k*hm1)*fh1(i,j));
        end
    end
end
f2=f2*(2*pi*G);
fg=f1+f2;
fg=fftshift(fg); fginv=ifft2(fg);
g=real(fginv(1+ydiff:ny+ydiff,1+xdiff:nx+xdiff));
