function [WAVE, scale]=wavelet_comp(X,No_of_scales)

Y=[zeros(1,6) X zeros(1,6)];

fs=250;
for j=1:100
    fx=(j-1)*0.5+1.5;
    scale(j)=fs/fx;
end

n=512;
dt=1/250;
pad=0;
dj=0.25;
J1=0;

mother='MORLET';
param=6;

for i=1:No_of_scales

    s0=dt*scale(i)*0.9680133;
[WAVE(i,:),PERIOD(i,:),SCALE(i,:),COI(i,:)]=wavelet(Y,dt,pad,dj,s0,J1,mother,param);

end