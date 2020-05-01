%DRX MODEL: written by Cram DG (2012)
%---growth.m---
%1. AVERAGE DISLOCATION DENSITY OVER SURFACE AREA
sumpdsq=sum((p(a).*((d(a)./2).^2))); %sum of dislocations*r^2
sumdsq=sum((d(a)./2).^2); %sum of r^2
sumd=sum(d(a)./2); %sum of r
avgdislosur=sumpdsq/sumdsq; %surface area average dislocation density
% COMPLEXION AREA
fco(a) = cplx_area(a)./(pi*d(a).^2);
fco(fco>1)=1;
total_cplx_area = sum(cplx_area(a));
%2. GRAIN VELOCITY + GRAIN SIZE
ov(a)=v(a); %saves old volume
if UseCompGrow 
    cmean_mob(a) = mob.*(1-fco(a))+ mco.*fco(a);
    vel(a)=cmean_mob(a).*(0.5*shearmod*(burger^2).*(avgdislosur-p(a)));
else
    vel(a)=mob*(0.5*shearmod*(burger^2)*(avgdislosur-p(a))); %grain velocity
end
d(a)=d(a)+(2*vel(a).*timeinc); %updating grain diameter
v(a)=(4/3)*pi()*((d(a)./2).^3); %updating volume
%3. UPDATING GRAIN DISLOCATIONS + DESTROYING GRAINS
%a) for negative velocity
%if grain shrinks: no change to the dislocation density
des=0; %resets count: number of destroyed grains in this increment de
ret=find(d(a)<grainsizelimit); %finds live grains that are too small (less than set limit)
ret=a(ret); %matrix containing the actual grain numbers which are too small
b(ret)=0; %sets grain status to dead
%resetting grain values to 0 or NaN (for identification purposes/plots)
nuc(ret)=0; %sets grain nucleation potential to 0
r(ret)=NaN; %sets grain size to NaN
y(ret)=NaN; %sets yield stress to NaN
sg(ret)=NaN; %sets subgrain size to NaN
rc(ret)=NaN; %sets critical subgrain size to NaN
d(ret)=NaN; %sets grain size to NaN
v(ret)=0; %sets grain volume to 0
desn=size(ret); %size of ret (retire) to obtain number grains destroyed/dead in strain step
des=desn(1,2); %number grains destroyed/dead in strain step
if (des>0) %if any grains have been destroyed then:
n=n-des; %adjusts the number of live grains
tdes=tdes+des; %adjusts the total number of grains that have ever been destroyed
a=find(b==1); %updates the live grain number matrix
end
%b) if grain velocity is positive
grow=find(vel(a)>0); %finds live grains whose velocity is positive
grow=a(grow); %matrix containing the actual grain numbers whose velocity is positive
p(grow)=p(grow).*(ov(grow)./v(grow)); %homogenise the dislocation density over the growing grain (ratio of old/new volume)
%4. SYSTEM VOLUME AND AVERAGE GRAIN SIZE
sumv=sum(v(a)); %sum of volume
sumgrain=sum(d(a)); %sum of the grain sizes
avggrain=sumgrain/n; %average grain size (max diameter summation)
expgrain=(pi/4)*avggrain; %average grain size (for experimental comparison)
%5. VOLUME FRACTION OF GRAINS
for i=1:max(gen); %find the highest DRX generation cycle
cyc=find(gen==i); %find grain numbers of same generation
genvol(step,i)=sum(v(cyc)); %sum up generation volume and record in each strain step
end
drxfrac=sum(v(ig+1:all))/sumv; %DRX fraction: volume of all nucleated grains summed up over the system's volume