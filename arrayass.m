%DRX MODEL: written by Cram DG (2012)
%---arrayass.m---
%1. CREATING PROPERTY ARRAYS FOR EACH GRAIN
a=zeros(1,ig); %grain numbers that are 'live'
b=zeros(1,100000); %parameter which states if the grain is live or dead (1=live, 0=dead)
cplx_curve=zeros(maxcomplexions,ig); %complextion radii of grain
cplx_radius=zeros(maxcomplexions,ig); %radius of area of complextion patches
cplx_nuc=zeros(maxcomplexions,ig); %array of complexions being nucleated during a strain step
cplx_area=zeros(1,100000); %surface area of complexion patch
c=zeros(maxcomplexions,ig); %array of live complexions (0=dead)
d=zeros(1,100000); %diameter
v=zeros(1,100000); %volume
m=zeros(1,100000); %Taylor factor
p=zeros(1,100000); %dislocation density
f=zeros(1,100000); %flow stress
e=zeros(1,100000); %total strain
de=zeros(1,100000); %strain increment
r=zeros(1,100000); %strain-rate
g=zeros(1,100000); %g parameter
s=zeros(1,100000); %saturation stress
y=zeros(1,100000); %yield stress
pl=zeros(1,100000); %plastic stress
dpl=zeros(1,100000); %plastic stress increment
sg=zeros(1,100000); %subgrain size
stor=zeros(1,100000); %stored energy
rc=zeros(1,100000); %critical subgrain radius
nc=zeros(1,100000); %normalized critical subgrain radius
fsg=zeros(1,100000); %fraction of subgrains that can nucleate (F_sub)
nuc=zeros(1,100000); %number of subgrains that can nucleate
nn=zeros(1,100000); %number of subgrains that have nucleated from this grain
vel=zeros(1,100000); %velocity of grain boundary
nosub=zeros(1,100000); %number of subgrains on the surface of grain
fsglost=zeros(1,100000); %fraction of large subgrains lost during nucleation (F_lost)
fsgcap=zeros(1,100000); %fraction of subgrains large enough to nucleate (exhausted subgrains truncated) (F_sub-F_lost)
trun=zeros(1,100000); %fraction of subgrains bigger than the grain size (F_trun)
fsgrcap=zeros(1,100000); %fraction of subgrains large enough to nucleate (non-physical subgrains truncated) (F_sub-F_trun)
fsgcapsum=zeros(1,100000); %fraction of subgrains big enough and available for nucleation, (either equal to fsgcap or fsgrcap), never goes below 0
sgplinc=zeros(1,100000); %subgrain size increment: stress dependence
sgpdec=zeros(1,100000);
ov=zeros(1,100000); %previous step volume
sgo=zeros(1,100000); %previous step subgrain size: stress dependence
sgtempinc=zeros(1,100000); %subgrain size increment: temperature dependence
gen=ones(1,100000); %DRX generation (cycle)
jlat=zeros(1,100000); %lattice flux
jexc=zeros(1,100000); %excess flux
vellat=zeros(1,100000); %subgrain velocity from lattice flux
velexc=zeros(1,100000); %subgrain velocity from excess flux
vac=zeros(1,100000); %vacancy concentration (equilibrium vacancies + excess vacancies)
%result/data storage arrays
res=NaN(10000,125); %results array, for recording at every strain step
genvol=zeros(100000,70); %volume of different grain generations
%2. ASSIGNING INITIAL VALUES TO ARRAYS
for i=1:ig
a(i)=i; %each live grain assigned a number
b(i)=1; %each initial grain is live, assigned 1
end
initgrainsize=initgrainsize*(4/pi); %converts the initial grain size from the maximum diameter, to the experimentally determined average grain diameter
avggrain=initgrainsize; %average grain size is equal to the initial grain size (necessary for initial k in isowork)
for i=1:ig
%grain size and orientation
d(i)=normrnd(initgrainsize,initgrainsize/3); %grain diameters distributed with a normal distribution: mean and standard deviation (D/3)
if d(i)<0 %if a negative diameter is assigned, it is replaced with one about the mean, standard deviation (D/8)
d(i)=normrnd(initgrainsize,initgrainsize/8);
end
d(i)=d(i)*1e-6; %diameter converted to microns
dinit(i)=d(i); %initial grain diameter saved
m(i)=normrnd(avgtaylor,taydist); %Taylor factor with mean and standard deviation
if (m(i)<2.3)||(m(i)>4.2) %if the Taylor factor lies outside the limits (2.3-4.2), then reassign
m(i)=normrnd(avgtaylor,taydist*0.75); %reassign with a new Taylor factor
end
%To ensure DRX systems with small number of initial grains has a
%relatively average and consistent representation for comparison, the first 7 grains are reassigned
d(1:7)=[80e-6*(4/pi),80e-6*(4/pi),80e-6*(4/pi),80e-6*(4/pi),80e-6*(4/pi),80e-6*(4/pi),80e-6*(4/pi)];
m(1:7)=[3,2.8,2.9,3.12,3.22,3.32,3.06];
v(i)=(4/3)*pi()*((d(i)/2)^3); %volume (assumed spherical)
%stress and dislocation density
p(i)=10^(9.5); %low initial dislocation density
pl(i)=alpha(metal)*m(i)*shearmod*burger*sqrt(p(i)); %plastic stress from dislocation density
y(i)=(polyval(yield_poly(1:3,metal),temp))-4e6+(kyield*((d(i)*(pi/4)).^(-0.5)))*10^6; %yield stress
f(i)=y(i)+pl(i); %flow stress (sum of yield and plastic components)
%nucleation parameters
sg(i)=0.5*ksub(metal)*burger*(shearmod/pl(i))^(msub); %subgrain size
sgo(i)=sg(i); %subgrain size of previous step
vac(i)=vaceq; %initial vacancy concentration is the equilibrium concentration
end

%complextion size and distribution
complexions_initial=0; 
for i=1:ig
    avg_cplx_r = cplx_size_ratio*0.5*initgrainsize*1e-6; %average complextion radius
    for j=1:maxcomplexions
        cplx_curve(j,i) = d(i)/2; %array of complexion radii of curvature
        cplx_radius(j,i) = normrnd(avg_cplx_r, 0.05*avg_cplx_r); %array of complexion radii on surface of grain
        if sum(pi.*cplx_radius(:,i).^2) > max_cplx_area*pi*d(i)^2 %limits complexions to a specified fraction of grain surface
            cplx_radius(j,i)=0;
            break
        end
        complexions_initial = complexions_initial+1; %initial number of complexions in system
        c(j,i)=i; %assign complexion as live
    end
    cplx_area(i) = sum(pi*cplx_radius(:,i).^2);
end
total_cplx_nuc=0;
errorc=0; %complexions that tried to nucleate but were larger than grain
