%DRX MODEL: written by Cram DG (2012)
%---isowork.m---
%1. CHOOSING A K VALUE
%calculation of average dislocation density over grains
sumv=sum(v(a)); %sum of total volume
sumpv=sum(v(a).*p(a)); %sum of total volume*dislocations
avgdislovol=sumpv/sumv; %volume average dislocation density
k=approxstraininc*(alpha(metal)*3*shearmod*burger*((avgdislovol)^(0.5))+((polyval(yield_poly(1:3,metal),temp))-4e6+(kyield*((avggrain*(pi/4)).^(-0.5)))*10^6)); %k calculation which that ensures the correct strain increment
%2. UPDATING STRAIN IN EACH GRAIN AND THE MACROSCOPIC STRAIN
de(a)=k./f(a); %isowork to determine strain increment in each grain
e(a)=e(a)+de(a); %update strain in each grain
sumve=sum(v(a).*de(a)); %sum of volume*strain increment
strainmacinc=sumve/sumv; %macroscopic strain increment
strainmac=strainmac+strainmacinc; %update macroscopic strain
%3. TIMESTEP
timeinc=strainmacinc/extsrate; %time increment
vactimeinc=timeinc/vactime; %time increment for vacancy calculation
time=time+timeinc; %time
%4. CONSTITUITIVE LAW: UPDATE STRESS IN EACH GRAIN
r(a)=de(a)./timeinc; %strain-rate in each grain
g(a)=((boltz*temp)./(shearmod*((burger)^(3))))*(log(1e7./r(a))); %'g' parameter
s(a)=polyval(sat_lin(1:2,metal),g(a)); %linear fit from mastercurve: log(sat/u)
s(a)=(10.^(s(a))*shearmod); %saturation stress in each grain (yield stress not inc.)
dpl(a)=(initshr(metal)*(1-(pl(a)./s(a)))).*de(a); %update plastic stress increment (Voce law)
y(a)=(polyval(yield_poly(1:3,metal),temp))-4e6+(kyield*((d(a)*(pi/4)).^(-0.5)))*10^6; %yield stress
pl(a)=pl(a)+dpl(a); %update plastic stress
f(a)=y(a)+pl(a); %flow stress
p(a)=((pl(a))./(alpha(metal)*m(a)*shearmod*burger.*((r(a)./srateeff).^(sratesens)))).^2; %update disloation density
sumvpl=sum(v(a).*pl(a)); %volume*plastic stress sum
sumvm=sum(v(a).*m(a)); %volume*taylor factor sum
plasticmac=sumvpl/sumv; %macroscopic plastic stress (volume average)
yieldmac=polyval(yield_poly(1:3,metal),temp); %using experimental value for macroscopic yield stress
taylormac=sumvm/sumv; %average Taylor factor (volume average)
stressmac=yieldmac+plasticmac; %macroscopic stress