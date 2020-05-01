%DRX MODEL: written by Cram DG (2012)
%---para.m---
%1. VARIABLES TO SET/CHANGE
strain=0.6; %strain which the model will iterate till
approxstraininc=1e-4; %average macroscopic strain increment
grainsizelimit=5E-7; %grain size limit, before grain is too small and eliminated from matrix
turnbull=(1/2); %Turnbull's factor for grain boundary mobility
binding_gb=15000/6.022e23; %binding energy between Sn and grain boundary
avgtaylor=3.06; %average Taylor factor
taydist=0.2; %standard deviation of Taylor factor
A_fac=75; %scaling factor for subgrain growth
%Sn dependent variables
%helps find Sn dependent parameters:
if wtsn==0
metal=1; %pure Cu
end
if wtsn==0.2
metal=2; %Cu-0.2Sn
end
if wtsn==2
metal=3; %Cu-2Sn
end
if wtsn==5
metal=4; %Cu-5Sn
end
initshr=[1.1e9,1.3e9,1.5e9,1.5e9]; %initial strain-hardening rate
alpha=[0.25,0.3,0.35,0.35]; %dislocation junction strength
ksub=[30,35,40,40]; %subgrain size parameter
%2. MODEL PARAMETERS (initial values)
n=ig; %number of grains in system
all=ig; %total number grains that have existed
des=0; %number of grains are being destroyed/eliminated in the current step
tdes=0; %total number of grains that have ever been destroyed/eliminated
new=0; %total number of nucleated grains that have entered the matrix
step=0; %main for loop step
time=0; %time
nucstep=0; %how many grains to nucleate in each step
strainmac=0; %macroscopic strain
plasticmac=0; %macroscopic stress
vactime=20; %factor which shortens the time for iteration when calculating the vacancy evolution
%3. MATERIAL CONSTANTS IN COPPER
tempm=1356; %melting temperature (K)
burger=2.55e-10; %Burger's vector
boltz=1.38065e-23; %Boltzmann's constant
sratesens=0.0222; %strain-rate sensitivity
srateeff=4.5e-7; %reference strain-rate
shearmod=42.1e9*(1+(((temp-300)/1356.)*-0.5)); %shear modulus
msub=1.; %subgrain size parameter
interf=0.625; %interfacial energy of grain boundary
kyield=0.04; %Hall-Petch constant for Cu in each individual grain (found from exp. at elevated temp)
gas=8.3145; %gas constant
molar=7.11e-6; %molar volume
atvol=1.18e-29; %atomic volume
gbdiffcoe=2.35e-5; %grain boundary diffusion coefficient
gbthick=1e-9; %grain boundary thickness
gbacteng=107.2e3; %grain boundary diffusion activation energy
subwidth=1e-9; %subgrain width
snacteng=195.4e3; %activation energy for Sn diffusing in Cu
sndiffcoe=8.2e-5; %pre-exponential diffusion coefficient for Sn diffusing in Cu
latdiffexp=1.6e-5; %pre-exponential to lattice diffusion
%4. VACANCIES AND BINDING ENERGIES - SOLUTE
atsn=((((1-(0.2/100))*118.71)./((wtsn./100)*63.546))+1).^(-1); %atomic fraction of Sn
Sn_dis_bind=0.38; %binding energy of Sn and dislocation (in eV)
Sn_vac_bind=0.38; %binding energy of Sn and vacancy (in eV)
Qsndis=Sn_dis_bind*96.485*1000; %binding energy of Sn and dislocation (in J/mol)
snconc=(atsn.*exp(Qsndis./(gas.*temp)))./(1+atsn.*exp(Qsndis./(gas.*temp))); %Sn concentration fraction around a dislocation
Qm=0.7; %migration activation energy (in eV)
Qf=1.28; %formation activation energy - pure Cu (in eV)
Eavg_esc=snconc*Sn_vac_bind; %average escape energy for vacancy from Sn (in eV)
Qbulk=(Qm+Qf)*96.485*1000; %bulk activation energy - pure Cu (J/mol)
Qbulk_eff=(Qm+Qf+Eavg_esc)*96.485*1000; %bulk activation energy - with Sn (J/mol)
Qvbind=Sn_vac_bind*96.485*1000; %binding energy of Sn and vacancy (in J/mol)
Qf_eff=Qf+Eavg_esc; %formation energy (effective) - with Sn (in eV) (used in excess flux)
Qm=(Qm)*96.485*1000; %migration activation energy (in J/mol)
Qf=(Qf)*96.485*1000; %formation activation energy - pure Cu (in J/mol)
vaceq=exp(-(Qf)/(gas*temp)); %equilibrium vacancy concentration
vacdiffexp=1.6e-5; %pre-exponential to vacancy diffusion
vacdiff=vacdiffexp*exp(-Qm/(gas*temp)); %vacancy migration diffusion
latdiff=latdiffexp*exp(-Qbulk_eff/(gas*temp)); %vacancy lattice diffusion
vac=vaceq; %vacancy concentration set to the equilibrium value
subdeg=3*(pi/180); %average subgrain misorientation
%5. GRAIN BOUNDARY CONSTANTS
gbdiff=gbdiffcoe*exp(-gbacteng/(gas*temp)); %grain boundary diffusion of Cu in Cu
mobpure=(turnbull)*((gbthick*gbdiff*molar)/(burger*burger*gas*temp)); %mobility of all grains (Turnbull's estimate)
sndiff=sndiffcoe*exp(-snacteng/(gas*temp)); %lattice diffusion of Sn in Cu
sndiff=sndiff*10; %diffusion across boundary ~x10 lattice diffusion
alphamob=((gbthick*(6.022e23/molar)*(boltz*temp)^2)/(binding_gb*sndiff))*(sinh(binding_gb/(boltz*temp))-(binding_gb/(boltz*temp))); %alpha term in Cahn's linear approx
mob=((1/mobpure)+alphamob*atsn)^(-1); %new mobility with solute drag
%6. EXPERIMENTAL DATA FILES
sat_lin=[-2.2713,-2.467,-2.477,-2.7779;-1.654,-1.3876,-1.2339,-1.1176;]; %2x4 matrix containing the linear fit coefficients for the saturation stress mastercurves
yield_poly=[104.236003413538,-42.0530490395261,-4.46383288908386,-121.579942767894;-199307.683939167,20308.1532197913,-45379.2528688980,88078.2606873771;98788228.3850782,33011992.9518360,84258615.3221549,86797716.7848746;]; %3x4 matrix containing the polyfit coefficients for the yield stress
sat=polyval(sat_lin(1:2,metal),((boltz*temp)/(shearmod*((burger)^(3))))*(log(1e7/extsrate))); %saturation stress calculation
sat=(10^(sat)*shearmod); %saturation stress
yield=polyval(yield_poly(1:3,metal),temp); %yield stress

%7. Complextion parameters
maxcomplexions=10; %maximum number of complexions that will be initiated per grain
max_cplx_area=0.75; %maximum fraction of grain surface area complextions can cover at intialization
cplx_size_ratio = 0.2; %ratio of average complextion radius to average grain radius
mco_fact = 0.5;
gammac_fact = 2;
mco = mco_fact*mob; %complextion mobility;
interf_comp = gammac_fact*interf; %interfacial energy of grain boundary complextion