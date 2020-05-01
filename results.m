%DRX MODEL: written by Cram DG (2012)
%---results.m---
%PRINTING RESULTS INTO ARRAY
%MACROSCOPIC
res(step,1)=time; %time
res(step,2)=strainmac; %macroscopic strain
res(step,3)=stressmac; %macroscopic stress
res(step,4)=expgrain; %average model grain size (to compare with experiment)
res(step,5)=drxfrac; %DRX volume fraction
res(step,6)=yieldmac; %macroscopic yield stress
res(step,7)=plasticmac; %macroscopic plastic stress
res(step,8)=taylormac; %macroscopic Taylor factor
res(step,9)=sumv; %volume of the system
res(step,10)=nucs; %nucleation potential sum
res(step,11)=avgdislosur; %grain boundary surface average dislocation density
res(step,12)=avgdislovol; %volume average dislocation density
%NUMBER OF GRAINS
res(step,13)=n; %number of grains in system
res(step,14)=new; %total number of nucleated grains that have entered the matrix
res(step,15)=tdes; %total number of grains that have ever been destroyed/eliminated
%INDIVIDUAL GRAIN PROPERTIES (grain 1)
%size and stress
res(step,16)=f(1); %flow stress
res(step,17)=e(1); %strain
res(step,18)=p(1); %dislocation density
res(step,19)=r(1); %strain-rate
res(step,20)=d(1); %diameter
res(step,21)=vel(1); %velocity
res(step,22)=y(1); %yield tress
res(step,23)=pl(1); %plastic stress
res(step,24)=s(1); %saturation stress
%nucleation/subgrain size evolution
res(step,25)=rc(1); %critical radius
res(step,26)=sg(1); %subgrain size (radius)
res(step,27)=sgplinc(1); %subgrain size (radius) change increment, 1/stress
res(step,28)=sgtempinc(1); %subgrain growth increment from flux
res(step,29)=nuc(1); %nucleation potential
res(step,30)=rc(1)/sg(1); %critical radius/subgrain size ratio
res(step,31)=stor(1); %stored energy
%subgrain size distribution fractions
res(step,32)=fsg(1); %fraction of subgrains that are big enough to nucleate: larger than rc (F_sub)
res(step,33)=fsglost(1); %fraction of large subgrains lost during nucleation (F_lost)
res(step,34)=trun(1); %fraction of subgrains bigger than the grain size (F_trun)
res(step,35)=fsgcap(1); %fraction of subgrains large enough to nucleate (exhausted subgrains truncated) (F_sub-F_lost)
res(step,36)=fsgrcap(1); %fraction of subgrains large enough to nucleate (non-physical subgrains truncated) (F_sub-F_trun)
res(step,37)=fsgcapsum(1); %fraction of subgrains big enough and available for nucleation, (either equal to fsgcap or fsgrcap), never goes below 0
%vacancies
res(step,38)=vac(1); %vacancy evolution
res(step,39)=jlat(1); %lattice vacancy flux
res(step,40)=jexc(1); %excess vacancy flux
res(step,41)=vellat(1); %lattice vacancy flux velocity
res(step,42)=velexc(1); %excess vacancy flux velocity
%complexion data
res(step,43) = total_cplx_nuc; %total number of grains nucleated from complexions
res(step,44) = errorc; %total number of complexions that tried to nucleate but were larger than the grain
res(step,45) = total_cplx_area;%total surface area of complexions in the system

