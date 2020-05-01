%DRX MODEL: written by Cram DG (2012)
%---nucleation.m---
%1. SUBGRAIN SIZE
%stress dependence
sgplinc(a)=0.5*ksub(metal)*burger.*(shearmod./pl(a)).^(msub)-sgo(a); %subgrain size (radius) change increment, 1/stress
%sgpdec(a) = -0.5*ksub(metal).*burger.*(shearmod./pl(a).^2).*initshr(metal).*(1-pl(a)./s(a)).*e(a).*timeinc;
sgo(a)=0.5*ksub(metal)*burger*(shearmod./pl(a)).^(msub); %saving old subgrain size (radius)

% vacancy evolution
for i=1:vactime
vac(a)=vac(a)+(((0.1*pl(a).*atvol)./(Qf_eff*1.6e-19)).*r(a).*vactimeinc)-((vacdiff./((2*sg(a)).^2)).*(vac(a)-vaceq)).*vactimeinc; %concentration of vacancies
end
jlat(a)=(2*latdiff*subdeg.*pl(a))./(burger*boltz*temp.*m(a)); %lattice vacancy flux
jexc(a)=(vacdiff.*(vac(a)-vaceq))./(atvol.*sg(a)); %excess vacancy flux
vellat(a)=A_fac*jlat(a)*(burger^2)*subwidth; %lattice vacancy flux velocity
velexc(a)=A_fac*jexc(a)*(burger^2)*subwidth; %excess vacancy flux velocity
% temperature dependence
sgtempinc(a)=(vellat(a)+velexc(a)).*timeinc; %subgrain growth increment from flux

%subgrain size evolution
sg(a)=sg(a)+sgtempinc(a)+sgplinc(a); %subgrain size updated: addition of temperature and stress increments
%sg(a) = sg(a)+sgpdec(a)+sgplinc(a);
%COMPLEXION RADIUS EVOLUTION

%2. NUCLEATION POTENTIAL
stor(a)=0.5*p(a)*shearmod*burger*burger; %driving force (stored energy)
rc(a)=(2*interf)./(stor(a)); %critical radius
nc(a)= rc(a)./sg(a); %normalized critical radius
fsg(a)=exp((-pi*(nc(a)).^2)./4); %fraction of subgrains that are big enough to nucleate: larger than rc (F_sub)
fsgcap(a)=fsg(a)-fsglost(a); %fraction of subgrains large enough to nucleate, but truncated due to the loss of large subgrains from prior nucleation (F_sub-F_lost)
trun(a)=exp((-pi*((d(a)./2)./sg(a)).^2)./4); %truncated fraction of subgrains bigger than the grain size (F_trun)
fsgrcap(a)=fsg(a)-trun(a); %fraction of subgrains large enough to nucleate, but truncated due to being larger than grain size (F_sub-F_trun)
nosub(a)=((64*(d(a)./2).^2)./(pi*pi*(sg(a).^2))); %number of subgrains on surface
%fraction of subgrains available for nucleation summation, (set to either fsgcap or fsgrcap)
for j=1:n
i=a(j);
if fsglost(i)>trun(i)
fsgcapsum(i)=fsgcap(i);
else
fsgcapsum(i)=fsgrcap(i);
end
end
%nucleation potential in each grain
fsgcapsum(find(fsgcapsum<0))=0; %if a grain has a negative fraction of subgrains available, this is unphysical and hence is set to 0
nuc(a)=fsgcapsum(a).*nosub(a); %number of subgrains that nucleate in each grain (fraction*subgrains that lie on grain boundary surface)
%summing up nucleation potential
nucs=sum(nuc(a)); %summing the nucleation potential over all grains
onucstep=floor(nucs); %number of subgrains to be nucleated (rounded down)

sumpdsq=sum((p(a).*((d(a)./2).^2))); %sum of dislocations*r^2
sumdsq=sum((d(a)./2).^2); %sum of r^2
avgdislosur=sumpdsq/sumdsq; %surface area average dislocation density

%Evolve complexion curvature and select complexions to nucleate
numc=0;
start_cplx_nuc=0;
%Loop iterates through each "complexion number" for all grains
if UseCompNuc
    cmplx_nuc=zeros(maxcomplexions,ig);
    for q = 1:maxcomplexions
        does_cplx_nuc = false(1,ig);
        x = find(c(q,:)>0); %indices of active complextions
        cplx_curve(q,x) = cplx_curve(q,x) + (mco*(shearmod*burger^2*abs(avgdislosur-p(x))-2*interf_comp./cplx_curve(q,x)))*timeinc;
        does_cplx_nuc(x) = cplx_curve(q,x) < cplx_radius(q,x);
        cmplx_nuc(q,does_cplx_nuc) = 1; %array of nucleated complexions
        c(q,does_cplx_nuc) = 0; %eliminates nucleated complexions
        numc = numc + length(find(does_cplx_nuc==1)); %number of nucleations from complexions this strain step
    end
    if numc>0
        start_cplx_nuc=1;
    end
    total_cplx_nuc = (total_cplx_nuc + numc); %tracks total number of complexions nucleated
end

%nucleation of complexions
if start_cplx_nuc
    for q=1:maxcomplexions %goes into nucleation process 
        x=find(cmplx_nuc(q,:)>0);
        nucstep = length(x);
        nuc_process;
    end
end
start_cplx_nuc=0;

%Nucleation of subgrains
nucstep = onucstep;
nuc_process


    