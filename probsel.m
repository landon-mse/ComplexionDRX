%DRX MODEL: written by Cram DG (2012)
%---probsel.m---
i=0; %reset
select=zeros(1,n); %create a matrix with n live grains
select=cumsum(nuc(a)); %cumulatively sum up the potential from each live grain with each consecutive cell
seln=rand()*sum(nuc(a)); %get random number between 0 and the nucleation potential sum over the live grains.
i=min(find(select>seln)); %find the 'live' grain number whose nucleation potential lies on this random number. The larger the grain's nulceation potential, the most chance it has being selected.
i=a(i); %select the actual grain number.
if (2*rc(i)>d(i)) || (nuc(i)==0) %precautionary, if a grain which shouldn't be selected is, the grain with the highest nulceation potential is instead.
i=find(nuc==max(nuc)); %find the grain with the highest nucleation potential
end