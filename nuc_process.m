%3. NUCLEATION PROCESS
%adding a new grain into the system
while nucstep>0
    if start_cplx_nuc
        i=x(nucstep);
        myRadius = cplx_curve(q,i);
    else
        probsel %script that selects a grain to nucleate based on probability, proportional to each grain's nucleation potential
        myRadius = rc(i);
    end
    
    if (v(i)-(((4/3)*pi*(myRadius)^3)))<0 %if the subtraction of a subgrain will make the grain's volume negative in a previous step, nucleation won't occur, and probsel tries again to select a grain
        
        if start_cplx_nuc
            nucstep=nucstep-1;
            errorc=errorc+1;
            total_cplx_nuc = total_cplx_nuc-1;
            c(q,i)=i; %makes failed complexion availible to nucleate again
        else
            nuc(i)=0; %making this grain's nuc(i)=0, means this grain will no longer be selected for a nucleation event in this strain step
        end 
        
    else
        all=all+1; %add one to the number of grains that have existed
        arrayup %script which adds a new column to all existing grain property arrays due to new grain
        nn(i)=nn(i)+1; %add one to the number of subgrains that nucleated from parent grain
        n=n+1; %add one to the number of live grains
        if start_cplx_nuc==0
        new=new+1; %add one to the total number of grains that have nucleated
        end
        b(all)=1; %call new grain live
        a(n)=all; %updates the live grain array
        %assigning new grain's properties
        d(all)=2.*myRadius; %diameter equal to parent's critical subgrain size
        v(all)=(4/3)*pi()*((d(all)/2)^3); %volume
        gen(all)=gen(i)+1; %DRX grain generation 1 greater than parent
        %assigning taylor factor
        m(all)=normrnd(avgtaylor,taydist); %orientation selection by normal distribution
        if (m(all)<2.3)||(m(all)>4.2) %if the Taylor factor lies outside the limits (2.3-4.2), reassign
            m(all)=normrnd(avgtaylor,taydist*0.75);
        end
        %assigning stress properties
        p(all)=1e9; %the grain is low in dislocations (not completely zero for physical and calculation reasons)
        pl(all)=alpha(metal)*m(all)*shearmod*burger*sqrt(p(all)); %approximate value for low plastic stress
        y(all)=(polyval(yield_poly(1:3,metal),temp))-4e6+(kyield*((d(all)*(pi/4)).^(-0.5)))*10^6; %yield stress
        f(all)=y(all)+pl(all); %flow stress (plastic + yield)
        %assigning subgrain size
        sg(all)=0.5*ksub(metal)*burger*(shearmod/pl(all))^(msub); %starting subgrain size
        sgo(all)=0.5*ksub(metal)*burger*(shearmod/pl(all))^(msub); %starting subgrain size (saved for the step before)
        cplx_area(all)=0;
        
        %properties of the old (parent) grain that has nucleated
        v(i)=v(i)-v(all); %take away the subgrain volume from the grain
        d(i)=2*(((3*v(i))/(4*pi))^(1/3)); %assigns new diameter
        if start_cplx_nuc
            %cplx_area(all) = pi*myRadius^2;
            cplx_area(i) = cplx_area(i) - pi*myRadius^2;
        end
        if rc(i)>2*d(i) %if grain's critical subgrain size is larger than its grain size, set nucleation potential to 0 so it can't nucleate again in this strain step.
            nuc(i)=0;
        end
        %recalculate truncation for all grains
        if nucstep==1 %truncation calculations occur just after the final nucleation event for this strain step
            for j=1:n
                i=a(j); %selecting only the live grain numbers
                if fsgcapsum(i)>0 %if the grain was involved in the nucleation potential sum
                    fsglost(i)=((floor(nucs)/nucs)*(fsgcapsum(i)))+fsglost(i); %updating the new truncated area of each subgrain size distribution which was involved in the nucleation event in all grains
                    fsgcap(i)=fsg(i)-fsglost(i); %new fraction of subgrains available to nucleate
                end
            end
         end
         nucstep=nucstep-1; %one less nucleation event that needs to occur in the strain step
    end
end