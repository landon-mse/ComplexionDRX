%DRX MODEL: written by Cram DG (2012)
%Modified by Cordova LT (2020)

%---main.m---
clear %clear workspace
close all %closes old plots
modelrun=[0,725,0.002,80,75;0,775,0.001,80,50;0,775,0.002,80,100;0,775,0.027,80,25;0,775,0.2,80,25;0,875,0.002,80,1000;0,975,0.002,80,5000;0.2,725,0.002,80,7;0.2,775,0.001,80,7;0.2,775,0.002,80,7;0.2,775,0.027,80,7;0.2,775,0.2,80,7;0.2,875,0.002,80,96;0.2,975,0.002,80,384;2,725,0.002,80,7;2,775,0.001,80,7;2,775,0.002,80,7;2,775,0.027,80,7;2,775,0.2,80,7;2,875,0.002,80,24;2,975,0.002,80,96;5,725,0.002,80,7;5,775,0.001,80,7;5,775,0.002,80,7;5,775,0.027,80,7;5,775,0.2,80,7;5,875,0.002,80,24;5,975,0.002,80,96;0,775,0.002,45,500;0,775,0.002,150,50;0,775,0.002,300,40;];
%creates a 31x5 matrix which contains the 31 different model conditions run and presented in the model each column represents: Sn wt%, temperature, strain-rate, initial grain size, number of initial grains
%each row number represents the model number run (mr)
%1. pure Cu, 725, 0.002, 80, 75
%2. pure Cu, 775, 0.001, 80, 50
%3. pure Cu, 775, 0.002, 80, 100
%4. pure Cu, 775, 0.027, 80, 25
%5. pure Cu, 775, 0.2, 80, 25
%6. pure Cu, 875, 0.002, 80, 1000
%7. pure Cu, 975, 0.002, 80, 5000
%8. Cu-0.2Sn, 725, 0.002, 80, 7
%9. Cu-0.2Sn, 775, 0.001, 80, 7
%10. Cu-0.2Sn, 775, 0.002, 80, 7
%11. Cu-0.2Sn, 775, 0.027, 80, 7
%12. Cu-0.2Sn, 775, 0.2, 80, 7
%13. Cu-0.2Sn, 875, 0.002, 80, 96
%14. Cu-0.2Sn, 975, 0.002, 80, 384
%15. Cu-2Sn, 725, 0.002, 80, 7
%16. Cu-2Sn, 775, 0.001, 80, 7
%17. Cu-2Sn, 775, 0.002, 80, 7
%18. Cu-2Sn, 775, 0.027, 80, 7
%19. Cu-2Sn, 775, 0.2, 80, 7
%20. Cu-2Sn, 875, 0.002, 80, 24
%21. Cu-2Sn, 975, 0.002, 80, 96
%22. Cu-5Sn, 725, 0.002, 80, 7
%23. Cu-5Sn, 775, 0.001, 80, 7
%24. Cu-5Sn, 775, 0.002, 80, 7
%25. Cu-5Sn, 775, 0.027, 80, 7
%26. Cu-5Sn, 775, 0.2, 80, 7
%27. Cu-5Sn, 875, 0.002, 80, 24
%28. Cu-5Sn, 975, 0.002, 80, 96
%29. pure Cu, 775, 0.002, 45, 500
%30. pure Cu, 775, 0.002, 150, 50
%31. pure Cu, 775, 0.002, 300, 40
%1. SELECTING HOT-WORKING CONDITIONS
%a) PICK HOT-WORKING CONDITIONS AS USED IN THESIS:
mr=3; %selecting the model condition (pick value between 1 and 31)
%if batch run: replace mr=# with 'for mr=[1,2,3,etc...]'
wtsn=modelrun(mr,1); %weight percent Sn in Cu (limited to 0, 0.2, 2, 5)
temp=modelrun(mr,2); %temperature (kelvin)
extsrate=modelrun(mr,3); %applied external strain-rate (s^-1)
initgrainsize=modelrun(mr,4); %initial grain size (microns)
ig=modelrun(mr,5); %initial number of grains in system

UseCompGrow = 0; %Use complextion grain growth
UseCompNuc = 1; %Use complextion nucleation

%b) MANUALLY PICK HOT-WORKING CONDITIONS: (set manual_run==1)
manual_run=1; % (if manual_run=0, the value of mr is used to select the hot-working conditions, if manual_run=1, the hot-working conditions can be edited below)
if manual_run==1
    wtsn=0; %weight percent Sn in Cu (limited to 0, 0.2, 2, 5)
    temp=750; %temperature (kelvin)
    extsrate=2e-3; %applied strain-rate (s^-1)
    initgrainsize=75; %initial grain size (microns)
    ig=75; %initial number of grains in system
end

%2. ASSIGNING PARAMETERS AND ARRAYS
para %assigning parameters: materials properties
arrayass %creating grain arrays
disp(['HOT-WORKING CONDITIONS Sn%:',num2str(wtsn),' temp:',num2str(temp),' strain-rate:',num2str(extsrate),' initial grain size:',num2str(initgrainsize*(pi/4)),' initial number of grains:',num2str(ig)])

disp(['model output: ','STRAIN, ','FLOW STRESS, ','NUCLEATION POTENTIAL SUM, ','NUMBER OF GRAINS IN SYSTEM']); %key for the model's output results
startt=clock; %time at start of model run
dispt=clock; %time used for displaying output results in command window
%3. MODEL CALCULATION LOOP
while strainmac<strain %keep iterating until the final strain is reached
    step=step+1; %update number of strain steps
    isowork %polyphase plasticity model
    nucleation %nucleation model
    growth %growth model
    results %saves data in each strain step
    if ((etime(clock,dispt)>5)) %if display time>5sec, display the live model results
        disp([strainmac,plasticmac,nucs,n]); %display the live model results
        format short g %display numerical output format is short g
        dispt=clock; %reset display time
    end
end

%4. CRITICAL STRAIN/NUCLEATION RATE
crit=res(max(find(res(:,5)<0.005)),2); %critical strain when DRX fraction > 0.5%
nucavg=((all-ig)/((1-crit)*(1/extsrate)))/sumv; %model average nucleation rate dn/dt
nucr=zeros(step,1); %instantaneous nucleation rate matrix
pos=10; %a number that specifies the strains steps the instantaneous nucleation rate is averaged around
for i=pos+1:step-pos
    nucr(i,1)=((res(i+pos,14)-res(i-pos,14))/(res(i+pos,1)-res(i-pos,1)))/sumv; %instantaneous nucleation at every strain calculated by taking an average over 'pos'x2 strain steps around the strain
end

%5. MODEL TIME
modelrunt=etime(clock,startt); %model run time in sec
hrt=floor(modelrunt/3600); %model run time (hours) column
mint=floor((modelrunt-3600*hrt)/60); %model run time (min) column
sect=modelrunt-3600*hrt-60*mint; %model run time (sec) column
disp(['model runtime (hr,min,sec): ',num2str(hrt),', ',num2str(mint),', ',num2str(sect)]) %displays the mode run time in hours, minutes and seconds
%if batch run: save each results by this command and editing the folder: save(['C:\folder\',num2str(mr),'_',num2str(wtsn),'_',num2str(temp),'_',num2str(extsrate),'_',num2str(initgrainsize*(pi/4)),'.mat']) %Saves model's workspace and results
%add 'end' to complete for loop

%Plots and Reports
drxreports;
