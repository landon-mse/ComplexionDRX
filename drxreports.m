%%%Plots 
figure
plot(res(:,2),res(:,43)/complexions_initial)
title(['Initial # of Complexions: ', num2str(complexions_initial), ',   Strain Increment: ', num2str(approxstraininc, '%10.1e')])
ylabel('fraction of complexions nucleated')
xlabel('strain')

figure
plot(res(:,2),res(:,3))
title(['stress-strain  Initial # of grains: ', num2str(ig)])
ylabel('stress (MPa)')
xlabel('strain')
ax = gca;
ax.YAxis.Exponent = 6;
ylim([0 102e6])

% figure
% plot(res(:,2),res(:,44))
% ylabel('errors')
% xlabel('strain')
%  
% figure
% plot(res(:,2),res(:,45))
% ylabel('total complexion area')
% xlabel('strain')


figure
plot(res(:,2),res(:,14)./res(:,43))
title('ratio of number subgrain-DRXd grains to Complextion-DRXd grains')
ylabel('subgrain-DRXd / Complextion-DRXd grains')
xlabel('strain')

disp(['Mco: ', num2str(mco_fact), '   surface-energy_co: ', num2str(gammac_fact)]);