function [ ] = StrainImagescPlots(E, axis_settings)
% StrainImagescPlots plots the strain for all the
% subjects in an imagesc plot
figure
number_of_subjects = length(E);
subplot_factor = ceil(sqrt(number_of_subjects)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    h = subplot(subplot_factor , subplot_factor , i);
    imagesc(E{i}(:,:)')
    axis manual
    axis square
    set(h,'CLim', axis_settings)
    xlabel('time frames')
    ylabel('points')
    title(['Strain Subject: ', num2str(i)])
end
end