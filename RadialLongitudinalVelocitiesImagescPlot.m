function [ ] = RadialLongitudinalVelocitiesImagescPlot(V_rl, axis_settings)
% RadialLongitudinalVelocitiesImagescPlot plots the radial and longitudinal
% velocities for all the subjects in a imagesc plot
figure
number_of_subjects = length(V_rl);
small_subplot_factor = ceil(sqrt(number_of_subjects / 2)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    h = subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i - 1);
    imagesc(V_rl{i}(:,:,1)')
    axis manual
    axis square
    set(h,'CLim', axis_settings)
    xlabel('time frames')
    ylabel('points')
    title(['Radial velocity Subject: ', num2str(i)])
    h = subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i);
    imagesc(V_rl{i}(:,:,2)') 
    axis manual
    axis square
    xlabel('time frames')
    ylabel('points')
    title(['Longitudinal velocity Subject: ', num2str(i)])
    set(h,'CLim', axis_settings)
end
end