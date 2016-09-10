function [ ] = ImagescPlots(U_rl, radial_settings, longitudinal_settings)
% ImagescPlots plots the radial and longitudinal displacements for all the
% subjects in a imagesc plot
figure
number_of_subjects = length(U_rl);
small_subplot_factor = ceil(sqrt(number_of_subjects / 2)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    h = subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i - 1);
    imagesc(U_rl{i}(:,:,1)')
    axis manual
    axis square
    xlabel('time frames')
    ylabel('points')
    title(['Radial displacement Subject: ', num2str(i)])
    set(h,'CLim', radial_settings)
    h = subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i);
    imagesc(U_rl{i}(:,:,2)') 
    axis manual
    axis square
    xlabel('time frames')
    ylabel('points')
    title(['Longitudinal displacement Subject: ', num2str(i)])
    set(h,'CLim', longitudinal_settings)
end
end