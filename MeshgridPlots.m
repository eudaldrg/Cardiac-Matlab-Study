function [ ] = MeshgridPlots(V_rl)
% MeshgridPlots plots the radial and longitudinal velocities for all the
% subjects in a meshgrid plot
figure
number_of_subjects = length(V_rl);
small_subplot_factor = ceil(sqrt(number_of_subjects / 2)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    time = 1:length(V_rl{i}(:,1,1));
    points = 1:length(V_rl{i}(1,:,1));
    subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i - 1);
    meshgrid(time, points);
    surf(V_rl{i}(:,:,1));
    axis([1 length(points) 1 length(time) -150 150]);
    title(['Radial velocity Subject: ', num2str(i)])
    xlabel('points')
    ylabel('time frames')
    zlabel('velocity')
    subplot(small_subplot_factor, small_subplot_factor * 2, 2 * i);
    meshgrid(time, points);
    surf(V_rl{i}(:,:,2));
    axis([1 length(points) 1 length(time) -150 150]);
    title(['Longitudinal velocity Subject: ', num2str(i)])
    xlabel('points')
    ylabel('time frames')
    zlabel('velocity')
end
end