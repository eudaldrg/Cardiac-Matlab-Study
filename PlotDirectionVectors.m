function [ ] = PlotDirectionVectors( Subject, e_long, e_radial, axis_settings )
% PlotDirectionVectors it gets the Subject data and direction vectors and 
% plots each point with its initial direction vector.
figure
axis equal
axis(axis_settings);
number_of_subjects = length(Subject);
plots_per_dimension = ceil(sqrt(number_of_subjects)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    subplot(plots_per_dimension, plots_per_dimension, i);
    title(['Subject ', num2str(i), ': Directions of each point at the beginning of the cycle']);
    hold on
    quiver(Subject{i}.phi_x(1,:), Subject{i}.phi_y(i,:), e_long{i}(1,:,1), e_long{i}(1,:,2),'r')
    quiver(Subject{i}.phi_x(1,:), Subject{i}.phi_y(i,:), e_radial{i}(1,:,1), e_radial{i}(1,:,2),'b')
    hold off
    xlabel('x coordinate')
    ylabel('y coordinate')
end
legend('Longitudinal direction','Radial direction')
end

