function [ ] = PlotMyocardialShape(Subject, endSystole, endDiastole, axis_settings)
%PlotMyocardialShape it gets the Subject data and creates a subplot
% for each subject. Subject is the data for all the subjects (included
% in the Data.m file copied into the workspace. axis_settings are used for
% the plot. One should pass it as would normally do into the axis function.
figure
axis equal
axis(axis_settings);
number_of_subjects = length(Subject);
plots_per_dimension = ceil(sqrt(number_of_subjects)); % Guarantees enough space to show all the subjects
for i = 1:number_of_subjects
    subplot(plots_per_dimension, plots_per_dimension, i);
    title(['Subject ', num2str(i)]);
    hold on
    plot(Subject{i}.phi_x(1,:),Subject{i}.phi_y(1,:),'Marker','x');
    plot(Subject{i}.phi_x(endSystole,:),Subject{i}.phi_y(endSystole,:),'r');
    plot(Subject{i}.phi_x(endDiastole,:),Subject{i}.phi_y(endDiastole,:),'g');
    hold off
    xlabel('x coordinate')
    ylabel('y coordinate')
end
legend('t=1','End systole','End diastole');
end