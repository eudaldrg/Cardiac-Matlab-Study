%% SHAPE, MOTION AND DEFORMATION ANALYSIS OF 2D ECHOCARDIOGRAPHIC SEQUENCES
% Drag your data into the workspace to begin the analysis. 
clc;
close all; 
warning('off', 'MATLAB:singularMatrix')
axis_parameters = [-60, 40, -120, 0];
%% Myocardial shape at the beginning of the cycle, at end of the systole
%% and at end of the diastole for all the subjects
PlotMyocardialShape(Subject, endSystole, endDiastole, axis_parameters)

%% Total length of the entire myocardial shape at end-systole (les), at
%% end-diastole (led) and ratio. 
number_of_subjects = length(Subject);
les = zeros(1, number_of_subjects);
led = zeros(1, number_of_subjects);
for i = 1:number_of_subjects
    les(i) = LengthHeartAtTime(Subject{i}, endSystole); % length myocardium at end systole
    led(i) = LengthHeartAtTime(Subject{i}, endDiastole); %length myocardium at end diastole
end
ratio = (led - les) ./ led;

%% Radial and longitudinal direction of the myocardial wall for all the
%% subjects
e_long = cell(number_of_subjects,1);
e_radial = cell(number_of_subjects,1);
for i = 1:number_of_subjects
    [e_long{i}, e_radial{i}] = GetLongAndRadDirectiorVectors(Subject{i});
end
% Plot direction vectors of the 4 Subjects 
PlotDirectionVectors(Subject, e_long, e_radial, axis_parameters)

%% Radial and Longitudinal Displacement Calculation 

U_rl = cell(number_of_subjects,1);
for i = 1:number_of_subjects
    U_rl{i} = RadialAndLongitudinalDisplacementAtTime(Subject{i}, endDiastole, e_long{i}, e_radial{i});
end

% RADIAL DISPLACEMENT OF 3 POINTS + LONGITUDINAL DISPLACEMENT OF 3 POINTS
figure;
subplot(2,1,1)
hold on;
axis([1,83,-5,16])
plot(1:83,U_rl{1}(:,6,1));
plot(1:83,U_rl{1}(:,19,1),'r');
plot(1:83,U_rl{1}(:,31,1),'m');
hold off;
title('Radial displacement at three different points of Subject 1')
legend('Point i=6','Point i=19','Point i=31')
xlabel('Number of frame')
ylabel('Displacement')
subplot(2,1,2)
hold on;
axis([1,83,-5,16])
plot(1:83,U_rl{1}(:,6,2));
plot(1:83,U_rl{1}(:,19,2),'r');
plot(1:83,U_rl{1}(:,31,2),'m');
hold off;
title('Longitudinal displacement at three different points of Subject 1')
legend('Point i=6','Point i=19','Point i=31')
xlabel('Number of frame')
ylabel('Displacement')

% Radial and longitudinal displacement for the 4 Subjects at mid-septal
% level
figure;
subplot(2,1,1)
hold on;
axis([1,83,-5,10])
plot(1:83,U_rl{1}(:,19,1));
plot(1:83,U_rl{2}(:,20,1),'r');
plot(1:83,U_rl{3}(:,18,1),'m');
plot(1:83,U_rl{4}(:,9,1),'g');
title('Radial displacement at the mid-septal level')
legend('Subject 1','Subject 2','Subject 3','Subject 4')
xlabel('Number of frame')
ylabel('Displacement')

subplot(2,1,2)
hold on;
axis equal;
axis([1,83,-5,10])
plot(1:83,U_rl{1}(:,19,2));
plot(1:83,U_rl{2}(:,20,2),'r');
plot(1:83,U_rl{3}(:,18,2),'m');
plot(1:83,U_rl{4}(:,9,2),'g');
title('Longitudinal displacement at the mid-septal level')
legend('Subject 1','Subject 2','Subject 3','Subject 4')
xlabel('Number of frame')
ylabel('Displacement')

% Imagesc plots
ImagescPlots(U_rl, [-2,2], [-5,5])

%%  Radial and Longitudinal Velocity Calculation

V_rl = cell(number_of_subjects, 1);
for i = 1:number_of_subjects
    number_of_points = length(Subject{i}.phi_x(1,:));
    number_of_measures = length(Subject{i}.phi_x(:, 1));
    V_rl{i} = zeros(number_of_measures - 1, number_of_points, 2);
    for j = 1:number_of_points
        V_rl{i}(:,j,1) = diff(U_rl{i}(:,j,1)) / timeInterval;
        V_rl{i}(:,j,2) = diff(U_rl{i}(:,j,2)) / timeInterval;
    end
end

% Radial velocity plot (spatiotemporal information using meshgrid)
MeshgridPlots(V_rl)

% Spatiotemporal radial & longitudinal velocity Subject of 1 (plotshaded)
x= 1:82;
y = zeros(2,82); % radial 
z = zeros(2,82); % longitudinal
for j = 1:82
        y(1,j) = (max(V_rl{1}(j,2:74,1)))';
        y(2,j) = (min(V_rl{1}(j,2:74,1)))';
        z(1,j) = (max(V_rl{1}(j,2:74,2)))';
        z(2,j) = (min(V_rl{1}(j,2:74,2)))';
end 
T_V_rl1 = zeros(82,2);
for i = 1:82
        T_V_rl1(i,1) = sum(V_rl{1}(i,2:74,1))/75; 
        T_V_rl1(i,2) = sum(V_rl{1}(i,2:74,2))/75;
end 
figure 
hold on;
axis([1,82,-250,250])
plotshaded(x,y,'r');
plotshaded(x,z,'b');
legend('Radial velocity','Longitudinal velocity')
plot(1:82,T_V_rl1(:,1),'r');
plot(1:82,T_V_rl1(:,2),'b');
title('Spatiotemporal representation of velocity Subject 1')
xlabel('Number of frame')
ylabel('velocity')

% Radial and longitudinal velocity plots (imagesc)
RadialLongitudinalVelocitiesImagescPlot(V_rl, [-50 50])

%% Strain Calculation 
E = cell(number_of_subjects, 1);
for i = 1:number_of_subjects
    number_of_measures = length(Subject{i}.phi_x(:, 1));
    number_of_points = length(Subject{i}.phi_x(1, :));
    E{i} = zeros(number_of_measures, number_of_points);
    for j = 1:number_of_measures
        for k = 2:(number_of_points - 1)
            diff_long_x = (Subject{i}.phi_x(j, k + 1) - Subject{i}.phi_x(j, k - 1)) / 2;
            diff_long_y = (Subject{i}.phi_y(j, k + 1) - Subject{i}.phi_y(j, k - 1)) / 2;
            diff_long_x_inicial = (Subject{i}.phi_x(1, k + 1) - Subject{i}.phi_x(1, k - 1)) / 2;
            diff_long_y_inicial = (Subject{i}.phi_y(1, k + 1) - Subject{i}.phi_y(1, k - 1)) / 2;
            diff_long = [diff_long_x, diff_long_y];
            diff_long_inicial = [diff_long_x_inicial, diff_long_y_inicial];
            E{i}(j, k) = (norm(diff_long) - norm(diff_long_inicial)) / norm(diff_long_inicial);
        end
    end
end

% Visualize the strain along the whole cycle for point i=19 for
% Subject 1, point i=20 for Subject 2, point i=18 for Subject 3 and point
% i=9 for Subject 4
figure;
hold on;
axis([1,83,-0.25,0.05]);
plot(1:83,E{1}(:,19));
plot(1:83,E{2}(:,20),'r');
plot(1:83,E{3}(:,18),'m');
plot(1:83,E{4}(:,9),'g');
title('Strain at the mid-septal level of all Subjects')
legend('Subject 1','Subject 2','Subject 3','Subject 4')
xlabel('Number of frame')
ylabel('Strain')

% Spatiotemporal strain plot of the 4 Subjects using imagesc function
StrainImagescPlots(E, [-0.2, 0.2])