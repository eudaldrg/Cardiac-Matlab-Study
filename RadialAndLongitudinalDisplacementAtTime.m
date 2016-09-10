function [ result ] = RadialAndLongitudinalDisplacementAtTime(subject, time, e_long, e_radial)
% RadialAndLongitudinalDisplacementAtTime computes the radial and
% longitudinal displacement vectors at time "time" for th subject "subject"
% given the correspondent direction vectors

number_of_points = length(subject.phi_x(1,:));
result = zeros(time,number_of_points,2);
for j = 2:time
    for i = 1:number_of_points
    tmp = [subject.phi_x(j,i) - subject.phi_x(1,i); subject.phi_y(j,i) - subject.phi_y(1,i)];
    P1 = [e_radial(j,i,1),e_long(j,i,1);e_radial(j,i,2),e_long(j,i,2)];
    % Ax = b. One way to solve this is with x = inv(A)*b. 
    % A better way, from both an execution time and numerical accuracy 
    % standpoint, is to use the matrix division operator x = A\b.
    result(j,i,:) = inv(P1)*tmp;
    end
end

end

