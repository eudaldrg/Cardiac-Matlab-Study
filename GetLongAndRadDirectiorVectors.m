function [ e_long, e_radial ] = GetLongAndRadDirectiorVectors(subject )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

number_of_measures = length(subject.phi_x(:, 1));
number_of_points = length(subject.phi_x(1, :));
e_long = zeros(number_of_measures, number_of_points, 2);
e_radial = zeros(number_of_measures, number_of_points, 2);
for j = 1:number_of_measures % time
    for i = 2:(number_of_points - 1) % points
        diff_long_x = subject.phi_x(j, i + 1) - subject.phi_x(j, i - 1);
        diff_long_y = subject.phi_y(j, i + 1) - subject.phi_y(j, i - 1);
        diff_long = [diff_long_x, diff_long_y];
        diff_long = diff_long / norm(diff_long);
        rotated_diff = [0, -1; 1, 0 ] * diff_long';
        if i < 75 / 2 
             e_long(j, i, :) = diff_long;
        else
            e_long(j, i, :) = -diff_long;
        end
        e_radial(j, i, :) = rotated_diff;
    end
    e_long(j, 1, :) = 0;
    e_radial(j, 1, :) = 0;
    e_long(j, number_of_points, :) = 0;
    e_radial(j, number_of_points, :) = 0;
end

end

