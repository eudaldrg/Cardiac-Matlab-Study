function [ len ] = LengthHeartAtTime(subject, time)
%LengthHeartAtTime this function gets a single subject (for better
%usability) and a time and returns the length of the myocardial shape at 
%the specified time
len = 0;

for i = 2:(length(subject.phi_x(1,:)) - 1) % points
    diff_long_x = (subject.phi_x(time,i+1) - subject.phi_x(time,i-1))/2;
    diff_long_y = (subject.phi_y(time,i+1) - subject.phi_y(time,i-1))/2;
    diff_long = [diff_long_x,diff_long_y];
    len = len + norm(diff_long);
end 

end

%Another aproximation of longitud with slyghtly different results
%len = 0;
%for i = 2:74 %puntos
%   diff_long_x = (subject.phi_x(time,i) - subject.phi_x(time,i-1));
%   diff_long_y = (subject.phi_y(time,i) - subject.phi_y(time,i-1));
%   diff_long = sqrt( (diff_long_x)^2 + (diff_long_y)^2);
%   len = len + diff_long;
%end