function [ angle ] = hei2ang_servo( height )
% Converts to degrees the centemeters of the step high
b=[4.1234    0.0028   19.7756   -2.4426]; % calculado de sin_fit

% h(i)=b(1)*sin(2*pi*b(2)*(ang(i)-b(3)))+b(4)

angle = round(asin((height-b(4))/b(1))/(2*pi*b(2)) + b(3));

end
