function [x,y] = Circle(centerX,centerY, r)
% This function plots a circle with
% r radius
% centerX, centerY center cordinates
t=-pi:0.001:pi;
x=r*cos(t)+centerX;
y=r*sin(t)+centerY;

plot(x,y,'k:')
end