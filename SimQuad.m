function [ X_end,Y_end ] = SimQuad( D,N )
% SimQuad calculates N number of lines with a length D on a sector
% quadrant y>0,x>1 and returns the endpoint coordinates as a output vector set.
% First line always is y = 0, x = D. When N => 2 the last line is y = D, x = 0.
% This is done so that the spacing between each line is always equal.
% 
% D     = Distance from origin. (m)
% N     = Number of endpoints.  (#)
% X_end = Vector of endpoint X-coordinates.
% Y_end = Vector of endpoint Y-coordinates.
% 
% Paxton Juuti
% TTY 24.07.2014
% 

if N>2,
    dTheta=90/(N-1); %The angle difference between lines
end

X_end(N)=0; %Initializing the endpoint vector
Y_end(N)=0;
X_end(1)=D; %First lines endpoint coordinates
Y_end(1)=0;

if N>=2,
    X_end(end)=0;  %Last lines endpoint coordinates if N~=1
    Y_end(end)=D;
    if N>2,
        for i=2:(N-1), %Loop for calculating the endpoints
            X_end(i)=D*cosd((i-1)*dTheta);
            Y_end(i)=D*sind((i-1)*dTheta);
        end
    end
end
end