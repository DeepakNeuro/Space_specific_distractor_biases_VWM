function [degree] = degreeCal(length,distance,dimx,pixelx)
%DEGREECAL - It can calculate the degree for a given screen, distance and
%dimension of the object
%   all dimensions are in cm
%   lenght is number of pixels
%   dimx is size of screen in cm
%   pixelx is number of pixels of the screen in x direction

    lengthCm=length*dimx/pixelx;
    degree=2*atan(lengthCm/(2*distance));
    degree=degree*180/pi;

end

