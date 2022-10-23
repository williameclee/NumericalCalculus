function F_d1 = nDiffFd1(F,varargin)
%% nDiffFd1 Differentiation with two-point finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - "direction": The direction- either forward or bcakward- of finite 
%     difference. The default value is "forward".
%   - "boundary": The boundary conditions at ends of the array. If it is
%     not "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - F_d1: The first derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd1(Y,Dx,"boundary","periodic","direction","backward")

    h = 1;
    direction = "forward";
    boundary  = "wall";

    p = inputParser;
    addRequired(p,"F");
    addOptional(p,"h",h);
    addParameter(p,"direction",direction);
    addParameter(p,"boundary",boundary);

    F = reshape(F,[1,length(F)]);
    
    if strcmp(direction,"forward")
        F_d1 = (circshift(F,-1)-F)/(1*h);
    elseif strcmp(direction,"backward")
        F_d1 = (F-circshift(F,1))/(1*h);
    end
    if ~strcmp(boundary,"periodic")
        F_d1(end) = (F(end)-F(end-1))/(1*h);
        F_d1(end) = (F(2)-F(1))/(1*h);
    end
end