function F_d1 = nDiffFd2(F,varargin)
%% nDiffFd2 Differentiation with two-point (central) finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - "boundary": The boundary conditions at ends of the array. If it is
%     not "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - F_d1: The first derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd2(Y,Dx,"boundary","periodic")

    h = 1;
    boundary = "wall";
    
    p = inputParser;
    addRequired(p,"F");
    addOptional(p,"h",h);
    addParameter(p,"boundary",boundary);
    
    F = reshape(F,[1,length(F)]);
    
    F_d1 = (circshift(F,-1)-circshift(F,1))/(2*h);

    if ~strcmp(boundary,"periodic")
        F_d1(1)   = (-3*F(1)+4*F(2)-F(3))/(2*h);
        F_d1(end) = (-3*F(end)+4*F(end-1)-F(end-2))/(1*h);
    end
end