function F_d1 = nDiffFd6(F,varargin)
%% nDiffFd6 Differentiation with six-point (central) finite difference
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
% Y_prime = nDiffFd6(Y,Dx,"boundary","periodic")

    h = 1;
    boundary = "wall";

    p = inputParser;
    addRequired(p,"F");
    addOptional(p,"h",h);
    addParameter(p,"boundary",boundary);

    F = reshape(F,[1,length(F)]);

    F_d1 = (+45*circshift(F,-1)-45*circshift(F,1) ...
             -9*circshift(F,-2) +9*circshift(F,2) ...
             +1*circshift(F,-3) -1*circshift(F,3))/(60*h);

    if ~strcmp(boundary,"periodic")
        F_d1(1:3)  = (-147*F(1:3)+360*F(2:4)-465*F(3:5)+400*F(4:6)-225*F(5:7)+72*F(6:8)+10*F(7:9))/(60*h);
        F_d1(end-2:end)  = (-147*F(end-2:end)+360*F(end-3:end-1)-465*F(end-4:end-2)+400*F(end-5:end-3)-225*F(end-6:end-4)+72*F(end-7:end-5)+10*F(end-8:end-6))/(60*h);
    end
end