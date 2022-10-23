function Fd = nDiffFd2(F,varargin)
%% nDiffFd2 Differentiation with two-point (central) finite difference
% INPUT:
%   - F: The 1D array to be differentiated.
%   - h: The uniform spacing between each element. The default value is 1.
%   - order: The order of derivative. The default value is 1 (first 
%     derivative).
%   - "boundary": The boundary conditions at ends of the array. If it is
%     not "periodic", forward and backward differences will be utilised to
%     achieve the same order of truincation error. The default value is
%     "wall".
% OUTPUT:
%   - Fd: The derivative of the input array.
% EXAMPLE:
% Y_prime = nDiffFd2(Y,Dx,"order",2,"boundary","periodic")

    h     = 1;
    order = 1;
    boundary = "wall";
    
    p = inputParser;
        addRequired(p,'F');
        addOptional(p,'h',h);
        addOptional(p,'order',order);
        addParameter(p,'boundary',boundary);
        parse(p,F,varargin{:});
    F        = p.Results.F;
    h        = p.Results.h;
    order    = p.Results.order;
    boundary = p.Results.boundary;
    
    F = reshape(F,[1,length(F)]);
    
    if order == 1 % first derivative
        Fd = (circshift(F,-1)-circshift(F,1))/(2*h);
        if ~strcmp(boundary,"periodic")
            Fd(1)   = (-3*F(1)+4*F(2)-F(3))/(2*h);
            Fd(end) = (3*F(end)-4*F(end-1)+F(end-2))/(2*h);
        end
    elseif order == 2 % secend derivative
        Fd = (circshift(F,-1)-2*F+circshift(F,1))/(1*h^2);
        if ~strcmp(boundary,"periodic")
            Fd(1)   = (F(1)-2*F(2)+F(3))/(1*h^2);
            Fd(end) = (F(end)-2*F(end-1)+F(end-2))/(1*h^2);
        end
    end
end