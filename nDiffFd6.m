function Fd = nDiffFd6(F,varargin)
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
        Fd = (+ 45*circshift(F,-1) - 45*circshift(F,1) ...
              -  9*circshift(F,-2) +  9*circshift(F,2) ...
              +  1*circshift(F,-3) -  1*circshift(F,3))/(60*h);
        if ~strcmp(boundary,"periodic")
            Fd(1:3)       = (- 147*F(1:3) + 360*F(2:4) ...
                             - 450*F(3:5) + 400*F(4:6) ...
                             - 225*F(5:7) +  72*F(6:8) ...
                             - 10*F(7:9))/(60*h);
            Fd(end-2:end) = (+ 147*F(end-2:end)   - 360*F(end-3:end-1) ...
                             + 450*F(end-4:end-2) - 400*F(end-5:end-3) ...
                             + 225*F(end-6:end-4) -  72*F(end-7:end-5) ...
                             +  10*F(end-8:end-6))/(60*h);
        end
    elseif order == 2 % secend derivative
        Fd = (-490*F ...
              +270*circshift(F,-1)+270*circshift(F,1) ...
               -27*circshift(F,-2) -27*circshift(F,2) ...
                +2*circshift(F,-3)  +2*circshift(F,3))/(180*h^2);
        if ~strcmp(boundary,"periodic")
            Fd(1:3)       = (+  938*F(1:3) - 4014*F(2:4) ...
                             + 7911*F(3:5) - 9490*F(4:6) ...
                             + 7380*F(5:7) - 3618*F(6:8) ...
                             + 1019*F(7:9) -  126*F(8:10))/(180*h^2);
            Fd(end-2:end) = (+ 938*F(end-2:end)    - 4014*F(end-3:end-1) ...
                             + 7911*F(end-4:end-2) - 9490*F(end-5:end-3) ...
                             + 7380*F(end-6:end-4) - 3618*F(end-7:end-5) ...
                             + 1019*F(end-8:end-6) -  126*F(end-9:end-7))/(180*h^2);
        end
    end
end