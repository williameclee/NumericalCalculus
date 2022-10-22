function F_d1 = nDiffFd4(F,varargin)
    h = 1;
    boundary = "wall";

    p = inputParser;
    addRequired(p,"F");
    addOptional(p,"h",h);
    addParameter(p,"boundary",boundary);

    F = reshape(F,[1,length(F)]);

    F_d1 = (+8*circshift(F,-1)-8*circshift(F,1) ...
            -1*circshift(F,-2)+1*circshift(F,2))/(12*h);

    if ~strcmp(boundary,"periodic")
        F_d1(1:2)        = (-25*F(1:2)+60*F(2:3)-36*F(3:4)+16*F(4:5)-3*F(5:6))/(12*h);
        F_d1(end-1:end)  = (-25*F(end-1:end)+60*F(end-2:end-1)-36*F(end-3:end-2)+16*F(end-4:end-3)-3*F(end-5:end-4))/(12*h);
    end
end