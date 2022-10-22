function F_d1 = nDiffFd1(F,varargin)
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