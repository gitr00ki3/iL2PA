function [error, rms] = myrms(anchor, newAnchor, size)
    error = 0;
    rms = zeros(size, 1);
    for i=1:size
        rms(i) = (anchor(i,1) - newAnchor(i,1))^2 + (anchor(i,2) - newAnchor(i,2))^2;
        error = error +  rms(i);
    end;
    error = sqrt( error/size);
end