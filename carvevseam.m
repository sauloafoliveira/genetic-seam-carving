function [out] = carvevseam(img, seams, forceDouble) 
    if forceDouble
        img = double(img);
    end
    
    for i = 1:size(seams, 1)
        seam = seams(i, :);
        
        for j = 1:size(seam, 2)
            img(j, seam(j), :) = -1;
        end
       
    end
    
    s = size(img);
    s(2) = s(2) - size(seams, 1);
    out = zeros(s);
    for j = 1:size(img, 1)
        out(j, :, :) = img(j, img(j, :, 1) ~= -1, :);
    end

    if forceDouble
        out = uint8(out);
    end
end