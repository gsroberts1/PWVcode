--- Turn images indices into physical positions (in 3D space)
function[truePositions,rotationMatrix] = getTruePosition(info,sliceNum)
    if isfield(info,'ImageOrientationPatient') % if dicom
        originShift = [info.ImagePositionPatient;1]; % origin is top left corner of image
        xres = info.PixelSpacing(1);
        yres = info.PixelSpacing(2);
        zres = info.SliceThickness;
        matrixx = info.Width;
        matrixy = info.Height;
        
        % sometimes get extremely small values that should be 0, so round
        xVector = round(info.ImageOrientationPatient(1:3),8); % what direction rows run w/r/to x
        yVector = round(info.ImageOrientationPatient(4:6),8); % what direction the cols run w/r/to y
        zVector = [cross(xVector,yVector);0];
        
        xVector = [xVector;0];
        yVector = [yVector;0];
        rotationMatrix = [xres*xVector yres*yVector zres*zVector originShift]; % turn these vectors into matrices
    else 
        sx = info.sx;
        sy = info.sy;
        sz = info.sz;
        originShift = [sx;sy;sz;1];
        
        matrixx = info.matrixx;
        matrixy = info.matrixy;
        
        ix = info.ix;
        iy = info.iy;
        iz = info.iz;
        jx = info.jx;
        jy = info.jy;
        jz = info.jz;
        kx = info.kx;
        ky = info.ky;
        kz = info.kz;
        
        xVector = round([ix;iy;iz;0],8); % what direction rows run w/r/to x
        yVector = round([jx;jy;jz;0],8); % what direction the cols run w/r/to y
        zVector = round([kx;ky;kz;0],8); % what direction the cols run w/r/to y
        rotationMatrix = [xVector yVector zVector originShift]; % turn these vectors into matrices
    end 
    
    truePositions = zeros(matrixx,matrixy,3);
    for i=1:matrixx
        for j=1:matrixy
            arrayPosition = double([i;j;sliceNum;1]);
            thisPosition = rotationMatrix*arrayPosition;
            thisPosition(4) = []; % remove dummy dimension
            truePositions(i,j,:) = thisPosition;
        end 
    end 
end
