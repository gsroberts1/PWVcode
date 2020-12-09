imList = dir('IM*');
anatList = imList(end-7:end);
pcList = imList(1:end-8);
xxList = dir('XX*');
psList = dir('PS*');
numThrowOut = length(xxList) + length(psList);

magIter = 1;
cdIter = 1;
vIter = 1;
for i=1:length(pcList)
    dicomInfo = dicominfo(pcList(i).name);
    ImageType = dicomInfo.ImageType;
    
    if strcmp(ImageType,'ORIGINAL\PRIMARY\M_FFE\M\FFE')
        mag(:,:,magIter) = dicomread(pcList(i).name);
        magIter = magIter + 1;
    end 
    
    if strcmp(ImageType,'ORIGINAL\PRIMARY\M_PCA\M\PCA')
        cd(:,:,cdIter) = dicomread(pcList(i).name);
        cdIter = cdIter + 1;
    end 
    
    if strcmp(ImageType,'ORIGINAL\PRIMARY\PHASE CONTRAST M\P\PCA')
        v(:,:,vIter) = dicomread(pcList(i).name);
        vIter = vIter + 1;
    end  
end 

for j=1:length(anatList)
    sagittal(:,:,j) = dicomread(anatList(j).name);
end 


