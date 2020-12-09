dicomDir = dir();
for j=3:length(dicomDir)
    cd(dicomDir(j).name); %go into dicom folder
    tar([dicomDir(j).name '.tgz'],{'*.dcm'});
    delete('*.dcm');
    cd('..') %move back up to directory with all dicoms
end 