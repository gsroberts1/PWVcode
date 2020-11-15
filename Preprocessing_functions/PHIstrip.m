dcmDir = dir();
for i=3:length(dcmDir)
    cd(dcmDir(i).name)
    dcmsList = dir('*.dcm');
    for j=1:length(dcmsList)
        dcm = dicomread(dcmsList(j).name);
        info = dicominfo(dcmsList(j).name);
        
        %Strip Info
        info.SourceApplicationEntityTitle = [];
        info.FileModDate = [];
        info.AccessionNumber = [];
        info.Manufacturer = [];
        info.InstitutionName = [];
        info.ReferringPhysicianName = [];
        info.StationName = [];
        info.OperatorsName = [];
        info.NameOfPhysiciansReadingStudy = [];
        info.ManufacturerModelName = [];
        info.PatientGroupLength = [];
        info.PatientBirthDate = [];
        info.PatientName = [];
        info.PatientSex = [];
        info.PatientAge = [];
        info.PatientWeight = [];
        info.PatientID = [];
        info.AdditionalPatientHistory = [];
        info.DeviceSerialNumber = [];
        info.SoftwareVersions = [];
        info.ProtocolName = [];
        info.StudyInstanceUID = [];
        info.SeriesInstanceUID = [];
        info.PerformedLocation = [];
        info.PerformedStationName = [];
        info.PerformedProcedureStepStartDate = [];
        info.PerformedProcedureStepStartTime = [];
        
        fields = fieldnames(info);
        filteredInfo = rmfield(info,fields(contains(fields,'Private')));
        
        dicomwrite(dcm,dcmsList(j).name,filteredInfo);
    end 
    cd ..
end 