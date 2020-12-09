anatDatasets(1).Names = 'Sagittal';
anatDatasets(1).Info = dicominfo('S12549.MR.RAD_CMR.0004.0001.2019.03.28.18.38.57.565583.533245258.IMA');
anatDatasets(1).RootDir = 'F:\PWV\Data4Grant\Siemens\Set2\2D_flow\HASTE_SAG_DB_IPAT_0004';
anatDatasets(1).Data = single(Sagittal);

pcDatasets(1).Names = 'AAo';
pcDatasets(1).Info = dicominfo('S12549.MR.RAD_CMR.0040.0001.2019.03.28.18.38.57.565583.533709094.IMA');
pcDatasets(1).RootDir = 'F:\PWV\Data4Grant\Siemens\Set2\2D_flow\2DFLOW_PWV_200_AAO_RETRO_BH_0031';
pcDatasets(1).Data.MAG = single(mean(AAo(:,:,:,1),3));
pcDatasets(1).Data.VMEAN = single(mean(AAo(:,:,:,3),3));
pcDatasets(1).Data.CD = single(mean(AAo(:,:,:,2),3));
pcDatasets(1).Data.mag = single(AAo(:,:,:,1));
pcDatasets(1).Data.cd = single(AAo(:,:,:,2));
pcDatasets(1).Data.v = single(AAo(:,:,:,3));

pcDatasets(2).Names = 'AAo2';
pcDatasets(2).Info = dicominfo('S12549.MR.RAD_CMR.0043.0001.2019.03.28.18.38.57.565583.533718312.IMA');
pcDatasets(2).RootDir = 'F:\PWV\Data4Grant\Siemens\Set2\2D_flow\2DFLOW_PWV_200_DAO1_RETRO_BH_0034';
pcDatasets(2).Data.MAG = single(mean(AAo2(:,:,:,1),3));
pcDatasets(2).Data.VMEAN = single(mean(AAo2(:,:,:,3),3));
pcDatasets(2).Data.CD = single(mean(AAo2(:,:,:,2),3));
pcDatasets(2).Data.mag = single(AAo2(:,:,:,1));
pcDatasets(2).Data.cd = single(AAo2(:,:,:,2));
pcDatasets(2).Data.v = single(AAo2(:,:,:,3));

pcDatasets(3).Names = 'AbdAo';
pcDatasets(3).Info = dicominfo('S12549.MR.RAD_CMR.0046.0001.2019.03.28.18.38.57.565583.533725482.IMA');
pcDatasets(3).RootDir = 'F:\PWV\Data4Grant\Siemens\Set2\2D_flow\2DFLOW_PWV_200_DAO2_RETRO_BH_0037';
pcDatasets(3).Data.MAG = single(mean(AbdAo(:,:,:,1),3));
pcDatasets(3).Data.VMEAN = single(mean(AbdAo(:,:,:,3),3));
pcDatasets(3).Data.CD = single(mean(AbdAo(:,:,:,2),3));
pcDatasets(3).Data.mag = single(AbdAo(:,:,:,1));
pcDatasets(3).Data.cd = single(AbdAo(:,:,:,2));
pcDatasets(3).Data.v = single(AbdAo(:,:,:,3));
