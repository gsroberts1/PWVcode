% Clear the command window and all MATLAB windows
clc
close all

%%%%%%%%%%PC Tool GUI%%%%%%%%%%%%%%%%%%%
%%% Initials by  Kevin Johnson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Global Data                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flags determining processing
global visual_flag; visual_flag = 0;
global wss_flag; wss_flag = 0;
global pressure_flag; pressure_flag = 0;
global flow_flag; flow_flag = 0;
global save_flag; save_flag = 0;
global frame_by_frame; frame_by_frame = 0;
global pwv_flag; pwv_flag = 0; %alw

% Visual Processing flag (save them in case we need them later)
global vis_alpha vis_thresh;
global vis_axis; vis_axis = 1;
global old_ind_step; old_ind_step = 0;

%WSS Processing
global wss_axis; wss_axis = 1;
global visc;

%Pressure Processing
global press_axis; press_axis = 1;
global press_axis_plot; press_axis_plot = 2;

%Raw Data
global MAG CD MASK;
global VELX VELY VELZ;
global VELXt VELYt VELZt;
global sMASK sCD sMAG PRESSURE;
global plist;

% Masking parameters
global m_alpha m_beta;
global m_iter;
global m_xstart m_xstop;
global m_ystart m_ystop;
global m_zstart m_zstop;
global m_xlength m_ylength m_zlength;

% Data Acquisition Parameters
global tframes tres;
global rcxres rcyres rczres;
global xfov yfov zfov;
global delX delY delZ;

% Output paramters
global ensight_flag; ensight_flag =0;
global cine_anatomy_flag; cine_anatomy_flag=0;
global sCD_cine sMAG_cine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Run GUI's (Plans put in for now)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run Segmentation GUI
disp('Running Phase Contrast Processing Tool...');
uiwait(segment_gui());
disp('done.');

% Determines the processing after segmentation
if visual_flag == 1
    disp('Performing Visualization');
elseif pressure_flag == 1
    disp('Performing Pressure Processing');
elseif wss_flag == 1
    disp('Wall Shear Stress Processing');
elseif flow_flag == 1
    disp('Running the flow measurement tool...');
elseif save_flag == 1
    disp('Saving Data...');
elseif pwv_flag == 1
    disp('Running PWV tool...');
else
    disp('No processing options were selected.');
end

% Load Velocity Data
VELX=single(zeros(m_xlength,m_ylength,m_zlength));
VELY=single(zeros(m_xlength,m_ylength,m_zlength));
VELZ=single(zeros(m_xlength,m_ylength,m_zlength));

vx_name ='comp_vd_1.dat';
vy_name ='comp_vd_2.dat';
vz_name ='comp_vd_3.dat';

fid=fopen(vx_name,'r');
fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
fclose(fid);
VELX = TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

fid=fopen(vy_name,'r');
fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
fclose(fid);
VELY = TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

fid=fopen(vz_name,'r');
fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
fclose(fid);
VELZ = TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

if tframes ~= 0 && frame_by_frame == 0
    VELXt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    VELYt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    VELZt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    for time=0:(tframes-1)

        disp(['Read Frame ',int2str(time+1),' of ',int2str(tframes)]);%alw added 1

        if tframes == 3
            vx_name = 'comp_vd_1.dat';
            vy_name = 'comp_vd_2.dat';
            vz_name = 'comp_vd_3.dat';
        else
            vx_name = sprintf('ph_%03d_vd_1.dat',time);
            vy_name = sprintf('ph_%03d_vd_2.dat',time);
            vz_name = sprintf('ph_%03d_vd_3.dat',time);
        end

        fid=fopen(vx_name,'r');
        fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
        fclose(fid);
        VELXt(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

        fid=fopen(vy_name,'r');
        fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
        fclose(fid);
        VELYt(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

        fid=fopen(vz_name,'r');
        fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
        fclose(fid);
        VELZt(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);
    end

end

global turb_flag; turb_flag = 0;

if turb_flag==1
    sTURB_cine=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    for time=0:(tframes-1)
        disp(['Read Frame ',int2str(time),' of ',int2str(tframes)]);
        mag_name = sprintf('ph_%03d_turb.dat',time);
        
        fid=fopen(mag_name,'r');
        fseek(fid,4*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'float'),[rcxres rcyres m_zlength]);
        fclose(fid);
        sTURB_cine(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);
    end
end

%%%%%%%%%%%%%export time resolved magnitude and complex difference
if cine_anatomy_flag==1 && frame_by_frame==0
    sCD_cine=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    sMAG_cine=single(zeros(m_xlength,m_ylength,m_zlength,tframes));

    for time=0:(tframes-1)
        disp(['Read Frame ',int2str(time),' of ',int2str(tframes)]);
        mag_name = sprintf('ph_%03d_mag.dat',time);
        cd_name = sprintf('ph_%03d_cd.dat',time);

        fid=fopen(mag_name,'r');
        fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
        fclose(fid);
        sMAG_cine(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

        fid=fopen(cd_name,'r');
        fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
        TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
        fclose(fid);
        sCD_cine(:,:,:,time+1)= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);
    end
end
clear TEMP;

%%Clean Up Data %%%%%%%%%%%%%%%%%
if(size(MASK,1)~=0)
    sMASK = MASK(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);
end
sMAG = MAG(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);
sCD  = CD(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);

clear global MAG; clear global CD; clear global MASK;

% Run the appropriate GUI based on the user's choice in the
% Processing Options panel
if pressure_flag == 1
    uiwait(pressure_gui());
elseif wss_flag == 1
    uiwait(wss_gui());
elseif flow_flag == 1
    uiwait( flow_gui());
elseif pressure_flag == 1
    grad_export_flag = 1;
elseif pwv_flag == 1
    pwv_tool;
else
    grad_export_flag = 0;
    disp('No processing options were selected.');
end

global GRADx GRADy GRADz;

% Linux Mode
linux = 0;
if linux == 1
    slash = '/'
else
    slash = '\'
end

if ensight_flag == 1

    %% --------------------------------------------------------
    %% parameters, need to be manually entered for each data set
    %% ---------------------------------------------------------
    patientName     = 'patient';

    % EnSight conversion flags, needed for velocity sign & orientation
    % with respect to anatomical (magnitude) images
    % data is accessed in image (x,y,z) coordinates as loaded by matlab
    % x: L/R for coronal & axial, A/P for sagittal
    %    read & PE directions
    % y: S/I for coronal & sagittal, A/P for  axial
    %    read and PE directions
    % z: through slab, A/P for coronal, R/L for sagittal, S/I for axial
    %    slice/partition direction (numSlices)
    %  orientation    = 'tra';        % main slab / slice orientation, other options: 'tra','cor'
    %  peDir          = 'y';          % phase encoding direction, other option: 'y'
    % get image size from MR data

    numSlices      = size(sCD,3);   % number of slices / partitions
    szy            = size(sCD,1);
    szx            = size(sCD,2);


    dx = delX;
    dy = delY;
    dz = delZ;

    for phase = 1:tframes

        %% -------------------------------------------------
        %% Init EnSight file conversion
        %% -------------------------------------------------

        if phase == 1
            % update text output
            outStr = ['Dicom --> EnSight, generating case & geo files ... '];
            disp(outStr);
            drawnow;

            % create new directory for EnSight & matlab files
            ensightDirPath  = sprintf('%s%s','EnSight_',patientName);
            [s,mess,messid] = mkdir(ensightDirPath);

            % generate case file
            % ------------------------------------------------
            casePathName  = sprintf('%s%s%s%s%s',ensightDirPath,slash,'EnSight_',patientName,'.case');
            geoPathName   = sprintf('%s%s%s%s%s',ensightDirPath,slash,'EnSight_',patientName,'.geo');
            dataPathName  = sprintf('%s%s%s%s%s',ensightDirPath,slash,'EnSight_',patientName,'_');

            geoFileName   = sprintf('EnSight_%s%s%',patientName,'.geo');
            dataFileName  = sprintf('EnSight_%s%s%',patientName,'_');


            timeStamps   =   (1:tframes)*tres; %time in ms

            % open text file an write data
            fidCase  = fopen(casePathName, 'wt');
            fprintf(fidCase,'FORMAT\n');
            fprintf(fidCase,'type:	 ensight gold\n');
            fprintf(fidCase,'GEOMETRY\n');
            fprintf(fidCase,'model:	 %s\n',geoFileName);
            fprintf(fidCase,'VARIABLE\n');
            fprintf(fidCase,'scalar per node:	 Magnitude	 %s**.mag\n',dataFileName );
            fprintf(fidCase,'scalar per node:	 ComplexDifference	 %s**.cd\n',dataFileName );
            
            if(size(sMASK,1)~=0)
                fprintf(fidCase,'scalar per node:	 Mask	 %s**.mask\n',dataFileName );
            end
            if cine_anatomy_flag==1
                fprintf(fidCase,'scalar per node:	 MagnitudeCINE	 %s**.mag_cine\n',dataFileName );
                fprintf(fidCase,'scalar per node:	 AngiogramCINE	 %s**.cd_cine\n',dataFileName );
            end
            if turb_flag==1
                fprintf(fidCase,'scalar per node:	 MagnitudeTURB	 %s**.turb_cine\n',dataFileName );
            end
            
            
            if pressure_flag ==1
                fprintf(fidCase,'scalar per node:	 Pressure	 %s**.pressure\n',dataFileName );
            end
            if grad_export_flag ==1
                fprintf(fidCase,'vector per node:	 P_Gradient	 %s**.pgrad\n',dataFileName );
            end
            fprintf(fidCase,'vector per node:	 Velocity	 %s**.vel\n',dataFileName );
            fprintf(fidCase,'TIME\n');
            fprintf(fidCase,'time set:		 1\n');
            fprintf(fidCase,'number of steps:	 %s\n',num2str(tframes));
            fprintf(fidCase,'filename start number:	 0\n');
            fprintf(fidCase,'filename increment:	 1\n');
            fprintf(fidCase,'time values:\n');
            fprintf(fidCase,'%.3f\n',timeStamps);
            fclose(fidCase); % close file
            % end generate case file


            % generate geo file
            % ------------------------------------------------
            % open binary file an write data
            fidGeo  = fopen(geoPathName, 'w', 'ieee-be');
            matSize = szx * szy;
            part    = 1;

            % generate and write header
            geoheaderStr(1:8) = 'C Binary';
            geoheaderStr(9:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:21) = 'Ensight Geometry File';
            geoheaderStr(22:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:44) = 'Created by MATLAB routine, (c) M. Markl 2006';
            geoheaderStr(45:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:14) = 'node id assign';
            geoheaderStr(15:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:17) = 'element id assign';
            geoheaderStr(18:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:4) = 'part';
            geoheaderStr(5:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');

            fwrite(fidGeo,part,'int');

            geoheaderStr(1:12) = 'Total Volume';
            geoheaderStr(13:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');
            geoheaderStr(1:5) = 'block';
            geoheaderStr(6:80) = ' ';
            fwrite(fidGeo,geoheaderStr,'char');

            %%KMJ Swap x/y
            fwrite(fidGeo,szy,'int');
            fwrite(fidGeo,szx,'int');
            fwrite(fidGeo,numSlices,'int');

            % needed for coordiante definition (KMJ New)
            [Xcor Ycor Zcor] = meshgrid( (0:m_ylength-1)*dx,(0:m_xlength-1)*dy,(0:m_zlength-1)*dz);  %position in mm

            %             for l=1:szy
            %                 Xtmp(l,1:szx) = [1:szx];
            %             end;
            %
            %             for l=1:szx
            %                 Ytmp(1:szy,l) = fliplr(1:szy);
            %             end;
            %
            %             for k = 1:numSlices
            %                 Xcor((k-1)*matSize+1:k*matSize) = Xtmp(:) * dx;
            %                 Ycor((k-1)*matSize+1:k*matSize) = Ytmp(:) * dy;
            %                 Zcor((k-1)*matSize+1:k*matSize) = k * dz;
            %             end

            fwrite(fidGeo,Xcor,'float'); % x-coordinate
            fwrite(fidGeo,Ycor,'float'); % y-coordinate
            fwrite(fidGeo,Zcor,'float'); % z-coordinate
            fclose(fidGeo); % close geo file
            % end generate geo file

            clear Xcor;
            clear Ycor;
            clear Zcor;

        end %% end phase = 1

        %% ------------------------------------------------------
        %% EnSight file conversion, generate EnSight data files
        %% ------------------------------------------------------

        %%%KMJ Values are based on logical flow axis
        signVx = -1;
        signVy = -1;
        signVz = -1;


        %         % arrange / mirror velocity data according to the main
        %         % orientation and in-plane phase encoding direction
        %         if orientation == 'cor'
        %             signVx = 1;
        %             signVy = -1;
        %             signVz = -1;
        %         elseif orientation == 'sag'
        %             signVx = 1;
        %             signVy = 1;
        %             signVz = 1;
        %         elseif orientation == 'tra'
        %             signVx = -1;
        %             signVy = 1;
        %             signVz = -1;
        %         else
        %             outStr = ['Dicom --> EnSight Error: Main orientation not recognized - aborting !! '];
        %             display('outStr');
        %             drawnow;
        %             return;
        %         end

        %%for m = 1:numPhases
        m = phase;

        % update text output
        outStr = ['Dicom --> EnSight, generating EnSight data files , phase ',num2str(m),' ...'];
        disp(outStr);
        drawnow;

        % generate data files
        % ------------------------------------------------
        % open binary file an write data
        dataPathMag  = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.mag');
        dataPathCD  = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.cd');
        dataPathFlow = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.vel');
        dataPathMask  = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.mask');
        dataPathPressure= sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.pressure');
        dataPathPGrad= sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.pgrad');
        dataPathMagCine= sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.mag_cine');
        dataPathCDCine= sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.cd_cine');
        dataPathTURBCine= sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.turb_cine');

        fidMag       = fopen(dataPathMag, 'w', 'ieee-be');
        fidCD        = fopen(dataPathCD, 'w', 'ieee-be');
        fidFlow      = fopen(dataPathFlow, 'w', 'ieee-be');
        if(size(sMASK,1)~=0)
            fidMask       = fopen(dataPathMask, 'w', 'ieee-be');
        end
        if cine_anatomy_flag==1
            fidMagCine       = fopen(dataPathMagCine, 'w', 'ieee-be');
            fidCDCine        = fopen(dataPathCDCine, 'w', 'ieee-be');
        end
        if turb_flag==1
            fidTURBCine       = fopen(dataPathTURBCine, 'w', 'ieee-be');
        end
        if pressure_flag==1
            fidPressure  = fopen(dataPathPressure, 'w', 'ieee-be');
        end
        if grad_export_flag==1
            fidPGrad  = fopen(dataPathPGrad, 'w', 'ieee-be');
        end
        part         = 1;

        % write data header
        [dummy szStr]= size(sprintf('%s%s','Magnitude, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','Magnitude, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidMag,dataheaderStr,'char');

        [dummy szStr]= size(sprintf('%s%s','ComplexDifference, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','ComplexDifference, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidCD,dataheaderStr,'char');

        if(size(sMASK,1)~=0)
        [dummy szStr]= size(sprintf('%s%s','Mask, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','Mask, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidMask,dataheaderStr,'char');
        end
        
        [dummy szStr]= size(sprintf('%s%s','Velocity, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','Velocity, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidFlow,dataheaderStr,'char');

        if cine_anatomy_flag ==1
            [dummy szStr]= size(sprintf('%s%s','MagnitudeCINE, time ',num2str(m-1)));
            dataheaderStr(1:szStr) = sprintf('%s%s','MagnitudeCINE, time ',num2str(m-1));
            dataheaderStr(szStr+1:80) = ' ';
            fwrite(fidMagCine,dataheaderStr,'char');

            [dummy szStr]= size(sprintf('%s%s','MagnitudeCINE, time ',num2str(m-1)));
            dataheaderStr(1:szStr) = sprintf('%s%s','MagnitudeCINE, time ',num2str(m-1));
            dataheaderStr(szStr+1:80) = ' ';
            fwrite(fidCDCine,dataheaderStr,'char');
        end
        
        if turb_flag==1
            [dummy szStr]= size(sprintf('%s%s','TurbCINE, time ',num2str(m-1)));
            dataheaderStr(1:szStr) = sprintf('%s%s','TurbCINE, time ',num2str(m-1));
            dataheaderStr(szStr+1:80) = ' ';
            fwrite(fidTURBCine,dataheaderStr,'char');
        end
        

        if pressure_flag == 1
            [dummy szStr]= size(sprintf('%s%s','Pressure, time ',num2str(m-1)));
            dataheaderStr(1:szStr) = sprintf('%s%s','Pressure, time ',num2str(m-1));
            dataheaderStr(szStr+1:80) = ' ';
            fwrite(fidPressure,dataheaderStr,'char');
        end

        if grad_export_flag == 1
            [dummy szStr]= size(sprintf('%s%s','P_Gradient, time ',num2str(m-1)));
            dataheaderStr(1:szStr) = sprintf('%s%s','P_Gradient, time ',num2str(m-1));
            dataheaderStr(szStr+1:80) = ' ';
            fwrite(fidPGrad,dataheaderStr,'char');
        end

        dataheaderStr(1:4) = 'part';
        dataheaderStr(5:80) = ' ';
        fwrite(fidMag,dataheaderStr,'char');
        fwrite(fidCD,dataheaderStr,'char');
        fwrite(fidFlow,dataheaderStr,'char');
        if(size(sMASK,1)~=0)
            fwrite(fidMask,dataheaderStr,'char');
        end
        if cine_anatomy_flag == 1
            fwrite(fidMagCine,dataheaderStr,'char');
            fwrite(fidCDCine,dataheaderStr,'char');
        end
        if turb_flag == 1
            fwrite(fidTURBCine,dataheaderStr,'char');
        end
        
        if pressure_flag == 1
            fwrite(fidPressure,dataheaderStr,'char');
        end
        if grad_export_flag == 1
            fwrite(fidPGrad,dataheaderStr,'char');
        end

        fwrite(fidMag,part,'int');
        fwrite(fidFlow,part,'int');
        fwrite(fidCD,part,'int');
        if(size(sMASK,1)~=0)
         fwrite(fidMask,part,'int');
        end
        
        if cine_anatomy_flag == 1
            fwrite(fidMagCine,part,'int');
            fwrite(fidCDCine,part,'int');
        end
        if turb_flag == 1
            fwrite(fidTURBCine,part,'int');
        end
        
        if pressure_flag ==1
            fwrite(fidPressure,part,'int');
        end
        if grad_export_flag ==1
            fwrite(fidPGrad,part,'int');
        end

        geoheaderStr(1:5) = 'block';
        geoheaderStr(6:80) = ' ';
        fwrite(fidMag,geoheaderStr,'char');
        fwrite(fidCD,geoheaderStr,'char');
        fwrite(fidFlow,geoheaderStr,'char');
        if(size(sMASK,1)~=0)
            fwrite(fidMask,geoheaderStr,'char');
        end
        if cine_anatomy_flag == 1
            fwrite(fidMagCine,geoheaderStr,'char');
            fwrite(fidCDCine,geoheaderStr,'char');
        end
        if turb_flag == 1
            fwrite(fidTURBCine,geoheaderStr,'char');
        end
        if pressure_flag == 1
            fwrite(fidPressure,geoheaderStr,'char');
        end
        if grad_export_flag == 1
            fwrite(fidPGrad,geoheaderStr,'char');
        end

        %%%%%%%%KMJ No longer Temporarily Stored%%%%%
        if(frame_by_frame==0)

            fwrite(fidMag,sMAG(:),'float');
            fwrite(fidCD,sCD(:),'float');
            if(size(sMASK,1)~=0)
            fwrite(fidMask,sMASK(:),'float'); 
            end
            %%%%KMJ new format
            imageMxFlow(:,:,:,2) = squeeze(VELXt(:,:,:,phase)) * signVx /1000; % convert to mm/msec
            imageMxFlow(:,:,:,1) = squeeze(VELYt(:,:,:,phase)) * signVy /1000; % convert to mm/msec
            imageMxFlow(:,:,:,3) = squeeze(VELZt(:,:,:,phase)) * signVz /1000; % convert to mm/msec

            % write velocity data (m/s)
            fwrite(fidFlow,imageMxFlow,'float');
            clear imageMxFlow;

            if cine_anatomy_flag==1
                fwrite(fidMagCine,sMAG_cine(:,:,:,phase),'float');
                fwrite(fidCDCine,sCD_cine(:,:,:,phase),'float');
            end
            if turb_flag==1
                fwrite(fidTURBCine,sTURB_cine(:,:,:,phase),'float');
            end
            

            %%%Write Pressure Gradient
            if grad_export_flag ==1
                imagePress = zeros(size(sCD));
                imagePress(plist) = GRADx(:,phase);
                imageMxFlow(:,:,:,1) = imagePress* signVx;

                imagePress(plist) = GRADy(:,phase);
                imageMxFlow(:,:,:,2) = imagePress* signVy;

                imagePress(plist) = GRADz(:,phase);
                imageMxFlow(:,:,:,3) = imagePress* signVz;

                fwrite(fidPGrad,imageMxFlow(:),'float');
                clear imageMxFlow;
                clear imagePress;
            end

            if pressure_flag==1

                imagePress = zeros(size(sCD));
                imagePress(plist) = 0.007501*PRESSURE(:,phase);

                % write pressure (mmHg)
                fwrite(fidPressure,imagePress(:),'float');
                fclose(fidPressure);

                clear imagePress;
            end
        else
            
            %%allocate these matrices
            if phase==1 
                VELXt=single(zeros(m_xlength,m_ylength,m_zlength));
                VELYt=single(zeros(m_xlength,m_ylength,m_zlength));
                VELZt=single(zeros(m_xlength,m_ylength,m_zlength));
                sCD_cine=single(zeros(m_xlength,m_ylength,m_zlength));
                sMAG_cine=single(zeros(m_xlength,m_ylength,m_zlength));
                TEMP = single(zeros(rcxres,rcyres,m_zlength));
                TURB = single(zeros(rcxres,rcyres,m_zlength));
            end
                
            %%%Read Frames on Demand
            fwrite(fidMag,sMAG(:),'float');
            fwrite(fidCD,sCD(:),'float');
            if(size(sMASK,1)~=0)
                fwrite(fidMask,sMASK(:),'float'); 
            end
            outStr = ['Dicom --> EnSight, generating EnSight data files , phase Frame by Frame Read...'];
            disp(outStr);

            vx_name = sprintf('ph_%03d_vd_1.dat',phase-1);
            vy_name = sprintf('ph_%03d_vd_2.dat',phase-1);
            vz_name = sprintf('ph_%03d_vd_3.dat',phase-1);

            fid=fopen(vx_name,'r');
            fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
            TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
            fclose(fid);
            VELXt= signVx/1000*TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

            fid=fopen(vy_name,'r');
            fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
            TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
            fclose(fid);
            VELYt= signVy/1000*TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

            fid=fopen(vz_name,'r');
            fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
            TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
            fclose(fid);
            VELZt= signVz/1000*TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

            %%%%Write out data here
            fwrite(fidFlow,VELYt,'float');
            fwrite(fidFlow,VELXt,'float');
            fwrite(fidFlow,VELZt,'float');

            if cine_anatomy_flag==1

                %%%%%%%%%%%%%export time resolved magnitude and complex difference
                mag_name = sprintf('ph_%03d_mag.dat',phase-1);
                cd_name = sprintf('ph_%03d_cd.dat',phase-1);

                fid=fopen(mag_name,'r');
                fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
                TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
                fclose(fid);
                sMAG_cine= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

                fid=fopen(cd_name,'r');
                fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
                TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
                fclose(fid);
                sCD_cine= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

                fwrite(fidMagCine,sMAG_cine,'float');
                fwrite(fidCDCine,sCD_cine,'float');
            end
            
            if turb_flag ==1
                 %%%%%%%%%%%%%export time resolved magnitude and complex difference
                mag_name = sprintf('ph_%03d_turb.dat',phase-1);
               
                fid=fopen(mag_name,'r');
                fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
                TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
                fclose(fid);
                sTURB_cine= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

                fwrite(fidTURBCine,sTURB_cine,'float');
            end
        end
        fclose(fidCD);
        fclose(fidMag);
        fclose(fidFlow);
        if(size(sMASK,1)~=0)
            fclose(fidMask);
        end
    end; % end phase loop
end %%ensight conversion