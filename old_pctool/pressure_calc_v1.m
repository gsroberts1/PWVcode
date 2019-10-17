%clear
clc
close all

%%%%%%%%%%PRESSURE GUI %%%%%%%%%%%%%%%%%%%
%%% By Kevin Johnson

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Global Data                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Flags determining processing
global visual_flag;
global wss_flag;
global pressure_flag;
global flow_flag;
global save_flag;
visual_flag=0;
wss_flag=0;
pressure_flag=0;
flow_flag=0;
save_flag=0;

%Visual Porcessing flag (save them incase need them later(
global vis_alpha;
global vis_thresh;
global vis_axis;
global old_ind_step;
old_ind_step = 0;
vis_axis = 1;

%WSS Processing
global wss_axis;
wss_axis = 1;
global visc;

%Pressure Processing
global press_axis;
press_axis = 1;
global press_axis_plot;
press_axis_plot=2;

%Raw Data
global MAG;
global CD;
global MASK;
global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;
global sMASK;
global sCD;
global sMAG;

%Masking parameters
global m_alpha;
global m_beta;
global m_iter;
global m_xstart;
global m_xstop;
global m_ystart;
global m_ystop;
global m_zstart;
global m_zstop;
global m_xlength;
global m_ylength;
global m_zlength;

%%%Data Acquisition Parameters
global tframes;
global tres;
global rcxres;
global rcyres;
global rczres;
global xfov;
global yfov;
global zfov;
global delX;
global delY;
global delZ;

%%Output paramters
global ensight_flag;
ensight_flag =0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Run GUI's (Plans put in for now)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Running Phase Contrast Processing tool');

%%%RUN Options GUI
uiwait(options_gui());
if(visual_flag==1) disp('Performing Visualization');end
if(pressure_flag==1) disp('Performing Pressure Processing');end
if(wss_flag==1) disp('Wall Shear Stress Processing');end
if(flow_flag==1) disp('Performing Flow Measurements');end
if(save_flag==1) disp('Saving Data');end


%%%RUN Segmentation GUI
uiwait(segment_gui());

%%Load Velocity Data%%%%%%%%%%
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
VELX= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

fid=fopen(vy_name,'r');
fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
fclose(fid);
VELY= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

fid=fopen(vz_name,'r');
fseek(fid,2*rcxres*rcyres*(m_zstart-1),'bof');
TEMP= reshape(fread(fid,m_zlength*rcxres*rcyres,'short'),[rcxres rcyres m_zlength]);
fclose(fid);
VELZ= TEMP(m_xstart:m_xstop,m_ystart:m_ystop,:);

if tframes ~= 0
    VELXt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    VELYt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));
    VELZt=single(zeros(m_xlength,m_ylength,m_zlength,tframes));


    for time=0:(tframes-1)

        disp(['Read Frame ',int2str(time),' of ',int2str(tframes)]);

        if tframes ==3
            vx_name ='comp_vd_1.dat';
            vy_name ='comp_vd_2.dat';
            vz_name ='comp_vd_3.dat';
        else
            vx_name = sprintf('ph_%03d_vd_1.dat',time);
            vy_name = sprintf('ph_%03d_vd_2.dat',time);
            vz_name = sprintf('ph_%03d_vd_3.dat',time);
            %             if time <10
            %                 vx_name =['ph_ ',int2str(time),'_vd_1.dat'];
            %                 vy_name =['ph_ ',int2str(time),'_vd_2.dat'];
            %                 vz_name =['ph_ ',int2str(time),'_vd_3.dat'];
            %             else
            %                 vx_name =['ph_',int2str(time),'_vd_1.dat'];
            %                 vy_name =['ph_',int2str(time),'_vd_2.dat'];
            %                 vz_name =['ph_',int2str(time),'_vd_3.dat'];
            %             end
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


%%Clean Up Data %%%%%%%%%%%%%%%%%
if(size(MASK,1)~=0)
    sMASK = MASK(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);
end
sMAG  = MAG(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);
sCD   = CD(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop);

clear global MAG;
clear global CD;
clear global MASK;

%%%RUN Pressure Calulation GUI
if pressure_flag ==1
    uiwait( pressure_gui());
end

%%%RUN Wall Shear Stress GUI
if wss_flag == 1
    uiwait( wss_gui());
end

%%%RUN Visualization GUI
if visual_flag == 1
    uiwait( visual_gui());
end

%%%RUN Flow Maasurement GUI
if flow_flag ==1
    uiwait( flow_gui());
end

%%%RUN Data Out
%uiwait( output_gui());


if ensight_flag ==1

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
            casePathName  = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'.case');
            geoPathName   = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'.geo');
            dataPathName  = sprintf('%s%s%s%s%s',ensightDirPath,'\','EnSight_',patientName,'_');

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

            fwrite(fidGeo,szy,'int');
            fwrite(fidGeo,szx,'int');
            fwrite(fidGeo,numSlices,'int');

            % needed for coordiante definition (KMJ New)
            [Xcor Ycor Zcor] = meshgrid( (0:m_xlength-1)*dx,(0:m_ylength-1)*dy,(0:m_zlength-1)*dz);  %position in mm
            
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
        dataPathFlow = sprintf('%s%s%s',dataPathName,num2str(m-1,'%02d'),'.vel');
        fidMag       = fopen(dataPathMag, 'w', 'ieee-be');
        fidFlow      = fopen(dataPathFlow, 'w', 'ieee-be');
        part         = 1;

        % write data header
        [dummy szStr]= size(sprintf('%s%s','Magnitude, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','Magnitude, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidMag,dataheaderStr,'char');
        [dummy szStr]= size(sprintf('%s%s','Velocity, time ',num2str(m-1)));
        dataheaderStr(1:szStr) = sprintf('%s%s','Velocity, time ',num2str(m-1));
        dataheaderStr(szStr+1:80) = ' ';
        fwrite(fidFlow,dataheaderStr,'char');
        dataheaderStr(1:4) = 'part';
        dataheaderStr(5:80) = ' ';
        fwrite(fidMag,dataheaderStr,'char');
        fwrite(fidFlow,dataheaderStr,'char');

        fwrite(fidMag,part,'int');
        fwrite(fidFlow,part,'int');

        geoheaderStr(1:5) = 'block';
        geoheaderStr(6:80) = ' ';
        fwrite(fidMag,geoheaderStr,'char');
        fwrite(fidFlow,geoheaderStr,'char');



        %%%%%%%%KMJ No longer Temporarily Stored%%%%%
        %             load tempDataStorageMag.mat
        %             % write magnitude data
        %             %dataTmpMag(:,:,:) = squeeze(imageMxMAG(:,:,:,m));
        fwrite(fidMag,sCD(:),'float');
        %            clear dataTmpMag;
        %            clear imageMxMAG;

        %            load tempDataStorageFlow.mat
        % construct velocity data, adjust sign & orientation

        %             imageMxFlow(:,:,:,1) = squeeze(imageMxFlow(:,:,:,1)) * signVx /1000; % convert to mm/msec
        %             imageMxFlow(:,:,:,2) = squeeze(imageMxFlow(:,:,:,2)) * signVy /1000; % convert to mm/msec
        %             imageMxFlow(:,:,:,3) = squeeze(imageMxFlow(:,:,:,3)) * signVz /1000; % convert to mm/msec
        %%%%KMJ new format
        imageMxFlow(:,:,:,2) = squeeze(VELXt(:,:,:,phase)) * signVx /1000; % convert to mm/msec
        imageMxFlow(:,:,:,1) = squeeze(VELYt(:,:,:,phase)) * signVy /1000; % convert to mm/msec
        imageMxFlow(:,:,:,3) = squeeze(VELZt(:,:,:,phase)) * signVz /1000; % convert to mm/msec

        % write velocity data (m/s)
        fwrite(fidFlow,imageMxFlow(:),'float');
        clear imageMxFlow;

        fclose(fidMag);
        fclose(fidFlow);

    end; % end phase loop
end %%ensight conversion




