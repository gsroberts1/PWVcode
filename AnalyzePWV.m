function varargout = AnalyzePWV(varargin)
    % ANALYZEPWV MATLAB code for AnalyzePWV.fig
    %      ANALYZEPWV, by itself, creates a new ANALYZEPWV or raises the existing
    %      singleton*.
    %
    %      H = ANALYZEPWV returns the handle to a new ANALYZEPWV or the handle to
    %      the existing singleton*.
    %
    %      ANALYZEPWV('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in ANALYZEPWV.M with the given input arguments.
    %
    %      ANALYZEPWV('Property','Value',...) creates a new ANALYZEPWV or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before AnalyzePWV_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to AnalyzePWV_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help AnalyzePWV

    % Last Modified by GUIDE v2.5 23-Jul-2020 22:54:21
    
    
    
    % Developed by Grant S Roberts, University of Wisconsin-Madison, 2019
    %   Used by: LoadPWV.m
    %   Dependencies: NONE
    

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @AnalyzePWV_OpeningFcn, ...
                       'gui_OutputFcn',  @AnalyzePWV_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT


function AnalyzePWV_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to AnalyzePWV (see VARARGIN)

    % Choose default command line output for AnalyzePWV
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % Get anatomical, pc, and bSSFP data from LoadPWV GUI
    handles.anatDatasets = varargin{1}; 
    handles.pcDatasets = varargin{2}; 
    handles.magDatasets = varargin{3}; 
    handles.pcDatasets(1).ROI = []; %initialize ROI field
    handles.magDatasets(1).ROI = []; %initialize ROI field
    handles.anatDatasets(1).Centerline = []; %initialize centerline field
    handles.pcDatasets(1).ROIdata = []; %throgh-plane flow curve data
    handles.pcDatasets(1).ROIdataGaussian = []; %flow smoothed with Gauss.
    handles.pcDatasets(1).ROIdataSpline = []; %flow smoothed with spline
    
    % For analysis on images if they are zoomed in 
    handles.global.zoomed = 0; %flag for if PC image is zoomed in GUI
    handles.global.spanUD = [0,0]; %if zoomed, what is the new span of rows
    handles.global.spanLR = [0,0]; %if zoomed, what is the new span of cols
    handles.global.zoomedAnat = 0; %flag for if anatomical image is zoomed
    handles.global.spanUDAnat = [0,0];
    handles.global.spanLRAnat = [0,0];
    
    handles.global.interpType = 'Gaussian'; %interpolation type (i.e. Gaussian)
    handles.global.pgShift = 0; %flag for shifting flow curve
    handles.global.showErrorBars = 0; %flag for showing error bars
    handles.global.startAnalyzing = 0; %flag to begin PWV calculations
    handles.global.homeDir = pwd;
    handles.global.totalROIs = 0;
    
    if numel(handles.magDatasets)==0 %if we didn't load bSSFP data, turn off handle
        set(handles.magRadio,'Enable','off')
    end 
    set(handles.pcPlanePopup,'String',{handles.pcDatasets.Names}); %list of all planes (AAo, AbdAo, etc.)
    set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data)); %list of all datasets (CD, MAG, v, etc.)
    set(handles.interpolatePopup,'String',{'Gaussian','Spline','None'}); %set all possible interpolation types
    set(handles.anatListbox,'String',{handles.anatDatasets.Names}); %list of all anatomical datasets (Axial, Coronal, etc.)
    set(handles.interpolatePopup,'Enable','off'); %initialize radios and buttons
    set(handles.showPCPlanesRadio,'Enable','off'); 
    set(handles.loadCenterline,'Enable','off'); 
    set(handles.drawCenterlineButton,'Enable','off');
    set(handles.deleteCenterlineButton,'Enable','off'); 
    set(handles.computePWVButton,'Enable','off'); 
    set(handles.ttpointRadio,'Value',1); 
    set(handles.ttuRadio,'Value',1);
    set(handles.ttfRadio,'Value',1);
    set(handles.xcorrRadio,'Value',1);
        
    guidata(hObject, handles);
    updatePCImages(handles); %update PC images on top plot
    updateAnatImages(handles); %update Anatomical images on bottom plot

% --- Outputs from this function are returned to the command line.
function varargout = AnalyzePWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    
%%%%%%%%%%%% LOAD 2DPC PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PLANE PLOT - CREATE FUNCTION
function pcPlanePlot_CreateFcn(hObject, eventdata, handles)


% --- PLANE SLIDER - CALLBACK
function pcSlider_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images if slider is moved

% --- PLANE SLIDER - CREATE FUNCTION
function pcSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

    
% --- PC RADIO - CALLBACK
function pcRadio_Callback(hObject, eventdata, handles)
    % If we push this radio, update popups and images accordingly.
    set(handles.pcPlanePopup,'String',{handles.pcDatasets.Names}); %list names of all pc slices (AAo,AbdAo)
    set(handles.pcDatasetPopup,'Enable','on'); %unlock listbox for selecting datasets
    set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data)); %list all datasets for 1st slice (mag,cd,v,MAG,etc.)
    set(handles.magRadio,'Value',0); %turn off mag radio (only one should be selected at once)
    updatePCImages(handles); %update images on PC plot 

    
% --- MAG RADIO - CALLBACK
function magRadio_Callback(hObject, eventdata, handles)
    set(handles.pcPlanePopup,'String',{handles.magDatasets.Names}); %list names of all bSSFP slices (AAo,AbdAo)
    set(handles.pcDatasetPopup,'Enable','off'); %no need for dataset popup (we only have mangitude bSSFP data\
    set(handles.pcDatasetPopup,'String',' '); %remove any text
    set(handles.pcRadio,'Value',0); %turn off pc radio (only one should be selected at once)
    updatePCImages(handles); %update images on PC plot 
    

% --- PLANE DROPDOWN - CALLBACK
function pcPlanePopup_Callback(hObject, eventdata, handles)   
    handles.global.zoomed = 0; %set image to full size
    handles.global.spanUD = 0; 
    handles.global.spanLR = 0;
    
    planeNum = get(handles.pcPlanePopup,'Value'); 
    if get(handles.pcRadio,'Value') %if we are looking at PC images
        set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(planeNum).Data)); %list all datasets for slice of interest
    end 
    set(handles.pcDatasetPopup,'Value',1); %set dataset to 1st in list (usually time-resolved complex difference)
    
    guidata(hObject, handles);
    updatePCImages(handles); %update images on PC plot anytime we click on a new plane

% --- PLANE DROPDOWN - CREATE FUNCTION
function pcPlanePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- DATASET DROPDOWN - CALLBACK
function pcDatasetPopup_Callback(hObject, eventdata, handles)
    updatePCImages(handles); %update images on PC plot anytime we click on a new dataset

% --- DATASET DROPDOWN - CREATE FUNCTION
function pcDatasetPopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ZOOM BUTTON - CALLBACK
function pcZoomButton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on PC plot (not anatomical plot on bottom)
    set(handles.messageBar,'String','Draw an ROI to zoom in.');
    rect = drawrectangle; %draw a rectangular field of view
    positions = round(rect.Position); %get positions of rectangle edges
    spanUD = handles.global.spanUD; %get current span of rows
    spanLR = handles.global.spanLR; %get current span of columns
    
    % Update row span with rectangle vertices (note that we can only zoom in (hence the plus)
    handles.global.spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4)); 
    handles.global.spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));
    handles.global.zoomed = 1; %set flag to zoomed
    
    guidata(hObject, handles);
    updatePCImages(handles); %update to show zoomed in images
 
    
% --- UNZOOM BUTTON - CALLBACK
function pcUnzoomButton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on top PC plot
    handles.global.zoomed = 0; %set flag to unzoomed
    handles.global.spanUD = 0; %reset spans to 0
    handles.global.spanLR = 0;
    
    guidata(hObject,handles);
    updatePCImages(handles); %update to show full size image

    
% --- DRAWROI BUTTON - CALLBACK
function drawROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on top PC plot 
    if get(handles.pcRadio,'Value') %if we are interested in PC (not bSSFP)
        set(handles.magRadio,'Enable','off'); %turn off bSSFP (mag) radio to not create errors
    else 
        set(handles.pcRadio,'Enable','off'); %or turn off PC radio
    end 
    set(handles.drawROIbutton,'Enable','off'); %make it so we can't draw a second ROI
    set(handles.pcPlanePopup,'Enable','off'); %make it so we can't select another plane
    set(handles.pcZoomButton,'Enable','off'); %make it so we can't zoom
    set(handles.pcUnzoomButton,'Enable','off'); %make it so we can't unzoom
    
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane of interest (eg AAo)
    mydlg = warndlg('Press enter when the ROI is set'); %open dialog warning 
    waitfor(mydlg); %MUST PRESS ENTER TO PROCEED
    circle = drawcircle('FaceAlpha',0.1,'Color','g','LineWidth',1,'Deletable',0); %draw circle on PC image
    while true
        w = waitforbuttonpress; %wait for enter push ...
        switch w 
            case 1 % if it was a keyboard press.
              key = get(gcf,'currentcharacter'); %get key that was pressed
                  switch key
                      case 27 % 27 is the escape key
                          set(handles.messageBar,'String','User pressed the escape key. Deleting ROI.');
                          deleteROIbutton_Callback(hObject, eventdata, handles); %kill the current ROI
                          break % break out of the while loop
                      case 13 % 13 is the enter/return key 
                          set(handles.messageBar,'String','ROI selected.');
                          circle.InteractionsAllowed = 'none'; %freeze circle
                          break
                      otherwise 
                          %wait for a different command
                  end
       end
    end

    if get(handles.pcRadio,'Value') %if we're doing PC analysis..
        handles.pcDatasets(planeNum).ROI = circle; %temporarily hold circle (will get deleted once ROI is gone)
    else %if we're doing bSSFP area analysis..
        handles.magDatasets(planeNum).ROI = circle;
    end 

    guidata(hObject,handles);
    updatePCImages(handles); 

    
% --- DELETE ROI BUTTON - CALLBACK
function deleteROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.pcPlanePlot); %make sure we do this on top PC plot 
    hold off %hold on is in the updateImage functions 
    
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane (eg AAo)
    if get(handles.pcRadio,'Value') %if we are doing PC analysis
        handles.pcDatasets(planeNum).ROI = []; %remove ROI object
    else %if we are doing bSSFP area change analysis
        handles.magDatasets(planeNum).ROI = [];
    end 

    set(handles.drawROIbutton,'Enable','on'); %allow for drawing ROIs again
    set(handles.pcPlanePopup,'Enable','on'); %show plane listbox again
    if get(handles.pcRadio,'Value') %turn all PWV radio buttons back on
        set(handles.magRadio,'Enable','on'); 
    else 
        set(handles.pcRadio,'Enable','on'); 
    end 
    set(handles.pcZoomButton,'Enable','on'); %allow for zooming again
    set(handles.pcUnzoomButton,'Enable','on'); %allow for unzooming again

    guidata(hObject,handles);
    updatePCImages(handles); %update PC plot (removes green ROI circle)
    
% --- LOAD ROI BUTTON - CALLBACK
function loadROIbutton_Callback(hObject, eventdata, handles)
    planeNum = get(handles.pcPlanePopup,'Value'); %get current plane (eg AAo)
    handles.global.totalROIs = handles.global.totalROIs + 1; %add 1 to total ROI count
    if get(handles.pcRadio,'Value') %if we're doing pc analysis
        v = handles.pcDatasets(planeNum).Data.v; %grab time-resolved velocity
        if handles.global.zoomed %account for zoom factor
            v = v(handles.global.spanUD,handles.global.spanLR,:);
        end 
        circle = handles.pcDatasets(planeNum).ROI; %pull circle data from handles
        radius = circle.Radius; %get radius of circle
        center = round(circle.Center); %get center coordinates

        [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
        X = X-center(2); %shift coordinate grid
        Y = Y-center(1);
        roiMask = sqrt(X.^2+Y.^2)<=radius; %anything outside radius is ignored

        %%% Create Linear Interpolated Data
        tq = 1:0.1:size(v,3); %interpolate time dimension
        
        if isfield(handles.pcDatasets(planeNum).Info,'matrixx') %if we're dealing with radial data (pcvipr recon)
            matrixx = handles.pcDatasets(planeNum).Info.matrixx; %matrix size in x dimension
            fovx = handles.pcDatasets(planeNum).Info.fovx;  %field of view (mm)
            xres = fovx/matrixx; %resolution (mm). ASSUMED TO BE SAME IN Y DIMENSION
            frames = handles.pcDatasets(planeNum).Info.frames;
            timeres = handles.pcDatasets(planeNum).Info.timeres; %temporal resolution (ms)
        else 
            xres = handles.pcDatasets(planeNum).Info.PixelSpacing(1); %resolution (mm). ASSUMED TO BE SAME IN Y DIMENSION
            rrInterval = handles.pcDatasets(planeNum).Info.NominalInterval; %average RR interval (ms)
            frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages; %number of cardiac frames
            timeres = rrInterval/frames; %temporal resolution (ms)
        end 
        
        area = sum(roiMask(:))*(xres)^2; %ROI area (mm^2)
        for i=1:size(v,3)
            vTemp = v(:,:,i); %through-plane velocity in frame i
            roiDataRaw(:,i) = double(vTemp(roiMask)); %indexed velocities within mask
            meanROI(i) = mean(roiDataRaw(:,i)); %mean velocity in frame i (mm/s)
            maxROI(i) = max(double(vTemp(roiMask))); %max velocity in frame i (mm/s)
            minROI(i) = min(double(vTemp(roiMask))); %min velocity in frame i (mm/s)
            stdROI(i) = std(double(vTemp(roiMask))); %stdv of velocity in frame i (mm/s)
            flowROI(i) = area.*meanROI(i).*0.001; %flow in frame i (mm^3/s = mL/s)
        end 
        
        times = double(timeres.*(1:frames)); %original times
        timesInterp = double(timeres.*tq); %interpolated times
        
        %Linear interpolation (to get more points on flow curve)
        meanROIfit = interp1(times,meanROI,timesInterp,'linear');
        maxROIfit = interp1(times,maxROI,timesInterp,'linear');
        minROIfit = interp1(times,minROI,timesInterp,'linear');
        stdROIfit = interp1(times,stdROI,timesInterp,'linear');
        flowROIfit = interp1(times,flowROI,timesInterp,'linear');
        
        %Add data to roiStatistics structure
        roiStatistics.radius = radius; 
        roiStatistics.center = center;
        roiStatistics.roiMask = roiMask; 
        roiStatistics.roiDataRaw = roiDataRaw;
        roiStatistics.times = timesInterp;
        roiStatistics.meanROI = meanROIfit; 
        roiStatistics.maxROI = maxROIfit; 
        roiStatistics.minROI = minROIfit;
        roiStatistics.stdROI = stdROIfit;
        roiStatistics.flowROI = flowROIfit; 
        roiStatistics.Name = ''; %will get changed below
        roiStatistics.ROInumber = handles.global.totalROIs;

        %%% Create Interpolated Curve with Gaussian Smoothing           
        meanROIfit = interp1(times,smoothdata(meanROI,'gaussian',5),timesInterp,'linear');
        maxROIfit = interp1(times,smoothdata(maxROI,'gaussian',5),timesInterp,'linear');
        minROIfit = interp1(times,smoothdata(minROI,'gaussian',5),timesInterp,'linear');
        stdROIfit = interp1(times,smoothdata(stdROI,'gaussian',5),timesInterp,'linear');
        flowROIfit = interp1(times,smoothdata(flowROI,'gaussian',5),timesInterp,'linear');

        roiStatisticsGaussian.times = timesInterp;
        roiStatisticsGaussian.meanROI = meanROIfit; 
        roiStatisticsGaussian.maxROI = maxROIfit; 
        roiStatisticsGaussian.minROI = minROIfit;
        roiStatisticsGaussian.stdROI = stdROIfit;
        roiStatisticsGaussian.flowROI = flowROIfit; 

        %%% Create Cubic Spline Fit   
        meanROIfit = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
        maxROIfit = csaps([0 times times(end)+times(1)],[0 maxROI 0],0.0001,timesInterp);
        minROIfit = csaps([0 times times(end)+times(1)],[0 minROI 0],0.0001,timesInterp);
        stdROIfit = csaps([0 times times(end)+times(1)],[0 stdROI 0],0.0001,timesInterp);
        flowROIfit = csaps([0 times times(end)+times(1)],[0 flowROI 0],0.0001,timesInterp);

        roiStatisticsSpline.times = timesInterp;
        roiStatisticsSpline.meanROI = meanROIfit; 
        roiStatisticsSpline.maxROI = maxROIfit; 
        roiStatisticsSpline.minROI = minROIfit;
        roiStatisticsSpline.stdROI = stdROIfit;
        roiStatisticsSpline.flowROI = flowROIfit; 

        %%% Save all data into handles
        if isstruct(handles.pcDatasets(planeNum).ROIdata)
            handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
            handles.pcDatasets(planeNum).ROIdataGaussian(end+1) = roiStatisticsGaussian;
            handles.pcDatasets(planeNum).ROIdataSpline(end+1) = roiStatisticsSpline;
        else 
            handles.pcDatasets(planeNum).ROIdata = roiStatistics;
            handles.pcDatasets(planeNum).ROIdataGaussian = roiStatisticsGaussian;
            handles.pcDatasets(planeNum).ROIdataSpline = roiStatisticsSpline;
        end 
        
        dataDir = handles.pcDatasets(planeNum).RootDir; %directory in which plane data is located
        if dataDir(end)=='\' %kill the slash if it exists
            dataDir(end) = [];
        end 
        if ~exist([dataDir '\ROIimages'],'dir') %if the proposed directory doesn't exist
            mkdir([dataDir '\ROIimages']); %make it
            cd([dataDir '\ROIimages']); %move into it
            frame = getframe(handles.pcPlanePlot); %get a snapshot of the PC plane plot with ROI
            image = frame2im(frame); %make into image
            imwrite(image,[handles.pcDatasets(planeNum).Names '.png']) %write it out as PNG
        else
            cd([dataDir '\ROIimages']); %if ROIimages already exists, move into it
            frame = getframe(handles.pcPlanePlot);
            image = frame2im(frame);
            imwrite(image,[handles.pcDatasets(planeNum).Names '.png'])
        end 
        cd(handles.global.homeDir); %lets go back home
    end 
    % else we do analysis for bSSFP
    % will be implemented in the future
    
    %%% Label each ROI w/ names (helpful because there may be 2 ROIs/plane)
    for i=1:numel(handles.pcDatasets)
        if isstruct(handles.pcDatasets(i).ROIdata) %if we've made ROI data for this dataset
            if length(handles.pcDatasets(i).ROIdata)==1
                handles.pcDatasets(i).ROIdata.Name = handles.pcDatasets(planeNum).Names; %name ROI
            else
                for j=1:length(handles.pcDatasets(i).ROIdata)
                    name = handles.pcDatasets(i).Names; %get plane name
                    planeName = [name ' ROI ' num2str(j)]; %needed if more than one ROI/plane
                    handles.pcDatasets(i).ROIdata(j).Name = planeName; %name ROI
                end 
            end 
        end 
    end 

    set(handles.interpolatePopup,'Enable','on'); %turn on interpolate button
    plotVelocity(handles); %plot flow curves
    deleteROIbutton_Callback(hObject, eventdata, handles); %delete ROI info
    
    guidata(hObject,handles);
    updatePCImages(handles); %update images (to remove green ROI circle)
    axes(handles.pcPlanePlot); %make sure we're still on PC plot
    
    
% --- COMPLETE LOADING BUTTON - CALLBACK
function completeLoadingROI_Callback(hObject, eventdata, handles)
    if handles.global.totalROIs>1 %if we've made at least 2 ROIs, proceed
        set(handles.drawROIbutton,'Enable','off'); %disallow interactions with ROIs now
        set(handles.deleteROIbutton,'Enable','off'); 
        set(handles.loadROIbutton,'Enable','off'); 
        set(handles.drawCenterlineButton,'Enable','on'); %allow centerline drawing and loading
        set(handles.loadCenterline,'Enable','on'); 
        set(handles.showPCPlanesRadio,'Enable','on'); %turn on displaying of 2D slices on anatomical images
        set(handles.range1Edit,'Enable','on');
        set(handles.range2Edit,'Enable','on');
        set(handles.anatListbox,'Enable','on');
        set(handles.anatZoomButton,'Enable','on');
        set(handles.anatUnzoomButton,'Enable','on');   

        for i = 1:numel(handles.anatDatasets) %for each anatomical dataset
            images = handles.anatDatasets(i).Data; %grab images from handle
            dims(1) = size(images,1); %get dimensions
            dims(2) = size(images,2);
            dims(3) = size(images,3);
            handles.anatDatasets(i).dims = dims; 

            originShift = [handles.anatDatasets(i).Info.ImagePositionPatient;1]; % origin is top left corner of image
            xres = handles.anatDatasets(i).Info.PixelSpacing(1);
            yres = handles.anatDatasets(i).Info.PixelSpacing(2);
            zres = handles.anatDatasets(i).Info.SliceThickness;

            % sometimes get very small values that should be 0 (so I round)
            xVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(1:3) ,8); % what direction rows run w/r/to x
            yVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(4:6) ,8); % what direction the cols run w/r/to y
            zVector = [cross(xVector,yVector);0]; %take cross product to get zVector
            xVector = [xVector;0]; %add a dummy zero to the end
            yVector = [yVector;0]; 
            rotationMatrix = [xres*xVector yres*yVector zres*zVector originShift]; %create rotation matrix
            sliceDim = find(rotationMatrix(:,3)); %find our slice dimension (eg 2=Y,3=Z)
            colsRunningDir = sign(nonzeros(rotationMatrix(:,3))); %check if cols are running forward (+) or backwards (-)
            handles.anatDatasets(i).rotationMatrix = rotationMatrix;      
            handles.anatDatasets(i).sliceDim = sliceDim;

            topFrontLeft = rotationMatrix*[1;1;1;1]; %topfrontleft coords
            topFrontLeft(4) = []; % remove dummy dimension
            backBottomRight = rotationMatrix*[dims(1);dims(2);dims(3);1];
            backBottomRight(4) = [];
            spanZ(1) = topFrontLeft(3); %first Z point (mm)
            spanZ(2) = backBottomRight(3); %last Z point (mm)
            handles.anatDatasets(i).spanZ = spanZ; %physical span in Z direction in mm (head to toe)
            
            %Below is to calculate the row in which we need to display our PC slice over anatomical image.
            %This is important for drawing centerlines correctly.
            planeLineRows = [];
            for j=1:numel(handles.pcDatasets)
                if isstruct(handles.pcDatasets(j).ROIdata) %if we've made ROI data for this dataset
                    if isfield(handles.pcDatasets(j).Info,'matrixx') %if we're dealing with radial data (pcvipr recon)
                        zLoc = handles.pcDatasets(j).Info.sz; %get physical location of plane in Z (mm above isocenter)
                    else %if we're dealing with dicom data
                        zLoc = handles.pcDatasets(j).Info.ImagePositionPatient(3); %do the same as above
                    end 
                    planeLinePhysical = colsRunningDir*(zLoc-spanZ(1)); %physical distance from top of anat image to PC slice
                    planeLineRows(end+1) = abs(round(planeLinePhysical/yres))+1; %turn distance into row in image (to display PC slice over anatomy)
                    %we add one because we start at an index of 1
                end 
            end 
            handles.anatDatasets(i).planeLineRows = planeLineRows; %save into handles
        end   

        dirSplits = regexp(handles.pcDatasets(1).RootDir,'\','split'); %break up directory in which plane data is located
        dirSplits(end) = []; %chop off last folder
        baseDir = strjoin(dirSplits,'\'); %rejoin string to get name of folder one up from plane data
        if ~exist([baseDir '\FlowWaveformImages'],'dir') %if directory doesn't exist
            mkdir([baseDir '\FlowWaveformImages']); %make it
            cd([baseDir '\FlowWaveformImages']); %go to it
            frame = getframe(handles.velocityPlot); %get snapshot of velocity plot
            imwrite(frame2im(frame),'FlowWaveform.png') %write out to PNG
        else %or if the directory already exists
            cd([baseDir '\FlowWaveformImages']); %go to it
            frame = getframe(handles.velocityPlot);
            imwrite(frame2im(frame),'FlowWaveform.png')
        end 
        cd(handles.global.homeDir); %lets go home now
        guidata(hObject,handles);
    else
        set(handles.messageBar,'String','You need at least 2 ROIs for PWV measurements.');
    end 
    
    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function velocityPlot_CreateFcn(hObject, eventdata, handles)
    

% --- INTERPOLATE POPUP - CALLBACK
function interpolatePopup_Callback(hObject, eventdata, handles)
    interp = get(handles.interpolatePopup,'Value');
    switch interp
        case 1
            handles.global.interpType = 'Gaussian'; %set global flag
        case 2
            handles.global.interpType = 'Spline';
        case 3
            handles.global.interpType = 'None';
        otherwise
    end 
 
    guidata(hObject, handles);
    plotVelocity(handles) %replot our velocity with interpolated data
    
    if handles.global.startAnalyzing %if we're already analyzing PWVs
        computePWVButton_Callback(hObject, eventdata, handles); %recompute PWVs with interpolated data
    end 

% --- INTERPOLATE POPUP - CREATE FUNCTION
function interpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ERROR BAR RADIO - CALLBACK
function errorBarRadio_Callback(hObject, eventdata, handles)
    if get(handles.errorBarRadio,'Value')
        handles.global.showErrorBars = 1; %set global flag
    else 
        handles.global.showErrorBars = 0; 
    end 
    
    guidata(hObject,handles);
    plotVelocity(handles) %replot flow curves with errors bars
    
    
% --- PG SHIFT RADIO - CALLBACK
function pgShiftRadio_Callback(hObject, eventdata, handles)
    if get(handles.pgShiftRadio,'Value')
        handles.global.pgShift = 1; %set global flag
    else 
        handles.global.pgShift = 0;
    end 
    
    guidata(hObject,handles); 
    plotVelocity(handles) %replot flow curves with half cycle shift
    if handles.global.startAnalyzing %if we're already analyzing PWVs
        computePWVButton_Callback(hObject, eventdata, handles); %recompute PWVs with shifted curves
    end 


    
%%%%%%%%%%%% ANATOMICAL PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- ANATOMICAL PLOT - CREATE FUNCTION
function anatPlanePlot_CreateFcn(hObject, eventdata, handles)


% --- ANATOMICAL SLIDER - CALLBACK
function anatSlider_Callback(hObject, eventdata, handles)
    updateAnatImages(handles); %update anatomical images when slider is moved

% --- ANATOMICAL SLIDER - CREATE FUNCTION
function anatSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

 
% --- ANATOMICAL LISTBOX - CALLBACK
function anatListbox_Callback(hObject, eventdata, handles)
    showPCPlanesRadio_Callback(hObject, eventdata, handles) %if anat dataset is changed, replace planes by calling this fcn
    updateAnatImages(handles); %now update images

% --- ANATOMICAL LISTBOX - CREATE FUNCTION
function anatListbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SHOW PLANE RADIO - CALLBACK
function showPCPlanesRadio_Callback(hObject, eventdata, handles)
    guidata(hObject,handles);
    updateAnatImages(handles);
    
    
% --- ZOOM ANATOMICAL - CALLBACK
function anatZoomButton_Callback(hObject, eventdata, handles)
    axes(handles.anatPlanePlot); %make sure we do this on PC plot (not anatomical plot on bottom)
    set(handles.messageBar,'String','Draw an ROI to zoom in');
    rect = drawrectangle; %draw a rectangular field of view
    positions = round(rect.Position); %get positions of rectangle edges
    spanUDAnat = handles.global.spanUDAnat; %get current span of rows
    spanLRAnat = handles.global.spanLRAnat; %get current span of columns
    % Update row span with rectangle vertices (note that we can only zoom in (hence the plus)
    handles.global.spanUDAnat = spanUDAnat(1)+positions(2):(spanUDAnat(1)+positions(2)+positions(4));
    handles.global.spanLRAnat = spanLRAnat(1)+positions(1):(spanLRAnat(1)+positions(1)+positions(3));
    handles.global.zoomedAnat = 1; %set flag to zoomed

    guidata(hObject, handles);
    updateAnatImages(handles);  %update to show zoomed in image
    
    
% --- UNZOOM ANATOMICAL - CALLBACK
function anatUnzoomButton_Callback(hObject, eventdata, handles)
    axes(handles.anatPlanePlot) %make sure we do this on bottom Anatomical plot
    handles.global.zoomedAnat = 0; %set flag to unzoomed
    handles.global.spanUDAnat = 0; %reset spans to 0
    handles.global.spanLRAnat = 0; 
    
    guidata(hObject, handles);
    updateAnatImages(handles); %update to show full size image
    
    
% --- LOAD CENTERLINE - CALLBACK
function loadCenterline_Callback(hObject, eventdata, handles)
    [anatCLdataset, CenterlineData] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the case file (NAME_pwvcasefile_DATE.mat)');
    load([CenterlineData '\' anatCLdataset]); %load CL data
    
    % If we find a dataset in which anatCLdataset.RootDir matches
    % handles.anatDatasets.RootDir, we'll put the CL/distance info there.
    if sum(strcmp(anatCLdataset.RootDir,{handles.anatDatasets.RootDir}))>0
        idx = find(strcmp(anatCLdataset.RootDir,{handles.anatDatasets.RootDir}));
        handles.anatDatasets(idx).Centerline = anatCLdataset.Centerline;
        handles.anatDatasets(idx).Distances = anatCLdataset.Distances;
    else %else, we'll just put it in the first anat dataset
        handles.anatDatasets(1).Centerline = anatCLdataset.Centerline;
        handles.anatDatasets(1).Distances = anatCLdataset.Distances;
    end 
    cd(handles.global.homeDir);
    set(handles.showPCPlanesRadio,'Enable','off'); %force off flags to avoid errors
    set(handles.drawCenterlineButton,'Enable','off');
    set(handles.computePWVButton,'Enable','on');
    set(handles.anatSlider,'Enable','off');
    set(handles.anatListbox,'Enable','off');
    set(handles.anatZoomButton,'Enable','off');
    set(handles.anatUnzoomButton,'Enable','off');   
    set(handles.drawCenterlineButton,'Enable','off');
    set(handles.range1Edit,'Enable','off');
    set(handles.range2Edit,'Enable','off');
    set(handles.drawCenterlineButton,'Enable','off');
    set(handles.deleteCenterlineButton,'Enable','on');

    guidata(hObject,handles);
    updatePCImages(handles);
    updateAnatImages(handles);
    

% --- DRAW CENTERLINE - CALLBACK
function drawCenterlineButton_Callback(hObject, eventdata, handles)
    axes(handles.anatPlanePlot); %go to Anatomical plot
    datasetNum = get(handles.anatListbox,'Value'); %get current anatomical dataset number
    
    % Lets first make sure we have something typed in range1/range2 
    if isempty(get(handles.range1Edit,'String')) || isempty(get(handles.range1Edit,'String')) 
        set(handles.messageBar,'String','Range is required to draw centerlines.')
    elseif str2double(get(handles.range1Edit,'String'))<0 || str2double(get(handles.range1Edit,'String'))>size(handles.anatDatasets(datasetNum).Data,3) 
        set(handles.messageBar,'String','Range of slices needs to be between 1 and the number of slices.')
    elseif str2double(get(handles.range1Edit,'String'))>str2double(get(handles.range1Edit,'String'))
        set(handles.messageBar,'String','Ranges need to increase from left to right.')
    else %good to go
        set(handles.showPCPlanesRadio,'Value',1); %force show PC planes (needed for centerline drawing)
        updateAnatImages(handles); %update images to show yellow PC planes

        set(handles.showPCPlanesRadio,'Enable','off'); %force off flags to avoid errors
        set(handles.drawCenterlineButton,'Enable','off');
        set(handles.computePWVButton,'Enable','on');
        set(handles.anatSlider,'Enable','off');
        set(handles.anatListbox,'Enable','off');
        set(handles.anatZoomButton,'Enable','off');
        set(handles.anatUnzoomButton,'Enable','off');   
        set(handles.range1Edit,'Enable','on');
        set(handles.range2Edit,'Enable','on');

        range1 = str2double(get(handles.range1Edit,'String')); %get first slice for centerline tracing
        range2 = str2double(get(handles.range2Edit,'String')); %last slice for centerline tracing
        
        xres = handles.anatDatasets(datasetNum).Info.PixelSpacing(1); %get resolutions of anat dataset
        yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
        sliceres = handles.anatDatasets(datasetNum).Info.SliceThickness;

        images = handles.anatDatasets(datasetNum).Data(:,:,range1:range2); %create image stack from range1:range2
        x = []; y = []; z = []; %allocate for coordinates of drawn points
        currSlice = range1;
        for i = 1:size(images,3) %go slice-by-slice through each image
            set(handles.sliceNumText,'String',num2str(currSlice)); %show number of current slice
            image = images(:,:,i); %get slice
            if handles.global.zoomedAnat %if we are zoomed in 
               image = image(handles.global.spanUDAnat,handles.global.spanLRAnat); %make it so we span the desired width
            end
            imshow(image,[]); %show single image
            
            for j = 1:numel(handles.anatDatasets(datasetNum).planeLineRows) %plot PC lines
                pcRow = handles.anatDatasets(datasetNum).planeLineRows(j);
                if handles.global.zoomedAnat
                    topRow = min(handles.global.spanUDAnat);
                    bottomRow = max(handles.global.spanUDAnat);
                    leftColumn = min(handles.global.spanLRAnat);
                    rightColumn = max(handles.global.spanLRAnat);
                    if (pcRow > topRow) && (pcRow < bottomRow)
                        pcRow = pcRow - topRow; %shift our row to our current coordinates
                        %image(pcRow,:) = max(max(image))+100;
                        hold on; plot(1:size(image,1),pcRow.*ones(1,size(image,1)),'y'); %plot as yellow line at height of 2D plane
                    end 
                else
                    hold on; plot(1:size(image,1),pcRow.*ones(1,size(image,1)),'y'); %plot as yellow line at height of 2D plane
                end 
                
            end 
            [xTemp,yTemp] = getpts(); %draw points on image along aorta
            if handles.global.zoomedAnat
                xTemp = xTemp + leftColumn;
                yTemp = yTemp + topRow;
            end 
            zTemp = i.*(ones(size(xTemp,1),1)); %add z-coordinates for slice
            x = [x; xTemp]; %append points
            y = [y; yTemp];
            z = [z; zTemp];
            currSlice = currSlice + 1; %move to next slice
        end 
        
        % Draw polyline from sagittal angle
        fig1 = figure; scatter(y,x); %scatter plot of selected points
        for j = 1:numel(handles.anatDatasets(datasetNum).planeLineRows) %plot PC lines
            pcRow = handles.anatDatasets(datasetNum).planeLineRows(j);
            xline(pcRow); %plot each plane location
        end 
        line1 = drawpolyline; %now draw a polyline through these points
        
        % Draw another polyline from orthogonal coronal angle
        fig2 = figure; scatter(y,z); %scatter plot of same selected points
        for j = 1:numel(handles.anatDatasets(datasetNum).planeLineRows) %plot PC lines
            pcRow = handles.anatDatasets(datasetNum).planeLineRows(j);
            xline(pcRow); %plot each plane location
        end 
        line2 = drawpolyline; %draw second polyline through these points
        splinePositions1 = interppolygon(line1.Position,150); %interpolate sagittal polyline (150 points)
        splinePositions2 = interppolygon(line2.Position,150); %interpolate coronal polyline (150 points)
        splineLine(:,1) = splinePositions2(:,2)+range1; %x image coords of centerline (note we start at range1)
        splineLine(:,2) = splinePositions1(:,2); %y coords
        splineLine(:,3) = splinePositions1(:,1); %z coords
        
        fig3 = figure; scatter3(y,x,z+range1); hold on; %make 3D plot of points
        scatter3(splineLine(:,3),splineLine(:,2),splineLine(:,1),'filled'); %plot spline fit over points
        xlabel('x'); ylabel('y'); zlabel('slice');
        x = linspace( min(splineLine(:,2))-50, max(splineLine(:,2))+50, 512); %create space for '2D' planes
        y = linspace( min(splineLine(:,1))-2, max(splineLine(:,1))+2, 512);
        [X,Y] = meshgrid(x,y);
        for j = 1:numel(handles.anatDatasets(datasetNum).planeLineRows)
            pcRow = handles.anatDatasets(datasetNum).planeLineRows(j); %grab row of PC slice
            Z = pcRow.*ones(size(X));
            hold on; surf(Z,X,Y); %plot PC slice on 3D figure
        end 
        hold off;

        distances = zeros(1,length(splineLine)-1); %initialize distance along centerline vector
        for i=1:length(splineLine)-1
            xSquared = (sliceres.*(splineLine(i,1)-splineLine(i+1,1))).^2;
            ySquared = (xres.*(splineLine(i,2)-splineLine(i+1,2))).^2;
            zSquared = (yres.*(splineLine(i,3)-splineLine(i+1,3))).^2;
            distances(i) = sqrt( xSquared + ySquared + zSquared ); %get 3D Euclidean distance between centerline points
        end 
        distances(end+1) = 0; %make the last entry the end of our centerline
        handles.anatDatasets(datasetNum).Distances = distances; %save centerline info to handles
        handles.anatDatasets(datasetNum).Centerline = splineLine;
        
        dirSplits = regexp(handles.anatDatasets(datasetNum).RootDir,'\','split'); %get root directory of anatomical dicoms
        dirSplits(end) = []; dirSplits(end) = []; %pull back one directory (adds an extra blank directory at end of dirSplits)
        baseDir = strjoin(dirSplits,'\'); %rejoin cell, this should be directory of all dicom images
        if ~exist([baseDir '\CenterlineData'],'dir') %if we haven't made this directory yet
            mkdir([baseDir '\CenterlineData']); %make it
            cd([baseDir '\CenterlineData']); %move to it
            frame = getframe(handles.anatPlanePlot); %get the anatomical image in the main GUI
            imwrite(frame2im(frame),'anatomicalSlice.png') %write to PNG
            savefig(fig3,'Centerline3D'); %save our 3D Matlab figure
            view(10,65); %change the view angle
            frame = getframe(fig3); %get a 2D image of the 3D plot
            imwrite(frame2im(frame),'Centerline3D.png'); %write to PNG
            anatCLdataset = handles.anatDatasets(datasetNum);
            anatCLdataset.Data = [];
            save anatCLdataset.mat anatCLdataset
        else %or just move to the directory if we've made it already
            cd([baseDir '\CenterlineData']);
            frame = getframe(handles.anatPlanePlot);
            imwrite(frame2im(frame),'anatomicalSlice.png')
            savefig(fig3,'Centerline3D');
            view(10,65);
            frame = getframe(fig3);
            imwrite(frame2im(frame),'Centerline3D.png');
            anatCLdataset = handles.anatDatasets(datasetNum);
            anatCLdataset.Data = [];
            save anatCLdataset.mat anatCLdataset
        end 
        cd(handles.global.homeDir);
        close(fig2); close(fig1); %close centerline tracing figures (leave 3D plot)
        set(handles.drawCenterlineButton,'Enable','off');
        set(handles.deleteCenterlineButton,'Enable','on');
        
        guidata(hObject,handles);
        updatePCImages(handles);
        updateAnatImages(handles);
    end
    
% --- DELETE CENTERLINE - CALLBACK
function deleteCenterlineButton_Callback(hObject, eventdata, handles)
    cla(handles.anatPlanePlot,'reset'); %reset anatomical plot
    for i=1:numel(handles.anatDatasets)
        if ~isempty(handles.anatDatasets(i).Centerline)
            handles.anatDatasets(i).Distances = []; %remove CL distances
            handles.anatDatasets(i).Centerline = []; %remove CL points
            
            dirSplits = regexp(handles.anatDatasets(i).RootDir,'\','split'); %get root directory of anatomical dicoms
            dirSplits(end) = []; dirSplits(end) = []; %pull back one directory (adds an extra blank directory at end of dirSplits)
            baseDir = strjoin(dirSplits,'\'); %rejoin cell, this should be directory of all dicom images
            cd([baseDir '\CenterlineData']); %move to it
            delete anatCLdataset.mat
        end
    end

    set(handles.drawCenterlineButton,'Enable','on'); %turn back on flags
    set(handles.showPCPlanesRadio,'Enable','on');
    set(handles.anatListbox,'Enable','on');
    set(handles.anatZoomButton,'Enable','on');
    set(handles.anatUnzoomButton,'Enable','on');
    set(handles.anatSlider,'Enable','on');
    set(handles.showPCPlanesRadio,'Value',0);
    set(handles.range1Edit,'Enable','on');
    set(handles.range2Edit,'Enable','on');
    set(handles.range1Edit,'String',''); %reset slice ranges
    set(handles.range2Edit,'String','');

    guidata(hObject,handles);
    updateAnatImages(handles);
    
    
% --- SLICE NUMBER TEXT - CREATE FUNCTION
function sliceNumText_CreateFcn(hObject, eventdata, handles)
%%% Shows current slice number as slider moves
    

% --- RANGE 1 EDIT - CALLBACK
function range1Edit_Callback(hObject, eventdata, handles)
%%% This is a text entry to type in the desired range of slices for
%%% centerline tracing. Range1 represents the first slice.

% --- RANGE 1 EDIT - CREATE FUNCTION
function range1Edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- RANGE 2 EDIT - CALLBACK
function range2Edit_Callback(hObject, eventdata, handles)
%%% This is a text entry to type in the desired range of slices for
%%% centerline tracing. Range2 represents the last slice

% --- RANGE 2 EDIT - CREATE FUNCTION
function range2Edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
   
    
% --- COMPUTE PWV - CALLBACK
function computePWVButton_Callback(hObject, eventdata, handles)
    handles.global.startAnalyzing = 1; %turn on flag to state that we are ready for PWV analysis
    set(handles.drawCenterlineButton,'Enable','off'); %disable drawing more centerlines
    set(handles.deleteCenterlineButton,'Enable','off'); %disable deleting centerline (we're locked in baby)
    %set(handles.computePWVButton,'Enable','off');
    set(handles.showPCPlanesRadio,'Enable','on'); %reenable showing PC planes on anatomical images
    set(handles.anatListbox,'Enable','on'); %turn back on listbox (so we can scroll through images)
    set(handles.anatZoomButton,'Enable','on'); %turn zoom and unzoom back on 
    set(handles.anatUnzoomButton,'Enable','on');
    set(handles.anatSlider,'Enable','on'); %turn back on slider
    
    handles.flow = organizeFlowInfo(handles);

    distances = [handles.anatDatasets.Distances]; %grab centerline and distances from handles
    centerline = [handles.anatDatasets.Centerline];
    for i=1:numel(handles.anatDatasets)
        if ~isempty(handles.anatDatasets(i).Centerline)
            planeRows = handles.anatDatasets(i).planeLineRows; %get rows of PC planes
        end
    end
    
    z = centerline(:,3); %get z coordinates 
    difference = abs(z - planeRows); %subtract z coords from row of PC planes 
    minima = sum(islocalmin(difference),2); %now look for minima for both planes (intersections)
    roiIdxAlongLine = find(minima); %find the indices of the intersection along the CL
    
    totalROIs = handles.global.totalROIs;
    if sum(minima)<totalROIs %if we're still short on CL/PC plane intersections..
        if z(1)<planeRows(1) %did we start our centerline before the plane?
            minima(1,1) = 1; %if so, we can just start at the earliest CL point
        end 
        if z(end)<planeRows(end)<0 %did we end our centerline before reaching the plane?
            minima(end,2) = 1; %if so, we can end at the last CL point
        end 
    end 
    
    if sum(minima)<handles.global.totalROIs %if we're still short, there is a problem
        set(handles.messageBar,'String','Could not detect all centerline/slice intersections. Look at line 863 in AnalyzePWV.m')
    end 
    
    % COMPUTE DISTANCES
    % Note that we have 3 distances if we have 3 ROIs (AscAo-->DescAo,
    % DescAo-->AbdAo, and AscAo-->AbdAo). WE ONLY WANT FORWARD DISTANCES
    allROIs = 1:totalROIs;
    for i=1:totalROIs
        ROIs2compare = allROIs ~= i; %look at other ROIs
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); %indices of other ROIs
        distances2ROIs = [NaN,NaN,NaN]; %initialize distance array
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j); %other ROIs to compare to current
            if iterator>i %prevents looking backwards
                idx1 = roiIdxAlongLine(i); %index of upstream ROI
                idx2 = roiIdxAlongLine(iterator); %index of downstream ROI
                dist = sum(distances(idx1:idx2)); %distance between (mm)
                distances2ROIs(iterator) = dist; %add to distance array
            end 
        end 
        handles.flow(i).distances2ROIs = distances2ROIs;
    end 
    
    % COMPUTE TIME SHIFTS
    flow = computeTTs(handles.flow,handles.global);
    
    % COMPUTE PWVs
    distance = []; TTPoint = []; TTFoot = []; TTUpstroke = []; Xcorr = [];
    for i = 1:numel(flow)
        mask = ~isnan(flow(i).distances2ROIs);
        distance = [distance flow(i).distances2ROIs(mask)]; %make array of all distances
        TTPoint = [TTPoint flow(i).TTPoint(mask)]; %make array of all TTs
        TTFoot = [TTFoot flow(i).TTFoot(mask)];
        TTUpstroke = [TTUpstroke flow(i).TTUpstroke(mask)];
        Xcorr = [Xcorr flow(i).Xcorr(mask)];
    end
    
    numMethods = get(handles.ttpointRadio,'Value')+get(handles.ttfRadio,'Value') ...
        +get(handles.ttuRadio,'Value')+get(handles.xcorrRadio,'Value'); %add all PWV buttons turned on
    numCompares = numel(distance); %get number of time shift methods  
    average = zeros(1,numCompares); %initialize average array
    
    cla(handles.TimeVsDistance,'reset'); %reset PWV plot
    axes(handles.TimeVsDistance); hold on; %force axes to PWV plot
    xlabel('Distance (mm)'); ylabel('Time Shift (ms)'); %label axes
    sz = 30; %create marker sizes (filled in circles) of 30 pixels
    legendSet = {}; %initialize legend cell array
    if get(handles.ttpointRadio,'Value')
        scatter(distance,TTPoint,sz,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]); %orange 
        for i=1:numCompares
            average(i) = average(i)+TTPoint(i); %add all distances for each ROI location
        end 
        legendSet{end+1} = 'TTPoint';
    end 
    if get(handles.ttfRadio,'Value')
        scatter(distance,TTFoot,sz,'filled','MarkerFaceColor',[0.4940 0.1840 0.5560]); %purple 
        for i=1:numCompares
            average(i) = average(i)+TTFoot(i); %keep adding distances
        end 
        legendSet{end+1} = 'TTFoot';
    end 
    if get(handles.ttuRadio,'Value')
        scatter(distance,TTUpstroke,sz,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]); %green
        for i=1:numCompares
            average(i) = average(i)+TTUpstroke(i);
        end 
        legendSet{end+1} = 'TTUpstroke';
    end 
    if get(handles.xcorrRadio,'Value')
        scatter(distance,Xcorr,sz,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]); %red 
        for i=1:numCompares 
            average(i) = average(i)+Xcorr(i);
        end 
        legendSet{end+1} = 'XCorr';
    end 

    average = average/numMethods; %get average TT for each ROI
    scatter(distance,average,40,'black'); %open black circles (size 40 pixels)
    legendSet{end+1} = 'AVERAGE'; %add average to legend
    
    d = min(distance):max(distance);
    if numel(distance)>2 %if we have more than 2 ROIs, need to calculate average
        % Here, we perform linear regression to find best fit line for each
        % time-to (TT) method. Note: lineFit(1)=slope; lineFit(2)=intercept
        
        %TTPoint Method
        linePointFit = polyfit(distance,TTPoint,1); %fit line to TTPoint points
        linePoint = linePointFit(1)*d + linePointFit(2); %calculate line
        
        %TTFoot Method
        lineFootFit = polyfit(distance,TTFoot,1);
        lineFoot = lineFootFit(1)*d + lineFootFit(2);
        
        %TTUpstroke Method
        lineUpstrokeFit = polyfit(distance,TTUpstroke,1);
        lineUpstroke = lineUpstrokeFit(1)*d + lineUpstrokeFit(2);
        
        %Xcorr Method
        lineXcorrFit = polyfit(distance,Xcorr,1);
        lineXcorr = lineXcorrFit(1)*d + lineXcorrFit(2);
        
        %Average TT
        lineAverageFit = polyfit(distance,average,1);
        lineAverage= lineAverageFit(1)*d + lineAverageFit(2);
        
        PWVpoint = 1/linePointFit(1); %PWV = 1/slope (mm/ms = m/s)
        PWVfoot = 1/lineFootFit(1);
        PWVupstroke = 1/lineUpstrokeFit(1);
        PWVxcorr = 1/lineXcorrFit(1);
        PWVaverage = 1/lineAverageFit(1);
        
        hold on; 
        if get(handles.ttpointRadio,'Value') %if our ttpoint button is on, plot average ttp
            plot(d,linePoint,':','LineWidth',0.2,'MarkerFaceColor',[0.8500 0.3250 0.0980]); %orange 
            set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']); %write out PWV value in text field
        end 
        if get(handles.ttfRadio,'Value')
            plot(d,lineFoot,':','LineWidth',0.2,'MarkerFaceColor',[0.4940 0.1840 0.5560]); %purple
            set(handles.ttfData,'String',[num2str(round(PWVfoot,2)) ' m/s']);
        end 
        if get(handles.ttuRadio,'Value')
            plot(d,lineUpstroke,':','LineWidth',0.2,'MarkerFaceColor',[0.4660 0.6740 0.1880]); %green
            set(handles.ttuData,'String',[num2str(round(PWVupstroke,2)) ' m/s']);
        end 
        if get(handles.xcorrRadio,'Value')
            plot(d,lineXcorr,':','LineWidth',0.2,'MarkerFaceColor',[0.6350 0.0780 0.1840]); %red
            set(handles.xcorrData,'String',[num2str(round(PWVxcorr,2)) ' m/s']);
        end 
        plot(d,lineAverage,'-k','LineWidth',0.2);
        legend(legendSet,'Location','southeast');
        hold off;
    
    else %else we just have two ROIs, one PWV measure (eg AscAo --> DescAo)
        PWVpoint = distance/TTPoint; %PWV = dx/dt (no linear regression)
        PWVfoot = distance/TTFoot;
        PWVupstroke = distance/TTUpstroke;
        PWVxcorr = distance/Xcorr;
        PWVaverage = distance/average;
        
        if get(handles.ttpointRadio,'Value') %if our ttpoint button is on, plot single point
            set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']); %write out PWV value in text field
        end 
        if get(handles.ttfRadio,'Value')
            set(handles.ttfData,'String',[num2str(round(PWVfoot,2)) ' m/s']);
        end 
        if get(handles.ttuRadio,'Value')
            set(handles.ttuData,'String',[num2str(round(PWVupstroke,2)) ' m/s']);
        end 
        if get(handles.xcorrRadio,'Value')
            set(handles.xcorrData,'String',[num2str(round(PWVxcorr,2)) ' m/s']);
        end 
        legend(legendSet,'Location','southeast');
    end 
    
    if numel(distance)==1 %if 2 ROIs, we have 1 distance
        set(handles.distanceHeader1,'String',' ROI 1 --> 2'); %set 1st text field
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); %write out distance b/w ROIs
    elseif numel(distance)==3 %if 3 ROIs, we have 3 distances
        set(handles.distanceHeader1,'String',' ROI 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' ROI 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' ROI 2 --> 3');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']);       
    else %if 4 ROIs, we have 6 distances
        set(handles.distanceHeader1,'String',' ROI 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' ROI 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' ROI 1 --> 4');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']); 
        set(handles.distanceHeader4,'String',' ROI 2 --> 3');
        set(handles.distance4,'String',[num2str(round(distance(4),1)) ' mm']);
        set(handles.distanceHeader5,'String',' ROI 2 --> 4');
        set(handles.distance5,'String',[num2str(round(distance(5),1)) ' mm']); 
        set(handles.distanceHeader6,'String',' ROI 3 --> 4');
        set(handles.distance6,'String',[num2str(round(distance(6),1)) ' mm']);         
    end 
    
    if numMethods>0 %if we have at least one ttbutton on
        set(handles.averageData,'String',[num2str(round(PWVaverage,2)) ' m/s']); %set PWV text field
    else %if have not methods selected (all ttbuttons are off)
        set(handles.averageData,'String','0 m/s'); %set PWV text field to 0 m/s
    end
    
    handles.flow = flow;
    guidata(hObject, handles);
    
    
% --- COMPLETE RESET BUTTON - CALLBACK
function resetButton_Callback(hObject, eventdata, handles)
    answer = questdlg('Pressing YES will reset all ROI data. Pressing NO will cancel. Would you like to reset data?', ...
        'Complete Reset', ...
        'YES','NO','NO'); %prompt user with warning
        switch answer %handle response
            case 'YES'
                reset = 1; %if yes, reset flag is on
            case 'NO'
                reset = 0; %if no, will just pass through this function
        end
    
    if reset
        set(handles.ttpointData,'String',' '); %reset PWV text outputs
        set(handles.ttuData,'String',' ');
        set(handles.ttfData,'String',' ');
        set(handles.xcorrData,'String',' ');
        set(handles.averageData,'String',' ');

        set(handles.distanceHeader1,'String',' '); %reset distance text outputs
        set(handles.distance1,'String',' ');
        set(handles.distanceHeader2,'String',' ');
        set(handles.distance2,'String',' ');
        set(handles.distanceHeader3,'String',' ');
        set(handles.distance3,'String',' ');
        set(handles.distanceHeader4,'String',' ');
        set(handles.distance4,'String',' ');
        set(handles.distanceHeader5,'String',' ');
        set(handles.distance5,'String',' ');
        set(handles.distanceHeader6,'String',' ');
        set(handles.distance6,'String',' ');    

        if isfield(handles,'flow')
            handles = rmfield(handles,'flow'); %remove flow struct
        end 
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROI'); %remove ROI fields
        handles.magDatasets = rmfield(handles.magDatasets,'ROI');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdata');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdataGaussian');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdataSpline');
        
        if isfield(handles.anatDatasets,'Centerline') %remove centerline field
            handles.anatDatasets = rmfield(handles.anatDatasets,'Centerline');
        end 
        if isfield(handles.anatDatasets,'Distances') %remove distance field
            handles.anatDatasets = rmfield(handles.anatDatasets,'Distances');
        end 

        % Turn back into empty fields (will get errors that field doesn't exist)
        handles.pcDatasets(1).ROI = [];
        handles.magDatasets(1).ROI = [];
        handles.anatDatasets(1).Centerline = [];
        handles.anatDatasets(1).Distances = [];
        handles.pcDatasets(1).ROIdata = [];
        handles.pcDatasets(1).ROIdataGaussian = [];
        handles.pcDatasets(1).ROIdataSpline = [];
        
        handles.global.zoomed = 0; %set globals back to initial values
        handles.global.zoomedAnat = 0;
        handles.global.spanUD = [0,0];
        handles.global.spanLR = [0,0];
        handles.global.spanUDAnat = [0,0];
        handles.global.spanLRAnat = [0,0];
        handles.global.interpType = 'Gaussian';
        handles.global.showErrorBars = 0;
        handles.global.startAnalyzing = 0;
        handles.global.pgShift = 0;
        handles.global.totalROIs = 0;

        set(handles.drawROIbutton,'Enable','on'); %turn original buttons back to original state
        set(handles.deleteROIbutton,'Enable','on');
        set(handles.loadROIbutton,'Enable','on');
        set(handles.pcPlanePopup,'String',{handles.pcDatasets.Names}); %reset listbox strings (just to be safe)
        set(handles.pcDatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
        set(handles.interpolatePopup,'String',{'Gaussian','Spline','None'});
        set(handles.anatListbox,'String',{handles.anatDatasets.Names});
        set(handles.interpolatePopup,'Enable','off');
        if numel(handles.magDatasets)==0
            set(handles.magRadio,'Enable','off')
        end 
        set(handles.pcPlanePopup,'Value',1);
        set(handles.pcDatasetPopup,'Value',1);
        set(handles.showPCPlanesRadio,'Enable','off');
        set(handles.drawCenterlineButton,'Enable','off');
        set(handles.deleteCenterlineButton,'Enable','off');
        set(handles.computePWVButton,'Enable','off');
        set(handles.showPCPlanesRadio,'Value',0)
        set(handles.showPCPlanesRadio,'Enable','off')
        set(handles.ttpointRadio,'Value',1);
        set(handles.ttuRadio,'Value',1);
        set(handles.ttfRadio,'Value',1);
        set(handles.xcorrRadio,'Value',1);
        set(handles.interpolatePopup,'Value',1);
        set(handles.pgShiftRadio,'Value',0);
        set(handles.messageBar,'String',' ');
        
        cla(handles.velocityPlot,'reset') %reset all plots
        cla(handles.TimeVsDistance,'reset')
        cla(handles.anatPlanePlot,'reset')
        cla(handles.pcPlanePlot,'reset')
        guidata(hObject, handles); %save data
        updatePCImages(handles); %reset PC images
        updateAnatImages(handles); %reset anatomical images
    end 
    

 
%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- EXPORT ANALYSIS - CALLBACK
function exportAnalysisButton_Callback(hObject, eventdata, handles)
    set(handles.messageBar,'String','Exporting data...');
    interpTypes = {'Gaussian','Spline','None'};
    for t=1:length(interpTypes) %export data for each interp type
        handles.global.interpType = interpTypes{t}; %switch interp type
        flow = computeTTs(handles.flow,handles.global); %recalculate flow struct for each interp
        numROIs = numel(flow);
        
        Distances = []; TTPoint = []; TTFoot = []; TTUpstroke = []; Xcorr = []; TTaverage = [];
        for i = 1:numROIs
            mask = ~isnan(flow(i).distances2ROIs); %find which indices to use for PWV
            Distances = [Distances flow(i).distances2ROIs(mask)]; %make array of all distances
            TTPoint = [TTPoint flow(i).TTPoint(mask)]; %make array of all TTs
            TTFoot = [TTFoot flow(i).TTFoot(mask)];
            TTUpstroke = [TTUpstroke flow(i).TTUpstroke(mask)];
            Xcorr = [Xcorr flow(i).Xcorr(mask)];
        end
        TTaverage = mean([TTFoot; TTPoint; TTUpstroke; Xcorr],1); %get average time shift
        
        % Make first column for excel file
        numCompares = numel(Distances); %get number of PWV measurements
        if numCompares==1 %if 2 ROIs, one measurement
            PLANES{1} = [flow(1).Name ' --> ' flow(2).Name]; %get names for Excel
        elseif numCompares==3
            PLANES{1} = [flow(1).Name ' --> ' flow(2).Name];
            PLANES{2} = [flow(1).Name ' --> ' flow(3).Name];
            PLANES{3} = [flow(2).Name ' --> ' flow(3).Name];   
        else
            PLANES{1} = [flow(1).Name ' --> ' flow(2).Name];
            PLANES{2} = [flow(1).Name ' --> ' flow(3).Name];
            PLANES{3} = [flow(1).Name ' --> ' flow(4).Name];
            PLANES{4} = [flow(2).Name ' --> ' flow(3).Name]; 
            PLANES{5} = [flow(2).Name ' --> ' flow(4).Name]; 
            PLANES{6} = [flow(3).Name ' --> ' flow(4).Name];        
        end 
        
        PLANES{numCompares+1} = 'FIT_PWV'; %These rows be for PWV fit parameters
        PLANES{numCompares+2} = 'm (slope)';
        PLANES{numCompares+3} = 'b (y-intercept)';

        
        PWV_Point = NaN(numCompares+3,1); %add dummy rows to match PLANES size
        PWV_Foot = NaN(numCompares+3,1);
        PWV_Upstroke = NaN(numCompares+3,1);
        PWV_Xcorr = NaN(numCompares+3,1);
        PWV_Average = zeros(1,numCompares);
        
        for i=1:numCompares
            PWV_Point(i) = Distances(i)/TTPoint(i); %calculate pointwise PWVs (not fits)
            PWV_Foot(i) = Distances(i)/TTFoot(i);
            PWV_Upstroke(i) = Distances(i)/TTUpstroke(i);
            PWV_Xcorr(i) = Distances(i)/Xcorr(i);
            PWV_Average(i) = (PWV_Point(i)+PWV_Foot(i)+PWV_Upstroke(i)+PWV_Xcorr(i))/4; %get simple average of all PWVs
        end 
        
        
        if numel(Distances)>2  %if we have at least 3 ROIs (more than one data point), we can use linear regression               
            [linePointFit,~] = polyfit(Distances,TTPoint,1); %get linear regression fit for all points 
            PWV_Point(numCompares+1) = 1/linePointFit(1); %calculate PWV (=1/slope)
            PWV_Point(numCompares+2) = linePointFit(1); %get slope
            PWV_Point(numCompares+3) = linePointFit(2); %get y-intercept
            
            [lineFootFit,~] = polyfit(Distances,TTFoot,1);
            PWV_Foot(numCompares+1) = 1/lineFootFit(1);
            PWV_Foot(numCompares+2) = lineFootFit(1);
            PWV_Foot(numCompares+3) = lineFootFit(2);
            
            [lineUpstrokeFit,~] = polyfit(Distances,TTUpstroke,1);
            PWV_Upstroke(numCompares+1) = 1/lineUpstrokeFit(1);
            PWV_Upstroke(numCompares+2) = lineUpstrokeFit(1);
            PWV_Upstroke(numCompares+3) = lineUpstrokeFit(2);
            
            [lineXcorrFit,~] = polyfit(Distances,Xcorr,1);
            PWV_Xcorr(numCompares+1) = 1/lineXcorrFit(1);
            PWV_Xcorr(numCompares+2) = lineXcorrFit(1);
            PWV_Xcorr(numCompares+3) = lineXcorrFit(2);
            
            [lineAverageFit,~] = polyfit(Distances,TTaverage,1);
            PWV_Average(numCompares+1) = 1/lineAverageFit(1);
            PWV_Average(numCompares+2) = lineAverageFit(1);
            PWV_Average(numCompares+3) = lineAverageFit(2);
            PWV_Average = PWV_Average';
        else %if we have only 2 ROIs (one data point), we can't do linear regression
            PWV_Point(numCompares+1) = Distances/TTPoint; %do a simple division (should be same as PWV_Point(1))
            PWV_Point(numCompares+2) = NaN;
            PWV_Point(numCompares+3) = NaN;
            
            PWV_Foot(numCompares+1) = Distances/TTFoot;
            PWV_Foot(numCompares+2) = NaN;
            PWV_Foot(numCompares+3) = NaN;

            PWV_Upstroke(numCompares+1) = Distances/TTUpstroke;
            PWV_Upstroke(numCompares+2) = NaN;
            PWV_Upstroke(numCompares+3) = NaN;

            PWV_Xcorr(numCompares+1) = Distances/Xcorr;
            PWV_Xcorr(numCompares+2) = NaN;
            PWV_Xcorr(numCompares+3) = NaN;
            
            PWV_Average(numCompares+1) = Distances/TTaverage;
            PWV_Average(numCompares+2) = NaN;
            PWV_Average(numCompares+3) = NaN;
            PWV_Average = PWV_Average'; %flip to match other PWV arrays
        end 
        
        for n = 1:3 %add dummy rows here, couldn't do it above to calculate PWV fits
            Distances(numCompares+n) = NaN;
            TTPoint(numCompares+n) = NaN;
            TTFoot(numCompares+n) = NaN;
            TTUpstroke(numCompares+n) = NaN;
            Xcorr(numCompares+n) = NaN;
            TTaverage(numCompares+n) = NaN;
        end 
        
        % Invert so data is running down (each array will be a column)f
        PLANES = PLANES';
        Distances = Distances';
        TTPoint = TTPoint';
        TTFoot = TTFoot';
        TTUpstroke = TTUpstroke';
        Xcorr = Xcorr';
        TTaverage = TTaverage';
        
        % Make table for writing excel file
        pwvTable = table(PLANES,Distances,TTPoint,TTFoot,TTUpstroke,Xcorr,PWV_Point,PWV_Foot,PWV_Upstroke,PWV_Xcorr,PWV_Average);
        dirSplits = regexp(handles.pcDatasets(1).RootDir,'\','split'); %break up directory in which plane data is located
        dirSplits(end) = []; %chop off last folder
        baseDir = strjoin(dirSplits,'\'); %rejoin string to get name of folder one up from plane data
        date = datestr(now); %get current date/time
        chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)]; %chop date up
        if ~exist([baseDir '\DataAnalysis'],'dir') %if directory doesn't exist
            mkdir([baseDir '\DataAnalysis']); %make it
            cd([baseDir '\DataAnalysis']); %go to it
            writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet','Sheet',['Interpolation - ' interpTypes{t}]); %write excel sheet for each interp
            flow = handles.flow; %make variables for saving
            save('flow.mat','flow')
            save('pwvTable.mat','pwvTable')
        else %or if the directory already exists
            cd([baseDir '\DataAnalysis']); %go to it
            writetable(pwvTable,['Summary_' chopDate '.xlsx'],'FileType','spreadsheet','Sheet',['Interpolation - ' interpTypes{t}]);
            flow = handles.flow;
            save('flow.mat','flow')
            save('pwvTable.mat','pwvTable')
        end 
        cd(handles.global.homeDir); %go back home  
        clear PLANES %need to do this because PLANES will keep getting transposed
    end 
    
    cd([baseDir '\DataAnalysis']); %go to it
    frame = getframe(handles.TimeVsDistance); %get snapshot of PWV plot 
    imwrite(frame2im(frame),'PWVanalysisPlot.png'); %write out to PNG
    anatDatasets = handles.anatDatasets;
    pcDatasets = handles.pcDatasets;
    magDatasets = handles.magDatasets;
    save('anatDatasets.mat','anatDatasets'); %save each dataset for future reference (might be large)
    save('pcDatasets.mat','pcDatasets');
    save('magDatasets.mat','magDatasets');
    cd(handles.global.homeDir); %go back home  

    set(handles.messageBar,'String','Data exported successfully!');
    guidata(hObject, handles);
    
% --- MESSAGE BAR - CALLBACK
function messageBar_Callback(hObject, eventdata, handles)

% --- MESSAGE BAR - CREATE FUNCTION
function messageBar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
   

% --- TTPoint READOUT - CREATE FUNCTION
function ttpointData_CreateFcn(hObject, eventdata, handles)
% --- TTPoint RADIO - CALLBACK
function ttpointRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttpointRadio,'Value') %if we're turned off
        set(handles.ttpointData,'String',' '); %don't display PWV
    end 
    if handles.global.startAnalyzing %if we're analyzing PWVs
        computePWVButton_Callback(hObject, eventdata, handles); %reanalyze without TTpoint 
    end 

% --- TTUpstroke READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)
% --- TTUpstroke RADIO - CALLBACK
function ttuRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',' ');
    end 
    if handles.global.startAnalyzing
        computePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTFoot READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)
% --- TTFoot RADIO - CALLBACK
function ttfRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',' ');
    end 
    if handles.global.startAnalyzing
        computePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- Xcorr READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)
% --- Xcorr RADIO - CALLBACK
function xcorrRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',' ');
    end 
    if handles.global.startAnalyzing
        computePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- AVERAGE READOUT - CREATE FUNCTION
function averageData_CreateFcn(hObject, eventdata, handles)


% --- DISTANCE 1 PlANES - CREATE FUNCTION
function distanceHeader1_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 1 - CREATE FUNCTION
function distance1_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 2 PlANES - CREATE FUNCTION
function distanceHeader2_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 2 - CREATE FUNCTION
function distance2_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 3 PlANES - CREATE FUNCTION
function distanceHeader3_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 3 - CREATE FUNCTION
function distance3_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 4 PlANES - CREATE FUNCTION
function distanceHeader4_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 4 - CREATE FUNCTION
function distance4_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 5 PlANES - CREATE FUNCTION
function distanceHeader5_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 5 - CREATE FUNCTION
function distance5_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 6 PlANES - CREATE FUNCTION
function distance6_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 6 - CREATE FUNCTION
function distanceHeader6_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% --- Update Images in PLANE PLOT
function updatePCImages(handles)
    axes(handles.pcPlanePlot); %force axes to PC plot
    
    planeNum = get(handles.pcPlanePopup,'Value'); % get current plane index (eg AAo)
    datasetNum = get(handles.pcDatasetPopup,'Value'); %get current dataset (eg MAG)

    if get(handles.pcRadio,'Value') %if this is a PC image
        dataset = handles.pcDatasets;
    else %or if this is bSSFP image
        dataset = handles.magDatasets;
    end 

    if isstruct(dataset(planeNum).Data) %if Data is a structure
        imageSet = struct2cell(dataset(planeNum).Data); %convert from struct to cell
    else  %or if we already have an array
        imageSet = dataset(planeNum).Data; %just change the name
    end 

    if iscell(imageSet) %if our set of images are contained in a cell
        images = imageSet(datasetNum); %pull images for current dataset
        images = cell2mat(images); %turn to matrix
    else 
        images = imageSet; %otherwise, do nothing
    end 

    if ndims(images)<3 %if we are dealing with time-averaged images (ndim=2)
        maxSize = max(size(images,1),size(images,2));
        steps = [1 maxSize]; %set our slider as wide as possible so we can't slide
        set(handles.pcSlider,'SliderStep', steps);
        slice = images; %change name for consistency (see below)
    else %if we are dealing with time-resolved images (ndim=3)
        dim3size = size(images,3); %get the size of the third dimension
        steps = [1/(dim3size-1) 10/(dim3size-1)]; %set so one 'slide' moves to the next slice exactly
        set(handles.pcSlider,'SliderStep', steps);
        sliceNum = 1+round( get(handles.pcSlider,'Value').*(dim3size-1) ); %get slice number from slider
        slice = images(:,:,sliceNum); %pull slice from images
    end 

    if handles.global.zoomed %if we are zoomed in 
        slice = slice(handles.global.spanUD,handles.global.spanLR); %make it so we span the desired width
    end 

    if ~isempty(handles.pcDatasets(planeNum).ROI) %if we have an ROI placed
        hold on %make it so we can slide while keeping the ROI on the figure
        imshow(slice,[]);
    else 
        cla(handles.pcPlanePlot,'reset') %otherwise, reset the plot
        imshow(slice,[]) %then show the image
    end 
    
    
% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    axes(handles.anatPlanePlot); %force axes to anatomical plot
    
    datasetNum = get(handles.anatListbox,'Value'); %get our current dataset (eg Sagittal)
    images = handles.anatDatasets(datasetNum).Data; %pull our images from handles

    if ndims(images)<3 %if we only have one image (should rarely happen)
        maxSize = max(size(images,1),size(images,2));
        steps = [1 maxSize]; %set our slider to as wide as possible 
        set(handles.anatSlider,'SliderStep', steps);
        anatSlice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.anatSlider,'SliderStep', steps);
        sliceNum = 1+round( get(handles.anatSlider,'Value').*(dim3size-1) );  %get slice number from slider
        anatSlice = images(:,:,sliceNum); %get slice from images
    end 
    
    if handles.global.zoomedAnat %if we are zoomed in 
        anatSlice = anatSlice(handles.global.spanUDAnat,handles.global.spanLRAnat); %make it so we span the desired width
    end
    
    % Show PC planes overlaid on image (if showPCPlanesRadio is on and
    % we're not looking at an axial slice.
    if get(handles.showPCPlanesRadio,'Value') && handles.anatDatasets(datasetNum).sliceDim~=3
        imshow(anatSlice,[]); %show image
        for i=1:numel(handles.anatDatasets(datasetNum).planeLineRows) %for each 2D plane 
            pcRow = handles.anatDatasets(datasetNum).planeLineRows(i); %get the row of the plane
            %anatSlice(pcRow,:) = max(max(anatSlice))+100; %arbitrarily assign line to high image value
            topRow = min(handles.global.spanUDAnat);
            bottomRow = max(handles.global.spanUDAnat);
            
            if (pcRow > bottomRow) && (pcRow < topRow) && handles.global.zoomedAnat
                pcRow = pcRow - topRow; %shift our row to our current coordinates
                hold on; plot(1:size(anatSlice,1),pcRow.*ones(1,size(anatSlice,1)),'y'); %plot as yellow line at height of 2D plane
            end 
            if ~handles.global.zoomedAnat
               hold on; plot(1:size(anatSlice,1),pcRow.*ones(1,size(anatSlice,1)),'y'); %plot as yellow line at height of 2D plane
            end
            
        end 
    else
        %cla(handles.anatPlanePlot,'reset') %otherwise, reset the plot
        imshow(anatSlice,[]); %else we just show the plain ol image
    end 
    set(handles.sliceNumText,'String',sliceNum); %show our current slice number

    
% --- "Time to" calculations (TTPoint, TTUpstroke, TTFoot, Xcorr)
function flow = computeTTs(flow,globals)
%%% See the following article by Anas Dogui in JMRI:
% Measurement of Aortic Arch Pulse Wave Velocity in Cardiovascular MR:
% Comparison of Transit Time Estimators and Description of a New Approach

    numROIs = globals.totalROIs;
    for i=1:numROIs %for each ROI
        if globals.startAnalyzing %if we're analyzing (put here because we call computePWV in interpolatePopup)
            switch globals.interpType  %find the appropriate interp data
                case 'None'
                    flowTemp = flow(i).Data.meanROI;
                case 'Gaussian'
                    flowTemp = flow(i).Gaussian.meanROI;  
                case 'Spline'
                    flowTemp = flow(i).Spline.meanROI; 
                otherwise
            end 
        else %else, we just use uninterpolated data
            flowTemp = flow(i).Data.meanROI;
        end 
        
        if mean(flowTemp)<0 %if our flow curve is mainly negative
            flowTemp = -1*flowTemp; %flip it upside down so we can find time shifts
        end 
        
        if globals.pgShift %if we want to shift our waveform over
            flowTemp = circshift(flowTemp,round(length(flowTemp)/2)); %shift by half cycle
        end 
                
        times = flow(i).Data.times; %get time frames (ms)
        timeres = times(2)-times(1); %temporal resolution (ms)
        [maxPeakVel,maxPeakVelIdx] = max(flowTemp); %find max velocity value and its location
        upstroke = flowTemp(1:maxPeakVelIdx); %define 'upstroke' region of the flow curve
        curvePoints(i).maxPeakVelIdx = maxPeakVelIdx; %add max point to curvePoints struct
        curvePoints(i).maxPeakVel = maxPeakVel; %add max velocity to curvePoints struct
        
        [~,SeventyPointIdx] = min(abs(upstroke-0.7*maxPeakVel)); %get point at 70% max peak
        curvePoints(i).SeventyPointIdx = SeventyPointIdx; %add to curvePoints struct
        curvePoints(i).SeventyPoint = flowTemp(SeventyPointIdx); %add 70% flow value to curvePoints
        
        [~,ThirtyPointIdx] = min(abs(upstroke-0.3*maxPeakVel)); %get point at 30% max peak
        curvePoints(i).ThirtyPointIdx = ThirtyPointIdx; %add to curvePoints struct
        curvePoints(i).ThirtyPoint = flowTemp(ThirtyPointIdx);
        
        [~,FiftyPointIdx] = min(abs(upstroke-0.5*maxPeakVel)); %get point at 50% max peak    
        curvePoints(i).FiftyPointIdx = FiftyPointIdx;
        curvePoints(i).FiftyPoint = flowTemp(FiftyPointIdx);

        flows(i,:) = normalize(flowTemp,'range'); %normalize curve from here on
    end      
    
    % TTPoint - time to point calculation
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i; %look at other ROIs
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); %indices of other ROIs
        TTPoint = [NaN,NaN,NaN]; %initialize distance array
        for j=1:length(ROIs2compare) 
            iterator = ROIs2compare(j); %other ROIs to compare to current
            if iterator>i %prevents looking backwards
                TTPoint(iterator) = timeres.* ...
                    ( curvePoints(iterator).FiftyPointIdx - curvePoints(i).FiftyPointIdx ); %get time difference b/w 50%'s (ms) 
            end 
        end
        flow(i).TTPoint = TTPoint; %add to flow struct
    end    
    
    % TTUpstroke - time to upstroke calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTUpstroke = [NaN,NaN,NaN];  
        
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                [sigmoid1,upslope1,tUpstroke1,timeresSigmoid1] = sigFit(flows(i,:)); %see sigFit function below
                [sigmoid2,upslope2,tUpstroke2,timeresSigmoid2] = sigFit(flows(iterator,:));
                TTUpstroke(iterator) = timeres.* ...
                    (timeresSigmoid1*tUpstroke2 - timeresSigmoid2*tUpstroke1); %get time diff b/w sigmoid upstrokes
            end 
        end 
         flow(i).TTUpstroke = TTUpstroke;
    end   
    
    % TTFoot - time to foot calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTFoot = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                % m = (y2-y1)/(x2-x1); y1=70%max flow and y2=30%max flow
                % x1 and x2 are the times (indices) where these values occur
                m1 = (curvePoints(i).SeventyPoint - curvePoints(i).ThirtyPoint) ...
                    /(curvePoints(i).SeventyPointIdx - curvePoints(i).ThirtyPointIdx);
                % t1 is the x-intercept of this line. x0 = x1-y1/m
                t1 = curvePoints(i).ThirtyPointIdx - (curvePoints(i).ThirtyPoint/m1);
                
                % Do the same for the consequent flow curve
                m2 = (curvePoints(iterator).SeventyPoint - curvePoints(iterator).ThirtyPoint) ...
                    /(curvePoints(iterator).SeventyPointIdx - curvePoints(iterator).ThirtyPointIdx);
                t2 = curvePoints(iterator).ThirtyPointIdx - (curvePoints(iterator).ThirtyPoint/m2);
                TTFoot(iterator) = timeres.*(t2-t1); %get time difference b/w x-intercepts of fits
            end 
        end 
        flow(i).TTFoot = TTFoot;
    end 
    
    % XCorr - cross correlation calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        Xcorr = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            flows(i) = normalize(flows(i),'range');
            iterator = ROIs2compare(j);
            if iterator>i
                [Xcorrs,lags] = xcorr(flows(i,:),flows(iterator,:)); %perform cross correlation between flow curves
                [~,maxXcorrIdx] = max(Xcorrs); %get index of max Xcorr value
                shift = -1*lags(maxXcorrIdx); %find time lag of Xcorr peak
                Xcorr(iterator) = shift.*timeres; %get time shift
            end 
        end 
        flow(i).Xcorr = Xcorr;
    end 
    
 
% --- Turn PolyLine into SplineLine    
function Y = interppolygon(X,N)
    if nargin < 2 || N < 2
        N = 2; %if only one arg or too small N, just assume 2
    end
    nDim = size(X,2); %should be 2
    dx = 0;
    
    for dim = 1:nDim
        dx = dx + diff(X(:,dim)).^2 ; %get sum of squares in each dim
    end
    
    lengthBetweenPoints = sqrt(dx); %now get distance
    lengthLine = sum(lengthBetweenPoints);
    origMetric = [0; cumsum(lengthBetweenPoints/lengthLine)];
    
    interpMetric = (0:(1/(N-1)):1)';
    Y = interp1(origMetric,X,interpMetric,'makima'); %makima seems to work well
    %Y = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
    

% --- Plot Velocities    
function plotVelocity(handles)
    cla(handles.velocityPlot,'reset'); %reset axes
    axes(handles.velocityPlot); %make sure we plot on the right axis
    
    times = handles.pcDatasets(1).ROIdata(1).times;
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5); %line of y=0 (for visual reference)
    % Note that above we assume same time scale for each plane (MR scan)
    
    legendSet = {'Baseline'}; %add baseline to legend names
    for i=1:numel(handles.pcDatasets) %for all planes
        if isstruct(handles.pcDatasets(i).ROIdata) %if we've made ROI data for this dataset
            for j=1:length(handles.pcDatasets(i).ROIdata) %for each ROI
                legendSet{end+1} = handles.pcDatasets(i).ROIdata(j).Name; %add name of ROI to list
                switch handles.global.interpType %check what interpolation we are using
                    case 'None'
                        velocity = handles.pcDatasets(i).ROIdata(j).meanROI; %grab mean velocity
                        stdv = handles.pcDatasets(i).ROIdata(j).stdROI; %grab stdv of velocity
                    case 'Gaussian'
                        velocity = handles.pcDatasets(i).ROIdataGaussian(j).meanROI;
                        stdv = handles.pcDatasets(i).ROIdataGaussian(j).stdROI;
                    case 'Spline'
                        velocity = handles.pcDatasets(i).ROIdataSpline(j).meanROI;
                        stdv = handles.pcDatasets(i).ROIdataSpline(j).stdROI;
                    otherwise
                end 
                
                if handles.global.pgShift %check if we need to shift waveform
                    velocity = circshift(velocity,round(length(velocity)/2)); %shift by half cycle
                    stdv = circshift(stdv,round(length(stdv)/2));
                end 

                if mean(velocity)<0 %if we are mainly negative velocities (as in descending aorta)
                    if handles.global.showErrorBars %and if we want to show error bars
                        hold on; errorbar(times,-1*velocity,stdv); %plot inverted velocity and errors
                    else
                        hold on; plot(times,-1*velocity); %else, plot inverted velocity
                    end 
                else %otherwise, don't invert velocity (as in ascending aorta)
                    if handles.global.showErrorBars
                        hold on; errorbar(times,velocity,stdv);
                    else
                        hold on; plot(times,velocity);
                    end
                end 
            end 
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)'); %set axes labels
    %xlim([times(1),times(end)]); %chop limits to make curve full width
    

% --- Sigmoid Fit Function
function [sigmoid,upslope,t1,dt] = sigFit(meanROI)
%%% See the following article by Anas Dogui in JMRI:
% Measurement of Aortic Arch Pulse Wave Velocity in Cardiovascular MR:
% Comparison of Transit Time Estimators and Description of a New Approach

    [~,peak] = max(meanROI); %find max
    upslope = meanROI(1:peak); %find upslope region of flow curve
    t0 = 1:peak; %get times for upslope curve
    t = 1:0.1:peak; %interpolate even more
    upslope = interp1(t0,upslope,t); %interpolate upslope
    upslope = normalize(upslope,'range'); %normalize from 0 to 1
    dt = t(2)-t(1); %new temporal resolution (=0.1)
    
    % c1 = b, c2 = a, c3 = x0, c4 = dx
    % Note that we could assume the equation e^t/(1+e^(t-t0)) since c1=1 and
    % c2=0. However, will keep the same as the Dogui paper.
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ) - upslope;
    c0 = [max(upslope),min(upslope),peak/2,dt/2]; %initial params for upslope region
    opts = optimset('Display', 'off'); %turn off display output
    c = lsqnonlin(sigmoidModel,c0,[],[],opts); %get nonlinear LSQ solution
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ); %calculate our sigmoid fit with our new params

    % Could use diff function here, but diff removes one element from array
    for i=2:numel(upslope)-1
        dt = 0.5*(t(i+1)-t(i-1)); %first time derivative
        ddt = 4*(t(i+1)-2*t(i)+t(i-1)); %second time derivative
        dy = 0.5*(sigmoid(i+1)-sigmoid(i-1)); %first y derivative
        ddy = 4*(sigmoid(i+1)-2*sigmoid(i)+sigmoid(i-1)); %second y derivative
        curvature(i) = ddy*dt ./ ((dt^2 + dy^2).^(3/2)); %curvature function
    end 
    [~,t1] = max(curvature);
    %save upslope sigmoid curvature plots

% --- Condense and organize flow data obtained from ROIs
function flow = organizeFlowInfo(handles)
% This function is designed to pull apart handles.pcDatasets. It is
% difficult to do analysis on the handles structure because some slices
% have 2 ROIs. It is much easier to have a structure that pulls out each
% ROI. It only adds a bit of memory since we aren't saving the raw images.
    count = 1; %overall iterator
    for i=1:numel(handles.pcDatasets) %for each PC dataset
        if isstruct(handles.pcDatasets(i).ROIdata) %do we have data?
            if length(handles.pcDatasets(i).ROIdata)==1 %if we just have one ROI
                flow(count).Name = handles.pcDatasets(i).ROIdata.Name; %pull only relevant info for PWV calcs
                flow(count).Data = handles.pcDatasets(i).ROIdata;
                flow(count).Info = handles.pcDatasets(i).Info;
                flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian;
                flow(count).Spline = handles.pcDatasets(i).ROIdataSpline;
                flow(count).pcDatasetREF = [i,1];
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdata) %if we have more than one ROI
                    flow(count).Name = handles.pcDatasets(i).ROIdata(j).Name; %parse into individual ROIs in flow struct
                    flow(count).Data = handles.pcDatasets(i).ROIdata(j);
                    flow(count).Info = handles.pcDatasets(i).Info;
                    flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian(j);
                    flow(count).Spline = handles.pcDatasets(i).ROIdataSpline(j);
                    flow(count).pcDatasetREF = [i,j];
                    count = count+1;
                end 
            end 
        end 
    end 
  
