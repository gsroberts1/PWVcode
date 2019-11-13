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

    % Last Modified by GUIDE v2.5 05-Nov-2019 11:15:17

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

    
    handles.anatDatasets = varargin{1};
    handles.pcDatasets = varargin{2};
    handles.magDatasets = varargin{3};
    handles.pcDatasets(1).ROI = [];
    handles.magDatasets(1).ROI = [];
    handles.anatDatasets(1).Centerline = [];
    handles.pcDatasets(1).ROIdata = [];
    handles.pcDatasets(1).ROIdataMakima = [];
    handles.pcDatasets(1).ROIdataGaussian = [];
    handles.pcDatasets(1).ROIdataSpline = [];
    
    handles.global.zoomed = 0;
    handles.global.zoomedAnat = 0;
    handles.global.spanUD = [0,0];
    handles.global.spanLR = [0,0];
    handles.global.spanUDAnat = [0,0];
    handles.global.spanLRAnat = [0,0];
    handles.global.interpType = 'None';
    handles.global.pgShift = 0;
    handles.global.showErrorBars = 0;
    handles.global.startAnalyzing = 0;
    
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.InterpolatePopup,'String',{'None','Makima','Gaussian','Spline'});
    set(handles.AnatListbox,'String',{handles.anatDatasets.Names});
    set(handles.InterpolatePopup,'Enable','off');
    if numel(handles.magDatasets)==0
        set(handles.MAGRadio,'Enable','off')
    end 
    
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.DrawCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ComputePWVButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Enable','off')
    set(handles.ttpRadio,'Value',1);
    set(handles.ttpointRadio,'Value',1);
    set(handles.ttuRadio,'Value',1);
    set(handles.ttfRadio,'Value',1);
    set(handles.xcorrRadio,'Value',1);
        
    guidata(hObject, handles);
    updateImages(handles);
    updateAnatImages(handles);

% --- Outputs from this function are returned to the command line.
function varargout = AnalyzePWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    
   
    
%%%%%%%%%%%% LOAD 2DPC PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PLANE PLOT - CREATE FUNCTION
function PlanePlot_CreateFcn(hObject, eventdata, handles)


% --- PLANE SLIDER - CALLBACK
function PlaneSlider_Callback(hObject, eventdata, handles)
    updateImages(handles);

% --- PLANE SLIDER - CREATE FUNCTION
function PlaneSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

    
% --- PC RADIO - CALLBACK
function PCRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'Enable','on');
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.MAGRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    updateImages(handles);

    
% --- MAG RADIO - CALLBACK
function MAGRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.magDatasets.Names});
    set(handles.DatasetPopup,'Enable','off');
    set(handles.PCRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    set(handles.DatasetPopup,'String',' ');
    updateImages(handles);
    

% --- PLANE DROPDOWN - CALLBACK
function PlanePopup_Callback(hObject, eventdata, handles)   
    handles.global.zoomed = 0;
    handles.global.spanUD = 0;
    handles.global.spanLR = 0;
    
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(planeNum).Data));
    end 
    set(handles.DatasetPopup,'Value',1);
    
    clear planeNum
    guidata(hObject, handles);
    updateImages(handles);

% --- PLANE DROPDOWN - CREATE FUNCTION
function PlanePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- DATASET DROPDOWN - CALLBACK
function DatasetPopup_Callback(hObject, eventdata, handles)
    updateImages(handles);

% --- DATASET DROPDOWN - CREATE FUNCTION
function DatasetPopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ZOOM BUTTON - CALLBACK
function ZoomButton_Callback(hObject, eventdata, handles)
    axes(handles.PlanePlot)
    
    set(handles.MessageBar,'String','Draw an ROI to zoom in.');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUD = handles.global.spanUD;
    spanLR = handles.global.spanLR;
    handles.global.spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4));
    handles.global.spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));
    handles.global.zoomed = 1;
    
    clear rect positions spanUD spanLR
    guidata(hObject, handles);
    updateImages(handles);
 
    
% --- UNZOOM BUTTON - CALLBACK
function UnzoomButton_Callback(hObject, eventdata, handles)
    axes(handles.PlanePlot)
    
    handles.global.zoomed = 0; 
    handles.global.spanUD = 0; 
    handles.global.spanLR = 0;
    
    guidata(hObject,handles);
    updateImages(handles);

    
% --- DRAWROI BUTTON - CALLBACK
function DrawROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.PlanePlot);
    
    if get(handles.PCRadio,'Value')
        set(handles.MAGRadio,'Enable','off');
    else 
        set(handles.PCRadio,'Enable','off');
    end 
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.PlanePopup,'Enable','off');
    set(handles.ZoomButton,'Enable','off');
    set(handles.UnzoomButton,'Enable','off');
    
    planeNum = get(handles.PlanePopup,'Value');
    if isempty(handles.pcDatasets(planeNum).ROI)
        mydlg = warndlg('Press enter when the ROI is set');
        waitfor(mydlg);
        circle = drawcircle('FaceAlpha',0.1,'Color','g','LineWidth',1,'Deletable',0);
        while true
            w = waitforbuttonpress; 
            switch w 
                case 1 % (keyboard press) 
                  key = get(gcf,'currentcharacter'); 
                      switch key
                          case 27 % 27 is the escape key
                              set(handles.MessageBar,'String','User pressed the escape key. Deleting ROI.');
                              DeleteROIbutton_Callback(hObject, eventdata, handles)
                              break % break out of the while loop
                          case 13 % 13 is the return key 
                              set(handles.MessageBar,'String','ROI selected.');
                              circle.InteractionsAllowed = 'none';
                              break
                          otherwise 
                              % Wait for a different command. 
                      end
           end
        end

    if get(handles.PCRadio,'Value')
        handles.pcDatasets(planeNum).ROI = circle;
    else 
        handles.magDatasets(planeNum).ROI = circle;
    end 
    
    clear planeNum mydlg w key
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('An ROI has already been placed!\n');
    end 

    
% --- DELETE ROI BUTTON - CALLBACK
function DeleteROIbutton_Callback(hObject, eventdata, handles)
   axes(handles.PlanePlot)
    hold off
    
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        handles.pcDatasets(planeNum).ROI = [];
    else 
        handles.magDatasets(planeNum).ROI = [];
    end 

    set(handles.DrawROIbutton,'Enable','on');
    set(handles.PlanePopup,'Enable','on');
        if get(handles.PCRadio,'Value')
            set(handles.MAGRadio,'Enable','on');
        else 
            set(handles.PCRadio,'Enable','on');
        end 
    set(handles.ZoomButton,'Enable','on');
    set(handles.UnzoomButton,'Enable','on');

    clear planeNum 
    guidata(hObject,handles);
    updateImages(handles);
    
% --- LOAD ROI BUTTON - CALLBACK
function LoadROIbutton_Callback(hObject, eventdata, handles)
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        v = handles.pcDatasets(planeNum).Data.v;
        if handles.global.zoomed
            v = v(handles.global.spanUD,handles.global.spanLR,:);
        end 
        circle = handles.pcDatasets(planeNum).ROI;
        radius = circle.Radius; 
        center = round(circle.Center);

        [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
        X = X-center(2); %shift coordinate grid
        Y = Y-center(1);
        roiMask = sqrt(X.^2+Y.^2)<=radius;

        %%% Create Linear Interpolated Data
        tq = 1:0.1:size(v,3);
        
        if isfield(handles.pcDatasets(planeNum).Info,'matrixx')
            matrixx = handles.pcDatasets(planeNum).Info.matrixx;
            fovx = handles.pcDatasets(planeNum).Info.fovx;
            xres = fovx/matrixx;
            timeres = handles.pcDatasets(planeNum).Info.timeres;
        else 
            xres = handles.pcDatasets(planeNum).Info.PixelSpacing(1);
            bpm = handles.pcDatasets(planeNum).Info.HeartRate;
            frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages;
            rrInterval = (60*1000)/bpm;
            timeres = rrInterval/frames;
        end 
        
        for i=1:size(v,3)
            vTemp = v(:,:,i);
            roiDataRaw(:,i) = double(vTemp(roiMask));
            meanROI(i) = mean(double(vTemp(roiMask)));
            medianROI(i) = median(double(vTemp(roiMask)));
            maxROI(i) = max(double(vTemp(roiMask)));
            minROI(i) = min(double(vTemp(roiMask)));
            stdROI(i) = std(double(vTemp(roiMask)));
            area = sum(roiMask(:))*(xres)^2;
            flowROI(i) = area.*meanROI(i);
        end 
        dFlow = abs(derivative(flowROI));
        dMean = abs(derivative(meanROI));
        
        times = double(timeres.*(1:length(flowROI)));
        timesInterp = double(timeres.*tq);
        meanROIfit = interp1(times,meanROI,timesInterp,'linear');
        medianROIfit = interp1(times,medianROI,timesInterp,'linear');
        maxROIfit = interp1(times,maxROI,timesInterp,'linear');
        minROIfit = interp1(times,minROI,timesInterp,'linear');
        stdROIfit = interp1(times,stdROI,timesInterp,'linear');
        flowROIfit = interp1(times,flowROI,timesInterp,'linear');
        dFlowfit = interp1(times,dFlow,timesInterp,'linear');
        dMeanfit = interp1(times,dMean,timesInterp,'linear');
        

        roiStatistics.radius = radius; roiStatistics.center = center;
        roiStatistics.roiMask = roiMask; roiStatistics.roiDataRaw = roiDataRaw;
        roiStatistics.times = timesInterp;
        roiStatistics.meanROI = meanROIfit; roiStatistics.dMean = dMeanfit;
        roiStatistics.medianROI = medianROIfit;
        roiStatistics.maxROI = maxROIfit; roiStatistics.minROI = minROIfit;
        roiStatistics.stdROI = stdROIfit;
        roiStatistics.flowROI = flowROIfit; roiStatistics.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdata)
            handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
        else 
            handles.pcDatasets(planeNum).ROIdata = roiStatistics;
        end 

        %%% Create Makima Fit
        meanROIfit = interp1(times,meanROI,timesInterp,'makima');
        medianROIfit = interp1(times,medianROI,timesInterp,'makima');
        maxROIfit = interp1(times,maxROI,timesInterp,'makima');
        minROIfit = interp1(times,minROI,timesInterp,'makima');
        stdROIfit = interp1(times,stdROI,timesInterp,'makima');
        flowROIfit = interp1(times,flowROI,timesInterp,'makima');
        dFlowfit = interp1(times,dFlow,timesInterp,'makima');
        dMeanfit = interp1(times,dMean,timesInterp,'makima');

        roiStatisticsMakima.times = timesInterp;
        roiStatisticsMakima.meanROI = meanROIfit; roiStatisticsMakima.dMean = dMeanfit;
        roiStatisticsMakima.medianROI = medianROIfit;
        roiStatisticsMakima.maxROI = maxROIfit; roiStatisticsMakima.minROI = minROIfit;
        roiStatisticsMakima.stdROI = stdROIfit;
        roiStatisticsMakima.flowROI = flowROIfit; roiStatisticsMakima.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdataMakima)
            handles.pcDatasets(planeNum).ROIdataMakima(end+1) = roiStatisticsMakima;
        else 
            handles.pcDatasets(planeNum).ROIdataMakima = roiStatisticsMakima;
        end 

        %%% Create Interpolated Curve with Gaussian Smoothing           
        meanROIfit = interp1(times,smoothdata(meanROI,'gaussian',5),timesInterp,'linear');
        medianROIfit = interp1(times,smoothdata(medianROI,'gaussian',5),timesInterp,'linear');
        maxROIfit = interp1(times,smoothdata(maxROI,'gaussian',5),timesInterp,'linear');
        minROIfit = interp1(times,smoothdata(minROI,'gaussian',5),timesInterp,'linear');
        stdROIfit = interp1(times,smoothdata(stdROI,'gaussian',5),timesInterp,'linear');
        flowROIfit = interp1(times,smoothdata(flowROI,'gaussian',5),timesInterp,'linear');
        dFlowfit = interp1(times,smoothdata(dFlow,'gaussian',5),timesInterp,'linear');
        dMeanfit = interp1(times,smoothdata(dMean,'gaussian',5),timesInterp,'linear');

        roiStatisticsGaussian.radius = radius; roiStatisticsGaussian.center = center;
        roiStatisticsGaussian.roiMask = roiMask; roiStatisticsGaussian.roiDataRaw = roiDataRaw;
        roiStatisticsGaussian.times = timesInterp;
        roiStatisticsGaussian.meanROI = meanROIfit; roiStatisticsGaussian.dMean = dMeanfit;
        roiStatisticsGaussian.medianROI = medianROIfit;
        roiStatisticsGaussian.maxROI = maxROIfit; roiStatisticsGaussian.minROI = minROIfit;
        roiStatisticsGaussian.stdROI = stdROIfit;
        roiStatisticsGaussian.flowROI = flowROIfit; roiStatisticsGaussian.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdataGaussian)
            handles.pcDatasets(planeNum).ROIdataGaussian(end+1) = roiStatisticsGaussian;
        else 
            handles.pcDatasets(planeNum).ROIdataGaussian = roiStatisticsGaussian;
        end 

        %%% Create Cubic Spline Fit   
        meanROIfit = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
        medianROIfit = csaps([0 times times(end)+times(1)],[0 medianROI 0],0.0001,timesInterp);
        maxROIfit = csaps([0 times times(end)+times(1)],[0 maxROI 0],0.0001,timesInterp);
        minROIfit = csaps([0 times times(end)+times(1)],[0 minROI 0],0.0001,timesInterp);
        stdROIfit = csaps([0 times times(end)+times(1)],[0 stdROI 0],0.0001,timesInterp);
        flowROIfit = csaps([0 times times(end)+times(1)],[0 flowROI 0],0.0001,timesInterp);
        dFlowfit = csaps([0 times times(end)+times(1)],[0 dFlow 0],0.0001,timesInterp);
        dMeanfit = csaps([0 times times(end)+times(1)],[0 dMean 0],0.0001,timesInterp);

        roiStatisticsSpline.times = timesInterp;
        roiStatisticsSpline.meanROI = meanROIfit; roiStatisticsSpline.dMean = dMeanfit;
        roiStatisticsSpline.medianROI = medianROIfit;
        roiStatisticsSpline.maxROI = maxROIfit; roiStatisticsSpline.minROI = minROIfit;
        roiStatisticsSpline.stdROI = stdROIfit;
        roiStatisticsSpline.flowROI = flowROIfit; roiStatisticsSpline.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdataSpline)
            handles.pcDatasets(planeNum).ROIdataSpline(end+1) = roiStatisticsSpline;
        else 
            handles.pcDatasets(planeNum).ROIdataSpline = roiStatisticsSpline;
        end 
        if ~exist('ROIimages','dir')
            mkdir ROIimages
            cd ROIimages
            frame = getframe(handles.PlanePlot);
            image = frame2im(frame);
            imwrite(image,[handles.pcDatasets(planeNum).Names '.png'])
            cd ..
        else
            cd ROIimages
            frame = getframe(handles.PlanePlot);
            image = frame2im(frame);
            imwrite(image,[handles.pcDatasets(planeNum).Names '.png'])
            cd ..
        end 
        
        clear dFlow dMean radius center roiMask area planeNum timeres v
        clear roiDataRaw meanROI medianROI maxROI minROI stdROI flowROI  
        clear dFlowfit dMeanfit circle matrixx t tq vTemp X Y xres i fovx
        clear meanROIfit medianROIfit maxROIfit minROIfit stdROIfit flowROIfit  
        clear roiStatistics roiStatisticsMakima roiStatisticsGaussian roiStatisticsSpline
        clear times timesInterp
    end 

    %%% Plot Curves
    handles.flow = organizeFlowInfo(handles);
    plotVelocity(handles);

    set(handles.InterpolatePopup,'Enable','on');
    guidata(hObject,handles);
    DeleteROIbutton_Callback(hObject, eventdata, handles);
    updateImages(handles);
    axes(handles.PlanePlot);
    
% --- COMPLETE LOADING BUTTON - CALLBACK
function CompleteLoadingROI_Callback(hObject, eventdata, handles)
    if isfield(handles,'flow') && numel(handles.flow)>1
        set(handles.DeleteROIbutton,'Enable','off');
        set(handles.LoadROIbutton,'Enable','off');

        set(handles.DrawCenterlineButton,'Enable','on');
        set(handles.ShowPlanesRadio,'Enable','on');
        set(handles.ZoomAnatomicalButton,'Enable','on');
        set(handles.UnzoomAnatomicalButton,'Enable','on');

        zLocs = EnterZslices(handles.flow);

        for i=1:numel(handles.flow)
            handles.flow(i).zLocs = zLocs(i);
        end 

        for i =1:numel(handles.anatDatasets)
            images = handles.anatDatasets(i).Data;
            dims(1) = size(images,1);
            dims(2) = size(images,2);
            dims(3) = size(images,3);
            handles.anatDatasets(i).dims = dims;

            originShift = [handles.anatDatasets(i).Info.ImagePositionPatient;1]; % origin is top left corner of image
            xres = handles.anatDatasets(i).Info.PixelSpacing(1);
            yres = handles.anatDatasets(i).Info.PixelSpacing(2);
            zres = handles.anatDatasets(i).Info.SliceThickness;

            % sometimes get extremely small values that should be 0, so round
            xVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(1:3) ,8); % what direction rows run w/r/to x
            yVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(4:6) ,8); % what direction the cols run w/r/to y
            zVector = [cross(xVector,yVector);0];
            xVector = [xVector;0];
            yVector = [yVector;0];
            rotationMatrix = [xres*xVector yres*yVector zres*zVector originShift];
            sliceDim = find(rotationMatrix(:,3));
            handles.anatDatasets(i).rotationMatrix = rotationMatrix;     
            handles.anatDatasets(i).sliceDim = sliceDim;


            topFrontLeft = rotationMatrix*[1;1;1;1];
            topFrontLeft(4) = []; % remove dummy dimension
            backBottomRight = rotationMatrix*[dims(1);dims(2);dims(3);1];
            backBottomRight(4) = [];
            spanZ(1) = topFrontLeft(3);
            spanZ(2) = backBottomRight(3);
            handles.anatDatasets(i).spanZ = spanZ;
        end   
        
        if ~exist('FlowWaveformImages','dir')
            mkdir FlowWaveformImages
            cd FlowWaveformImages
            frame = getframe(handles.VelocityPlot);
            imwrite(frame2im(frame),'FlowWaveform.png')
            cd ..
        else
            cd FlowWaveformImages
            frame = getframe(handles.VelocityPlot);
            imwrite(frame2im(frame),'FlowWaveform.png')
            cd ..
        end 

        clear backBottomRight dims i images originShift rotationMatrix sliceDim
        clear spanZ zLocs topFrontLeft xres xVector yres yVector zres zVector
        guidata(hObject,handles);
    else
        set(handles.MessageBar,'String','You need at least 2 ROIs for PWV measurements.');
    end 
    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function VelocityPlot_CreateFcn(hObject, eventdata, handles)
    

% --- INTERPOLATE POPUP - CALLBACK
function InterpolatePopup_Callback(hObject, eventdata, handles)
    interp = get(handles.InterpolatePopup,'Value');
    switch interp
        case 1
            handles.global.interpType = 'None';
        case 2 
            handles.global.interpType = 'Makima';
        case 3 
            handles.global.interpType = 'Gaussian';
        case 4
            handles.global.interpType = 'Spline';
        otherwise
    end 
 
    clear interp
    guidata(hObject, handles);
    plotVelocity(handles)
    
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- INTERPOLATE POPUP - CREATE FUNCTION
function InterpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ERROR BAR RADIO - CALLBACK
function ErrorBarRadio_Callback(hObject, eventdata, handles)
    if get(handles.ErrorBarRadio,'Value')
        handles.global.showErrorBars = 1;
    else 
        handles.global.showErrorBars = 0;
    end 
    
    guidata(hObject,handles);
    plotVelocity(handles)
    
    
% --- PG SHIFT RADIO - CALLBACK
function PGshiftRadio_Callback(hObject, eventdata, handles)
    if get(handles.PGshiftRadio,'Value')
        handles.global.pgShift = 1;
    else 
        handles.global.pgShift = 0;
    end 
    
    guidata(hObject,handles);
    plotVelocity(handles)
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 


%%%%%%%%%%%% ANATOMICAL PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- ANATOMICAL PLOT - CREATE FUNCTION
function AnatomicalPlot_CreateFcn(hObject, eventdata, handles)


% --- ANATOMICAL SLIDER - CALLBACK
function SliderAnatomical_Callback(hObject, eventdata, handles)
    updateAnatImages(handles);

% --- ANATOMICAL SLIDER - CREATE FUNCTION
function SliderAnatomical_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

 
% --- ANATOMICAL LISTBOX - CALLBACK
function AnatListbox_Callback(hObject, eventdata, handles)
    ShowPlanesRadio_Callback(hObject, eventdata, handles)
    updateAnatImages(handles);

% --- ANATOMICAL LISTBOX - CREATE FUNCTION
function AnatListbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SHOW PLANE RADIO - CALLBACK
function ShowPlanesRadio_Callback(hObject, eventdata, handles) 
    if get(handles.ShowPlanesRadio,'Value')
        datasetNum = get(handles.AnatListbox,'Value');
        yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
        anatRotation = handles.anatDatasets(datasetNum).rotationMatrix;
        colsRunningDir = sign(nonzeros(anatRotation(:,3)));
        spanZ = handles.anatDatasets(datasetNum).spanZ;
        for i=1:numel(handles.flow)
            planeLinePhysical = colsRunningDir*(handles.flow(i).zLocs-spanZ(1));
            planeLineRow = round(planeLinePhysical/yres);
            handles.flow(i).planeLineRow = planeLineRow;
            guidata(hObject,handles);
        end 
    end 
    
    clear i datasetNum yres anatRotation colsRunningDir spanZ 
    clear planeLinePhysical planeLineRow
    updateAnatImages(handles);
    guidata(hObject,handles);
    
% --- ZOOM ANATOMICAL - CALLBACK
function ZoomAnatomicalButton_Callback(hObject, eventdata, handles)
    axes(handles.AnatomicalPlot)
    
    set(handles.MessageBar,'String','Draw an ROI to zoom in');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUDAnat = handles.global.spanUDAnat;
    spanLRAnat = handles.global.spanLRAnat;
    spanUDAnat = spanUDAnat(1)+positions(2):(spanUDAnat(1)+positions(2)+positions(4));
    spanLRAnat = spanLRAnat(1)+positions(1):(spanLRAnat(1)+positions(1)+positions(3));
    handles.global.zoomedAnat = 1;
    
    clear rect positions 
    guidata(hObject, handles);
    updateAnatImages(handles);
    
    
% --- UNZOOM ANATOMICAL - CALLBACK
function UnzoomAnatomicalButton_Callback(hObject, eventdata, handles)
    axes(handles.AnatomicalPlot)
    
    handles.global.zoomedAnat = 0; 
    handles.global.spanUDAnat = 0; 
    handles.global.spanLRAnat = 0;
    
    guidata(hObject, handles);
    updateAnatImages(handles);


% --- DRAW CENTERLINE - CALLBACK
function DrawCenterlineButton_Callback(hObject, eventdata, handles)
%         axes(handles.AnatomicalPlot); 
%         if isempty(get(handles.range1Edit,'String')) || isempty(get(handles.range1Edit,'String'))
%             set(handles.MessageBar,'String','Range is required to draw centerlines')
%         else 
        
            
            
        ShowPlanesRadio_Callback(hObject, eventdata, handles);
        
        set(handles.ShowPlanesRadio,'Value',1);
        set(handles.ShowPlanesRadio,'Enable','off');
        set(handles.DrawROIbutton,'Enable','off');
        set(handles.DeleteCenterlineButton,'Enable','on');
        set(handles.ComputePWVButton,'Enable','on');
        set(handles.SliderAnatomical,'Enable','off');
        set(handles.AnatListbox,'Enable','off');
        set(handles.ZoomAnatomicalButton,'Enable','off');
        set(handles.UnzoomAnatomicalButton,'Enable','off');
        set(handles.ShowPlanesRadio,'Value',1);        
        
        range1 = str2double(get(handles.range1Edit,'String'));
        range2 = str2double(get(handles.range2Edit,'String'));
        datasetNum = get(handles.AnatListbox,'Value');
        
        xres = handles.anatDatasets(datasetNum).Info.PixelSpacing(1);
        yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
        sliceres = handles.anatDatasets(datasetNum).Info.SliceThickness;
        
        anatRotation = handles.anatDatasets(datasetNum).rotationMatrix;
        colsRunningDir = sign(nonzeros(anatRotation(:,3)));
        spanZ = handles.anatDatasets(datasetNum).spanZ;
        for i=1:numel(handles.flow)
            planeLinePhysical = colsRunningDir*(handles.flow(i).zLocs-spanZ(1));
            planeLineRow = round(planeLinePhysical/yres);
            handles.flow(i).planeLineRow = planeLineRow;
            guidata(hObject,handles);
        end 
        
        
        
        
        
%         axes(handles.AnatomicalPlot);
%         images = handles.anatDatasets(datasetNum).Data(:,:,range1:range2);
%         x = []; y = []; z = [];
%         for i=1:size(images,3)
%             image = images(:,:,i);
%             imshow(image,[]);
%             [xTemp,yTemp] = getpts();
%             zTemp = i.*(ones(size(xTemp,1),1));
%             x = [x; xTemp];
%             y = [y; yTemp];
%             z = [z; zTemp];
%         end 
%         
%         figure; scatter(y,x); xline(127); xline(325);
%         line1 = drawpolyline;
%         figure; scatter(y,z);  xline(127); xline(325);
%         line2 = drawpolyline;
%                 
%         splinePositions1 = interppolygon(line1.Position,150);
%         splinePositions2 = interppolygon(line2.Position,150);
%         %NEED TO EDIT IF CORONAL
%         splineLine(:,1) = flipud(splinePositions2(:,2)); %x physical
%         splineLine(:,2) = splinePositions1(:,2); %y
%         splineLine(:,3) = splinePositions1(:,1); %z
%         figure; scatter3(y,x,z); hold on; 
%         scatter3(splineLine(:,3),splineLine(:,2),splineLine(:,1),'filled');
%         xlabel('x'); ylabel('y'); zlabel('slice');
%         x = linspace( min(splineLine(:,2))-50, max(splineLine(:,2))+50, 512);
%         y = linspace( min(splineLine(:,1))-2, max(splineLine(:,1))+2, 512);
%         [X,Y] = meshgrid(x,y);
%         for i=1:numel(handles.flow)
%             Z = (handles.flow(i).planeLineRow).*ones(size(X));
%             hold on; surf(Z,X,Y);
%         end 
%         hold off;

%         figure; scatter3(splineLine(:,3),splineLine(:,2),splineLine(:,1),30,'filled','g'); hold on;
%         xlabel('x'); ylabel('y'); zlabel('slice');
%         image1 = handles.pcDatasets(1).Data.MAG;
%         image1 = image1(100:412,100:412);
%         image2 = handles.pcDatasets(2).Data.MAG;
%         image2 = image2(100:412,100:412);
%         
%         x1 = (linspace( min(splineLine(:,2))-80, max(splineLine(:,2))+80, size(image1,1))*1.4)-110;
%         x2 = linspace( min(splineLine(:,2))-80, max(splineLine(:,2))+80, size(image2,1))+15;
%         
%         y1 = linspace( min(splineLine(:,1))-8, max(splineLine(:,1))+8, size(image1,1))-0.7;
%         y1 = repmat(y1,size(image1,1),1);
%         
%         y2 = linspace( min(splineLine(:,1))-8, max(splineLine(:,1))+8, size(image2,1));
%         y2 = repmat(y2,size(image2,1),1)-0.25;
%         
%         z1 = (handles.flow(1).planeLineRow).*ones(1,length(x1));
%         z2 = (handles.flow(3).planeLineRow).*ones(1,length(x2));
%         
%         surf(z1,x1,y1,image1); colormap('gray'); hold on;
%         surf(z2,x2,y2,image2); colormap('gray'); shading interp
% 
%         
%         distances = zeros(1,length(splinePositions1)-1);
%         for i=1:length(splinePositions1)-1
%             xSquare = (sliceres.*(splineLine(i,1)-splineLine(i+1,1))).^2;
%             ySquare = (xres.*(splineLine(i,2)-splineLine(i+1,2))).^2;
%             zSquare = (yres.*(splineLine(i,3)-splineLine(i+1,3))).^2;
%             distances(i) = sqrt( xSquare + ySquare + zSquare );
%         end 
%         distances(end+1)=0;
%         handles.anatDatasets(datasetNum).Distances = distances;
%         handles.anatDatasets(datasetNum).Centerline = splineLine;
% 
%         clear datasetNum yres anatRotation colsRunningDir spanZ i w mydlg 
%         clear key planeLinePhysical planeLineRow matrixx fovx xres matrixy 
%         clear fovy yres distances splinePositions
%         
%         guidata(hObject,handles);
%         updateImages(handles);

% % %         if isempty(handles.anatDatasets(datasetNum).Centerline)
% % %             mydlg = warndlg('Press enter when the centerline is drawn');
% % %             waitfor(mydlg);
% % %             line = drawpolyline('Color','r','LineWidth',1);
% % %             while true
% % %                 w = waitforbuttonpress; 
% % %                 switch w 
% % %                     case 1 % (keyboard press) 
% % %                       key = get(gcf,'currentcharacter'); 
% % %                           switch key
% % %                               case 27 % 27 is the escape key
% % %                                   set(handles.MessageBar,'String','User pressed the escape key. Deleting ROI.');
% % %                                   DeleteCenterlineButton_Callback(hObject, eventdata, handles)
% % %                                   break % break out of the while loop
% % %                               case 13 % 13 is the return key 
% % %                                   set(handles.MessageBar,'String','ROI selected.');
% % %                                   line.InteractionsAllowed = 'none';
% % %                                   break
% % %                               otherwise 
% % %                                   % Wait for a different command. 
% % %                           end
% % %                end
% % %             end
% % % 
% % %         splinePositions = interppolygon(line.Position,100);
% % %         line.Position = splinePositions;
% % % 
% % %         if isfield(handles.anatDatasets(datasetNum).Info,'matrixx')
% % %             matrixx = handles.anatDatasets(datasetNum).Info.matrixx;
% % %             fovx = handles.anatDatasets(datasetNum).Info.fovx;
% % %             xres = fovx/matrixx;
% % %             matrixy = handles.anatDatasets(datasetNum).Info.matrixy;
% % %             fovy = handles.anatDatasets(datasetNum).Info.fovy;
% % %             yres = fovy/matrixy;
% % %         else 
% % %             xres = handles.anatDatasets(datasetNum).Info.PixelSpacing(1);
% % %             yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
% % %         end
% % % 
% % %         distances = zeros(1,length(splinePositions)-1);
% % %         for i=1:length(splinePositions)-1
% % %             distances(i) = sqrt( (xres.*(splinePositions(i,1)-splinePositions(i+1,1))).^2 + (yres.*(splinePositions(i,2)-splinePositions(i+1,2))).^2 );
% % %         end 
% % %         distances(end+1)=0;
% % %         handles.anatDatasets(datasetNum).Distances = distances;
% % %         handles.anatDatasets(datasetNum).Centerline = line;
% % % 
% % %         clear datasetNum yres anatRotation colsRunningDir spanZ i w mydlg 
% % %         clear key planeLinePhysical planeLineRow matrixx fovx xres matrixy 
% % %         clear fovy yres distances splinePositions
% % %         guidata(hObject,handles);
% % %         updateImages(handles);
% % %         else
% % %             fprintf('A Centerline has already been placed!\n');
% % %         end 
%         end
    
    datasetNum = get(handles.AnatListbox,'Value');
    oldData = load('D:\PWV\volunteers\Mulan-Volunteer_PWV_01833_2019-11-02-h11\01833_00008_pwv-radial_Aao_PG\2500proj\Data-DataAnalysis\anatDatasets.mat');
    handles.anatDatasets(datasetNum).Distances = oldData.anatDatasets(datasetNum).Distances;
    handles.anatDatasets(datasetNum).Centerline = oldData.anatDatasets(datasetNum).Centerline;
    guidata(hObject,handles);
    clear oldData  
    
% --- DELETE CENTERLINE - CALLBACK
function DeleteCenterlineButton_Callback(hObject, eventdata, handles)
    cla(handles.AnatomicalPlot,'reset');
    
    datasetNum = get(handles.AnatListbox,'Value');
    handles.anatDatasets(datasetNum).Centerline = [];
    handles.anatDatasets(datasetNum).Length = [];

    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.DrawROIbutton,'Enable','on');
    set(handles.AnatListbox,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    set(handles.SliderAnatomical,'Enable','on');
    set(handles.ShowPlanesRadio,'Value',0);
    
    
    clear datasetNum
    hold off
    guidata(hObject,handles);
    updateAnatImages(handles);
    
    
   
% --- SLICE NUMBER TEXT - CREATE FUNCTION
function sliceNumText_CreateFcn(hObject, eventdata, handles)
    

% --- RANGE 1 EDIT - CALLBACK
function range1Edit_Callback(hObject, eventdata, handles)

% --- RANGE 1 EDIT - CREATE FUNCTION
function range1Edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- RANGE 2 EDIT - CALLBACK
function range2Edit_Callback(hObject, eventdata, handles)

% --- RANGE 2 EDIT - CREATE FUNCTION
function range2Edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
   
    
% --- COMPUTE PWV - CALLBACK
function ComputePWVButton_Callback(hObject, eventdata, handles)
    handles.global.startAnalyzing = 1;
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.AnatListbox,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    set(handles.SliderAnatomical,'Enable','on');

    distances = [handles.anatDatasets.Distances];
    centerline = [handles.anatDatasets.Centerline];
    planeRows = [handles.flow.planeLineRow];
    y = centerline(:,3);
    
    difference = ones(length(y),1);
    for i=1:length(planeRows)
        difference = difference.*(y-planeRows(i));
    end 
    minima = islocalmin(abs(difference));
    
    if y(1)<planeRows(1)
        minima(1) = 1;
    end 
    
    if y(end)<planeRows(end)<0
        minima(end) = 1;
    end 
    roiIdxAlongLine = nonzeros(minima.*(1:length(y))');
    
    numROIs = numel(handles.flow);  
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        distances2ROIs = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                idx1 = roiIdxAlongLine(i);
                idx2 = roiIdxAlongLine(iterator);
                dist = sum(distances(idx1:idx2));
                distances2ROIs(iterator) = dist;
            end 
        end 
        handles.flow(i).Distance = distances2ROIs;
    end 
    
    guidata(hObject,handles);
    flow = computeTTs(handles.flow);
    handles.flow = flow;
    guidata(hObject,handles);
    
    
    distance = []; TTPeak = []; TTPoint = []; TTFoot = []; TTUpstroke = []; Xcorr = [];
    for i = 1:numel(handles.flow)
        mask = ~isnan(handles.flow(i).Distance);
        distance = [distance handles.flow(i).Distance(mask)];
        TTPeak = [TTPeak handles.flow(i).TTPeak(mask)];
        TTPoint = [TTPoint handles.flow(i).TTPoint(mask)];
        TTFoot = [TTFoot handles.flow(i).TTFoot(mask)];
        TTUpstroke = [TTUpstroke handles.flow(i).TTUpstroke(mask)];
        Xcorr = [Xcorr handles.flow(i).Xcorr(mask)];
    end
     
    cla(handles.TimeVsDistance,'reset');
    axes(handles.TimeVsDistance); hold on; 
    xlabel('Distance (mm)'); ylabel('Time Shift (ms)');
    sz = 30;
    
    numMethods = get(handles.ttpRadio,'Value')+get(handles.ttpointRadio,'Value')...
    +get(handles.ttfRadio,'Value')+get(handles.ttuRadio,'Value')+get(handles.xcorrRadio,'Value');
    numCompares = numel(distance);
    
    legendSet = {};
    average = zeros(1,numCompares);
    if get(handles.ttpRadio,'Value')
        scatter(distance,TTPeak,sz,'filled','MarkerFaceColor',[0 0.4470 0.7410]);
        for i=1:numCompares
            average(i) = average(i)+TTPeak(i);
        end 
        legendSet{end+1} = 'TTPeak';
    end 
    if get(handles.ttpointRadio,'Value')
        scatter(distance,TTPoint,sz,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
        for i=1:numCompares
            average(i) = average(i)+TTPoint(i);
        end 
        legendSet{end+1} = 'TTPoint';
    end 
    if get(handles.ttfRadio,'Value')
        scatter(distance,TTFoot,sz,'filled','MarkerFaceColor',[0.4940 0.1840 0.5560]);
        for i=1:numCompares
            average(i) = average(i)+TTFoot(i);
        end 
        legendSet{end+1} = 'TTFoot';
    end 
    if get(handles.ttuRadio,'Value')
        scatter(distance,TTUpstroke,sz,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
        for i=1:numCompares
            average(i) = average(i)+TTUpstroke(i);
        end 
        legendSet{end+1} = 'TTUpstroke';
    end 
    if get(handles.xcorrRadio,'Value')
        scatter(distance,Xcorr,sz,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]);
        for i=1:numCompares
            average(i) = average(i)+Xcorr(i);
        end 
        legendSet{end+1} = 'XCorr';
    end 

    average = average/numMethods;
    scatter(distance,average,40,'black');
    legendSet{end+1} = 'AVERAGE';
    
    d = min(distance):max(distance);
    if numel(distance)>2
        linePeakFit = polyfit(distance,TTPeak,1);
        linePeak = linePeakFit(1)*d + linePeakFit(2);
        linePointFit = polyfit(distance,TTPoint,1);
        linePoint= linePointFit(1)*d + linePointFit(2);
        lineFootFit = polyfit(distance,TTFoot,1);
        lineFoot= lineFootFit(1)*d + lineFootFit(2);
        lineUpstrokeFit = polyfit(distance,TTUpstroke,1);
        lineUpstroke = lineUpstrokeFit(1)*d + lineUpstrokeFit(2);
        lineXcorrFit = polyfit(distance,Xcorr,1);
        lineXcorr = lineXcorrFit(1)*d + lineXcorrFit(2);
        lineAverageFit = polyfit(distance,average,1);
        lineAverage= lineAverageFit(1)*d + lineAverageFit(2);
        PWVpeak = 1/linePeakFit(1);
        PWVpoint = 1/linePointFit(1);
        PWVfoot = 1/lineFootFit(1);
        PWVupstroke = 1/lineUpstrokeFit(1);
        PWVxcorr = 1/lineXcorrFit(1);
        PWVaverage = 1/lineAverageFit(1);
        
        hold on; 
        if get(handles.ttpRadio,'Value')
            plot(d,linePeak,':','LineWidth',0.2,'MarkerFaceColor',[0 0.4470 0.7410]);
            set(handles.ttpData,'String',[num2str(round(PWVpeak,2)) ' m/s']);
        end 
        if get(handles.ttpointRadio,'Value')
            plot(d,linePoint,':','LineWidth',0.2,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
            set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']);
        end 
        if get(handles.ttfRadio,'Value')
            plot(d,lineFoot,':','LineWidth',0.2,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
            set(handles.ttfData,'String',[num2str(round(PWVfoot,2)) ' m/s']);
        end 
        if get(handles.ttuRadio,'Value')
            plot(d,lineUpstroke,':','LineWidth',0.2,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
            set(handles.ttuData,'String',[num2str(round(PWVupstroke,2)) ' m/s']);
        end 
        if get(handles.xcorrRadio,'Value')
            plot(d,lineXcorr,':','LineWidth',0.2,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
            set(handles.xcorrData,'String',[num2str(round(PWVxcorr,2)) ' m/s']);
        end 
        plot(d,lineAverage,'-k','LineWidth',0.2);
        legend(legendSet,'Location','southeast');
        hold off;
    
    else 
        PWVpeak = distance/TTPeak;
        PWVpoint = distance/TTPoint;
        PWVfoot = distance/TTFoot;
        PWVupstroke = distance/TTUpstroke;
        PWVxcorr = distance/Xcorr;
        PWVaverage = distance/average;
        
        if get(handles.ttpRadio,'Value')
            set(handles.ttpData,'String',[num2str(round(PWVpeak,2)) ' m/s']);
        end 
        if get(handles.ttpointRadio,'Value')
            set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']);
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
    
    
    

    
    if numel(distance)==1
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
    elseif numel(distance)==3
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' Plane 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' Plane 2 --> 3');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']);       
    else
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' Plane 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' Plane 1 --> 4');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']); 
        set(handles.distanceHeader4,'String',' Plane 2 --> 3');
        set(handles.distance4,'String',[num2str(round(distance(4),1)) ' mm']);
        set(handles.distanceHeader5,'String',' Plane 2 --> 4');
        set(handles.distance5,'String',[num2str(round(distance(5),1)) ' mm']); 
        set(handles.distanceHeader6,'String',' Plane 3 --> 4');
        set(handles.distance6,'String',[num2str(round(distance(6),1)) ' mm']);         
    end 
    
    
    if numMethods>0
        set(handles.averageData,'String',[num2str(round(PWVaverage,2)) ' m/s']);
    else 
        set(handles.averageData,'String','0 m/s');
    end
    
    if ~exist('CenterlineImages','dir')
        mkdir CenterlineImages
        cd CenterlineImages
        frame = getframe(handles.AnatomicalPlot);
        imwrite(frame2im(frame),'Centerline.png')
        cd ..
    else
        cd CenterlineImages
        frame = getframe(handles.AnatomicalPlot);
        imwrite(frame2im(frame),'Centerline.png')
        cd ..
    end 
    
    if ~exist('PWVanalysisPlotImages','dir')
        mkdir PWVanalysisPlotImages
        cd PWVanalysisPlotImages
        frame = getframe(handles.TimeVsDistance);
        imwrite(frame2im(frame),'PWVanalysisPlot.png')
        cd ..
    else
        cd PWVanalysisPlotImages
        frame = getframe(handles.TimeVsDistance);
        imwrite(frame2im(frame),'PWVanalysisPlot.png')
        cd ..
    end 
    
    clear allROIs average centerline d difference dist distance distances distances2ROIs
    clear i idx idx2 iterator j legendSet lineAverage lineAverageFit idx1
    clear lineFoot lineFootFit linePeak linePeakFit linePoint linePointFit
    clear lineUpstroke lineUpstrokeFit lineXcorr lineXcorrFit mask
    clear minima numCompares numMethods numROIs planeRows PWVaverage PWVfoot 
    clear PWVpeak PWVpoint PWVupstroke PWVxcorr sz startAnalyzing Xcorr
    clear TTFoot TTPeak TTPoint TTs TTUpstroke y ROIs2compare roiIdxAlongLine
    
% --- COMPLETE RESET BUTTON - CALLBACK
function resetButton_Callback(hObject, eventdata, handles)
    answer = questdlg('Pressing YES will reset all ROI data. Pressing NO will cancel. Would you like to reset data?', ...
        'Complete Reset', ...
        'YES','NO','NO');
    % Handle response
    switch answer
        case 'YES'
            reset = 1;
        case 'NO'
            reset = 0;
    end
    
    if reset
        set(handles.ttpData,'String',' ');
        set(handles.ttpointData,'String',' ');
        set(handles.ttfData,'String',' ');
        set(handles.ttuData,'String',' ');
        set(handles.xcorrData,'String',' ');
        set(handles.xcorrData,'String',' ');
        set(handles.averageData,'String',' ');

        set(handles.distanceHeader1,'String',' ');
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
            handles = rmfield(handles,'flow');
        end 
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROI');
        handles.magDatasets = rmfield(handles.magDatasets,'ROI');
        
        if isfield(handles.anatDatasets,'Centerline')
            handles.anatDatasets = rmfield(handles.anatDatasets,'Centerline');
        end 
        if isfield(handles.anatDatasets,'Centerline')
            handles.anatDatasets = rmfield(handles.anatDatasets,'Distances');
        end 
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdata');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdataMakima');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdataGaussian');
        handles.pcDatasets = rmfield(handles.pcDatasets,'ROIdataSpline');
        
        handles.pcDatasets(1).ROI = [];
        handles.magDatasets(1).ROI = [];
        handles.anatDatasets(1).Centerline = [];
        handles.anatDatasets(1).Distances = [];
        handles.pcDatasets(1).ROIdata = [];
        handles.pcDatasets(1).ROIdataMakima = [];
        handles.pcDatasets(1).ROIdataGaussian = [];
        handles.pcDatasets(1).ROIdataSpline = [];
        
        handles.global.zoomed = 0;
        handles.global.zoomedAnat = 0;
        handles.global.spanUD = [0,0];
        handles.global.spanLR = [0,0];
        handles.global.spanUDAnat = [0,0];
        handles.global.spanLRAnat = [0,0];
        handles.global.interpType = 'None';
        handles.global.showErrorBars = 0;
        handles.global.startAnalyzing = 0;
        handles.global.pgShift = 0;

        set(handles.DrawROIbutton,'Enable','on');
        set(handles.DeleteROIbutton,'Enable','on');
        set(handles.LoadROIbutton,'Enable','on');
        set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
        set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
        set(handles.InterpolatePopup,'String',{'None','Makima','Gaussian','Spline'});
        set(handles.AnatListbox,'String',{handles.anatDatasets.Names});
        set(handles.InterpolatePopup,'Enable','off');
        if numel(handles.magDatasets)==0
            set(handles.MAGRadio,'Enable','off')
        end 
        set(handles.PlanePopup,'Value',1);
        set(handles.DatasetPopup,'Value',1);

        set(handles.ShowPlanesRadio,'Enable','off');
        set(handles.DrawCenterlineButton,'Enable','off');
        set(handles.DeleteCenterlineButton,'Enable','off');
        set(handles.DeleteCenterlineButton,'Enable','off');
        set(handles.ComputePWVButton,'Enable','off');
        set(handles.ShowPlanesRadio,'Value',0)
        set(handles.ShowPlanesRadio,'Enable','off')
        set(handles.ttpRadio,'Value',1);
        set(handles.ttpointRadio,'Value',1);
        set(handles.ttuRadio,'Value',1);
        set(handles.ttfRadio,'Value',1);
        set(handles.xcorrRadio,'Value',1);
        set(handles.InterpolatePopup,'Value',1);
        set(handles.PGshiftRadio,'Value',0);
        
        clear mydlg w key 
        cla(handles.VelocityPlot,'reset')
        cla(handles.TimeVsDistance,'reset')
        cla(handles.AnatomicalPlot,'reset')
        cla(handles.PlanePlot,'reset')
        guidata(hObject, handles);
        updateImages(handles);
        updateAnatImages(handles);
    end 
    




   
%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- EXPORT ANALYSIS - CALLBACK
function ExportAnalysisButton_Callback(hObject, eventdata, handles)
    interpTypes = {'None', 'Makima', 'Gaussian', 'Spline'};

    for t=1:length(interpTypes)
        handles.global.interpType = interpTypes{t};
        flow = computeTTs(handles.flow);
        Distance = []; TTPeak = []; TTPoint = []; TTFoot = []; TTUpstroke = []; Xcorr = [];
        
        for i = 1:numel(flow)
            mask = ~isnan(flow(i).Distance);
            Distance = [Distance flow(i).Distance(mask)];
            TTPeak = [TTPeak flow(i).TTPeak(mask)];
            TTPoint = [TTPoint flow(i).TTPoint(mask)];
            TTFoot = [TTFoot flow(i).TTFoot(mask)];
            TTUpstroke = [TTUpstroke flow(i).TTUpstroke(mask)];
            Xcorr = [Xcorr flow(i).Xcorr(mask)];
        end
        

        if numel(Distance)==1
            PLANES{1} = 'Plane 1 --> 2';
        elseif numel(Distance)==3
            PLANES{1} = 'Plane 1 --> 2';
            PLANES{2} = 'Plane 1 --> 3';
            PLANES{3} = 'Plane 2 --> 3';    
        else
            PLANES{1} = 'Plane 1 --> 2';
            PLANES{2} = 'Plane 1 --> 3';
            PLANES{3} = 'Plane 1 --> 4';
            PLANES{4} = 'Plane 2 --> 3';
            PLANES{5} = 'Plane 2 --> 4';
            PLANES{6} = 'Plane 3 --> 4';        
        end 
        
        numMethods = get(handles.ttpRadio,'Value')+get(handles.ttpointRadio,'Value')...
        +get(handles.ttfRadio,'Value')+get(handles.ttuRadio,'Value')+get(handles.xcorrRadio,'Value');
        numCompares = numel(Distance);
        entryCount = numel(Distance);
        
        PWV_Peak = NaN(entryCount+3,1);
        PWV_Point = NaN(entryCount+3,1);
        PWV_Foot = NaN(entryCount+3,1);
        PWV_Upstroke = NaN(entryCount+3,1);
        PWV_Xcorr = NaN(entryCount+3,1);
        PWV_Average = zeros(1,numCompares);
        if get(handles.ttpRadio,'Value')
            for i=1:numCompares
                PWV_Peak(i) = Distance(i)/TTPeak(i);
                PWV_Average(i) = PWV_Average(i)+PWV_Peak(i);
            end 
        end 
        if get(handles.ttpointRadio,'Value')
            for i=1:numCompares
                PWV_Point(i) = Distance(i)/TTPoint(i);
                PWV_Average(i) = PWV_Average(i)+PWV_Point(i);
            end 
        end 
        if get(handles.ttfRadio,'Value')
            for i=1:numCompares
                PWV_Foot(i) = Distance(i)/TTPeak(i);
                PWV_Average(i) = PWV_Average(i)+PWV_Foot(i);
            end 
        end 
        if get(handles.ttuRadio,'Value')
            for i=1:numCompares
                PWV_Upstroke(i) = Distance(i)/TTUpstroke(i);
                PWV_Average(i) = PWV_Average(i)+PWV_Upstroke(i);
            end 
        end 
        if get(handles.xcorrRadio,'Value')
            for i=1:numCompares
                PWV_Xcorr(i) = Distance(i)/Xcorr(i);
                PWV_Average(i) = PWV_Average(i)+PWV_Xcorr(i);
            end 
        end
       
        PWV_Average = PWV_Average/numMethods;
        
        if numel(Distance)>2    
            if get(handles.ttpRadio,'Value')
                [linePeakFit,Speak] = polyfit(Distance,TTPeak,1);
                PWV_Peak(entryCount+1) = 1/linePeakFit(1);
                PWV_Peak(entryCount+2) = linePeakFit(1);
                PWV_Peak(entryCount+3) = linePeakFit(2);
                Speak
            end 
            
            if get(handles.ttpRadio,'Value')
                [linePointFit,Spoint] = polyfit(Distance,TTPoint,1);
                PWV_Point(entryCount+1) = 1/linePointFit(1);
                PWV_Point(entryCount+2) = linePointFit(1);
                PWV_Point(entryCount+3) = linePointFit(2);
                Spoint
            end 
            
            if get(handles.ttfRadio,'Value')
                [lineFootFit,Sfoot] = polyfit(Distance,TTFoot,1);
                PWV_Foot(entryCount+1) = 1/lineFootFit(1);
                PWV_Foot(entryCount+2) = lineFootFit(1);
                PWV_Foot(entryCount+3) = lineFootFit(2);
                Sfoot
            end 
            
            if get(handles.ttuRadio,'Value')
                [lineUpstrokeFit,Supstroke] = polyfit(Distance,TTUpstroke,1);
                PWV_Upstroke(entryCount+1) = 1/lineUpstrokeFit(1);
                PWV_Upstroke(entryCount+2) = lineUpstrokeFit(1);
                PWV_Upstroke(entryCount+3) = lineUpstrokeFit(2);
                Supstroke
            end 
            
            if get(handles.xcorrRadio,'Value')
                [lineXcorrFit,Sxcorr] = polyfit(Distance,Xcorr,1);
                PWV_Xcorr(entryCount+1) = 1/lineXcorrFit(1);
                PWV_Xcorr(entryCount+2) = lineXcorrFit(1);
                PWV_Xcorr(entryCount+3) = lineXcorrFit(2);
                Sxcorr
            end 
            
            [lineAverageFit,Saverage] = polyfit(Distance,PWV_Average,1);
            PWV_Average(entryCount+1) = 1/lineAverageFit(1);
            PWV_Average(entryCount+2) = lineAverageFit(1);
            PWV_Average(entryCount+3) = lineAverageFit(2);
            PWV_Average = PWV_Average';
            Saverage
        else 
            if get(handles.ttpRadio,'Value')
                PWV_Peak(entryCount+1) = Distance/TTPeak;
                PWV_Peak(entryCount+2) = NaN;
                PWV_Peak(entryCount+3) = NaN;
            end 
            
            if get(handles.ttpointRadio,'Value')
                PWV_Point(entryCount+1) = Distance/TTPoint;
                PWV_Point(entryCount+2) = NaN;
                PWV_Point(entryCount+3) = NaN;
            end 
            
            if get(handles.ttfRadio,'Value')
                PWV_Foot(entryCount+1) = Distance/TTFoot;
                PWV_Foot(entryCount+2) = NaN;
                PWV_Foot(entryCount+3) = NaN;
            end 
            
            if get(handles.ttuRadio,'Value')
                PWV_Upstroke(entryCount+1) = Distance/TTUpstroke;
                PWV_Upstroke(entryCount+2) = NaN;
                PWV_Upstroke(entryCount+3) = NaN;
            end 
            
            if get(handles.xcorrRadio,'Value')
                PWV_Xcorr(entryCount+1) = Distance/Xcorr;
                PWV_Xcorr(entryCount+2) = NaN;
                PWV_Xcorr(entryCount+3) = NaN;
            end  
            
            PWV_Average(entryCount+1) = Distance/PWV_Average;
            PWV_Average(entryCount+2) = NaN;
            PWV_Average(entryCount+3) = NaN;
            PWV_Average = PWV_Average';
            
        end 
        
        PLANES{entryCount+1} = 'FIT_PWV';
        PLANES{entryCount+2} = 'm (slope)';
        PLANES{entryCount+3} = 'b (y-intercept)';
        for n = 1:3
            Distance(entryCount+n) = NaN;
            TTPeak(entryCount+n) = NaN;
            TTPoint(entryCount+n) = NaN;
            TTFoot(entryCount+n) = NaN;
            TTUpstroke(entryCount+n) = NaN;
            Xcorr(entryCount+n) = NaN;
        end 
        
        PLANES = PLANES';
        Distance = Distance';
        TTPeak = TTPeak';
        TTPoint = TTPoint';
        TTFoot = TTFoot';
        TTUpstroke = TTUpstroke';
        Xcorr = Xcorr';
        
        pwvTable = table(PLANES,Distance,TTPeak,TTPoint,TTFoot,TTUpstroke,Xcorr,PWV_Peak,PWV_Point,PWV_Foot,PWV_Upstroke,PWV_Xcorr,PWV_Average);
        if ~exist('Data-DataAnalysis','dir')
            mkdir 'Data-DataAnalysis'
        end
        cd 'Data-DataAnalysis'
        writetable(pwvTable,'Summary.xlsx','FileType','spreadsheet','Sheet',['Interpolation - ' interpType]);
        cd ..
        clear PLANES Distance TTPeak TTPoint TTFoot TTUpstroke Xcorr 
        clear PWV_Peak PWV_Point PWV_Foot PWV_Upstroke PWV_Xcorr PWV_Average
    end 
    
    guidata(hObject, handles);
    
    cd 'Data-DataAnalysis'
    flow = handles.flow;
    anatDatasets = handles.anatDatasets;
    pcDatasets = handles.pcDatasets;
    magDatasets = handles.magDatasets;
    save('anatDatasets.mat','anatDatasets')
    save('pcDatasets.mat','pcDatasets')
    save('magDatasets.mat','magDatasets')
    save('flow.mat','flow')
    cd ..
    
    dataDir = uigetdir('C:\','Save Data Location');
    movefile('CenterlineImages',dataDir)
    movefile('FlowWaveformImages',dataDir);
    movefile('PWVanalysisPlotImages',dataDir);
    movefile('ROIimages',dataDir);
    movefile('Data-DataAnalysis',dataDir);

    
% --- MESSAGE BAR - CALLBACK
function MessageBar_Callback(hObject, eventdata, handles)

% --- MESSAGE BAR - CREATE FUNCTION
function MessageBar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
   

% --- TTPeak READOUT - CREATE FUNCTION
function ttpData_CreateFcn(hObject, eventdata, handles)
% --- TTPeak RADIO - CALLBACK
function ttpRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttpRadio,'Value')
        set(handles.ttpData,'String',' ');
    end 
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTPoint READOUT - CREATE FUNCTION
function ttpointData_CreateFcn(hObject, eventdata, handles)
% --- TTPoint RADIO - CALLBACK
function ttpointRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttpointRadio,'Value')
        set(handles.ttpointData,'String',' ');
    end 
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTUpstroke READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)
% --- TTUpstroke RADIO - CALLBACK
function ttuRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',' ');
    end 
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTFoot READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)
% --- TTFoot RADIO - CALLBACK
function ttfRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',' ');
    end 
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- Xcorr READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)
% --- Xcorr RADIO - CALLBACK
function xcorrRadio_Callback(hObject, eventdata, handles)
    if ~get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',' ');
    end 
    if handles.global.startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
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
function updateImages(handles)
    axes(handles.PlanePlot);
    
    planeNum = get(handles.PlanePopup,'Value');
    datasetNum = get(handles.DatasetPopup,'Value');
    
    isPC = get(handles.PCRadio,'Value');
    if isPC
        dataset = handles.pcDatasets;
    else
        dataset = handles.magDatasets;
    end 

    if isstruct(dataset(planeNum).Data)
        imageSet = struct2cell(dataset(planeNum).Data);
    else
        imageSet = dataset(planeNum).Data;
    end 

    if iscell(imageSet)
        images = imageSet(datasetNum);
        images = cell2mat(images);
    else 
        images = imageSet;
    end 

    if ndims(images)<3
        steps = [1 500];
        set(handles.PlaneSlider, 'SliderStep', steps);
        slice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.PlaneSlider, 'SliderStep', steps);
        sliceNum = 1+round( get(handles.PlaneSlider,'Value').*(dim3size-1) );
        slice = images(:,:,sliceNum);
    end 

    if handles.global.zoomed
        slice = slice(handles.global.spanUD,handles.global.spanLR);
    end 

    if ~isempty(handles.pcDatasets(planeNum).ROI)
        hold on
        imshow(slice,[]);
    else 
        cla(handles.PlanePlot,'reset')
        imshow(slice,[])
    end 
    
    clear dataset datasetNum images imageSet isPC planeNum slice steps 
    
    
% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    axes(handles.AnatomicalPlot); %show image
    
    datasetNum = get(handles.AnatListbox,'Value');
    images = handles.anatDatasets(datasetNum).Data;

    if ndims(images)<3
        steps = [1 1];
        set(handles.SliderAnatomical, 'SliderStep', steps);
        anatSlice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.SliderAnatomical, 'SliderStep', steps);
        sliceNum = 1+round( get(handles.SliderAnatomical,'Value').*(dim3size-1) );
        anatSlice = images(:,:,sliceNum);
    end 
    
    if handles.global.zoomedAnat
        anatSlice = anatSlice(handles.global.spanUDAnat,handles.global.spanLRAnat);
    end
    
    if get(handles.ShowPlanesRadio,'Value') && handles.anatDatasets(datasetNum).sliceDim~=3
        imshow(anatSlice,[]);
        for i=1:numel(handles.flow)
            pcRow = handles.flow(i).planeLineRow;
            anatSlice(pcRow,:) = max(max(anatSlice))+100;
            hold on; plot(1:size(anatSlice,1),pcRow.*ones(1,size(anatSlice,2)),'y');
        end 
    else 
       imshow(anatSlice,[]); 
    end 
    set(handles.sliceNumText,'String',sliceNum);
    
    clear anatSlice datasetNum dim3size i images pcRow sliceNum spanLRAnat spanUDAnat steps

    
% --- "Time to" calculations (TTPeak, TTPoint, TTUpstroke, TTFoot, Xcorr)
function flow = computeTTs(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if handles.global.startAnalyzing
            switch handles.global.interpType
                case 'None'
                    flowTemp = flow(i).Data.meanROI;
                case 'Makima'
                    flowTemp = flow(i).Makima.meanROI;
                case 'Gaussian'
                    flowTemp = flow(i).Gaussian.meanROI;  
                case 'Spline'
                    flowTemp = flow(i).Spline.meanROI; 
                otherwise
            end 
        else
            flowTemp = flow(i).Data.meanROI;
        end 
        
        if mean(flowTemp)<0
            flowTemp = -1*flowTemp;
        else 
            flowTemp = flowTemp;
        end 
        
        if handles.global.pgShift
            flowTemp = circshift(flowTemp,round(length(flowTemp)/2));
        end 
                
        times = flow(i).Data.times;
        timeres = times(2)-times(1);
        [maxPeakVel,maxPeakVelIdx] = max(flowTemp);
        upstroke = flowTemp(1:maxPeakVelIdx);
        
        [~,SeventyPointIdx] = min(abs(upstroke-0.7*maxPeakVel));
        curvePoints(i).SeventyPointIdx = SeventyPointIdx;
        curvePoints(i).SeventyPoint = flowTemp(SeventyPointIdx);
        
        [~,ThirtyPointIdx] = min(abs(upstroke-0.3*maxPeakVel));
        curvePoints(i).ThirtyPointIdx = ThirtyPointIdx;
        curvePoints(i).ThirtyPoint = flowTemp(ThirtyPointIdx);
        
        flowTemp = normalize(flowTemp,'range');
        [~,FiftyPointIdx] = min(abs(upstroke-0.5*maxPeakVel));
        curvePoints(i).maxPeakVelIdx = maxPeakVelIdx;   
        curvePoints(i).maxPeakVel = maxPeakVel;        
        curvePoints(i).FiftyPointIdx = FiftyPointIdx;
        curvePoints(i).FiftyPoint = flowTemp(FiftyPointIdx);

        flows(i,:) = flowTemp;
    end 
    
    % TTPeak - time to peak calculation
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTPeak = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                TTPeak(iterator) = timeres.*( curvePoints(iterator).maxPeakVelIdx - curvePoints(i).maxPeakVelIdx );
            end 
        end 
        flow(i).TTPeak = TTPeak;
    end      
    
    % TTPoint - time to point calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTPoint = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                TTPoint(iterator) = timeres.*( curvePoints(iterator).FiftyPointIdx - curvePoints(i).FiftyPointIdx );
            end 
        end
        flow(i).TTPoint = TTPoint;
    end    
    
    % TTUpstroke - time to upstroke calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTUpstroke = [NaN,NaN,NaN];  
        
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                [sigmoid1,upslope1,~,tUpstroke1,~,timeresSigmoid1] = sigFit(flows(i,:));
                [sigmoid2,upslope2,~,tUpstroke2,~,timeresSigmoid2] = sigFit(flows(iterator,:));
                TTUpstroke(iterator) = timeres.*(timeresSigmoid1*tUpstroke2 - timeresSigmoid2*tUpstroke1);
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
                % m = (y2-y1)/(x2-x1); y1=80%max flow and y2=20%max flow
                % x1 and x2 are the times (indices) where these values occur
                m1 = (curvePoints(i).SeventyPoint-curvePoints(i).ThirtyPoint)/(curvePoints(i).SeventyPointIdx-curvePoints(i).ThirtyPointIdx);
                % t1 is the x-intercept of this line. 
                % This can be solved analytically; x0 = x1-y1/m
                t1 = curvePoints(i).ThirtyPointIdx - (curvePoints(i).ThirtyPoint/m1);
                % Do the same for the consequent flow curve
                m2 = (curvePoints(iterator).SeventyPoint-curvePoints(iterator).ThirtyPoint)/(curvePoints(iterator).SeventyPointIdx-curvePoints(iterator).ThirtyPointIdx);
                t2 = curvePoints(iterator).ThirtyPointIdx - (curvePoints(iterator).ThirtyPoint/m2);
                TTFoot(iterator) = timeres.*(t2-t1);
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
                XcorrPlot = xcorr(flows(i,:),flows(iterator,:));
                [~,maxXcorrIdx] = max(XcorrPlot);
                shift = -1*(maxXcorrIdx - length(flows(i,:)));
                Xcorr(iterator) = shift.*timeres;
            end 
        end 
        flow(i).Xcorr = Xcorr;
    end 
    
 
% --- Condense and organize flow data obtained from ROIs
function flow = organizeFlowInfo(handles)
    count = 1;
    for i=1:numel(handles.pcDatasets)
        if isstruct(handles.pcDatasets(i).ROIdata)
            if length(handles.pcDatasets(i).ROIdata)==1
                flow(count).Name = handles.pcDatasets(i).Names;
                flow(count).Data = handles.pcDatasets(i).ROIdata;
                flow(count).Info = handles.pcDatasets(i).Info;
                flow(count).Makima = handles.pcDatasets(i).ROIdataMakima;
                flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian;
                flow(count).Spline = handles.pcDatasets(i).ROIdataSpline;
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdata)
                    name = handles.pcDatasets(i).Names;
                    planeName = [name ' Plane ' num2str(j)];
                    flow(count).Name = planeName;
                    flow(count).Data = handles.pcDatasets(i).ROIdata(j);
                    flow(count).Info = handles.pcDatasets(i).Info;
                    flow(count).Makima = handles.pcDatasets(i).ROIdataMakima(j);
                    flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian(j);
                    flow(count).Spline = handles.pcDatasets(i).ROIdataSpline(j);
                    count = count+1;
                end 
            end 
        end 
    end 
    clear count i name planeName

    
% --- Turn PolyLine into SplineLine    
function Y = interppolygon(X,N)
    if nargin < 2 || N < 2
        N = 2;
    end
    nDim = size(X,2);
    dx = 0;
    
    for dim = 1:nDim
        dx = dx + diff(X(:,dim)).^2 ;
    end
    
    lengthBetweenPoints = sqrt(dx);
    lengthLine = sum(lengthBetweenPoints);
    origMetric = [0; cumsum(lengthBetweenPoints/lengthLine)];
    
    interpMetric = (0:(1/(N-1)):1)';
    Y = interp1(origMetric,X,interpMetric,'makima');
    %Y = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
    clear interpMetric origMetric lengthLine lengthBetweenPoints dx nDim
    
    
% --- Plot Velocities    
function plotVelocity(handles)
    cla(handles.VelocityPlot,'reset')
    axes(handles.VelocityPlot);
    
    legendSet = {'Baseline',handles.flow.Name};
    count = length(legendSet)-1;
    times = handles.flow(1).Data.times;
    
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5);
    for i=1:count
        switch handles.global.interpType
            case 'None'
                flow = handles.flow(i).Data.meanROI;
                stdev = handles.flow(i).Data.stdROI;
            case 'Makima'
                times = handles.flow(i).Makima.times;
                flow = handles.flow(i).Makima.meanROI;
                stdev = handles.flow(i).Makima.stdROI;
            case 'Gaussian'
                times = handles.flow(i).Gaussian.times;
                flow = handles.flow(i).Gaussian.meanROI;  
                stdev = handles.flow(i).Gaussian.stdROI;
            case 'Spline'
                times = handles.flow(i).Spline.times;
                flow = handles.flow(i).Spline.meanROI; 
                stdev = handles.flow(i).Spline.stdROI;
            otherwise
        end
        
        if handles.global.pgShift
            flow = circshift(flow,round(length(flow)/2));
            stdev = circshift(stdev,round(length(stdev)/2));
        end 
        
        if mean(flow)<0
            if handles.global.showErrorBars
                hold on; errorbar(times,-1*flow,stdev);
            else
                hold on; plot(times,-1*flow);
            end 
        else 
            if handles.global.showErrorBars
                hold on; errorbar(times,flow,stdev);
            else
                hold on; plot(times,flow);
            end
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)');
    xlim([times(1),times(end)]);
    
    guidata(hObject,handles)
    clear legendSet count times flow stdev i
    

% --- Sigmoid Fit Function
function [sigmoid,upslope,curvature,t1,t2,timeres] = sigFit(meanROI)
    [~,peak] = max(meanROI);
    upslope = meanROI(1:peak);
    t = 1:peak;
    times = 1:0.1:peak;
    t2 = length(times);
    upslope = interp1(t,upslope,times);
    upslope = normalize(upslope,'range');

    timeres = times(2)-times(1);

    %c1 = b, c2 = a, c3 = x0, c4 = dx
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((times-c(3))./c(4)) ) - upslope;
    c0 = [max(upslope),min(upslope),peak/2,timeres/2];
    opts = optimset('Display', 'off');
    c = lsqnonlin(sigmoidModel,c0,[],[],opts);
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((times-c(3))./c(4)) );

    for i=2:numel(upslope)-1
        dt = 0.5*(times(i+1)-times(i-1));
        ddt = 4*(times(i+1)-2*times(i)+times(i-1));
        dy = 0.5*(sigmoid(i+1)-sigmoid(i-1));
        ddy = 4*(sigmoid(i+1)-2*sigmoid(i)+sigmoid(i-1));
        curvature(i) = (ddy*dt - ddt*dy)./((dt^2 + dy^2).^(3/2));
    end 
    [~,t1] = max(curvature);
    %save upslope sigmoid curvature plots
    clear peak t times sigmoidModel c0 opts c 
    clear dt ddt dy ddy 
    
    
% --- Calculate derivative   
function fPrime = derivative(f)
    f(end+1) = f(1);
    fPrime = diff(f);
   
