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

    % Last Modified by GUIDE v2.5 18-Oct-2019 21:52:16

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

    global zoomed spanUD spanLR
    handles.pcDatasets = varargin{1};
    handles.magDatasets = varargin{2};
    handles.pcDatasets(1).ROI = [];
    handles.magDatasets(1).ROI = [];
    handles.pcDatasets(1).ROIdata = [];
    zoomed = 0;
    spanUD = [0,0];
    spanLR = [0,0];
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    guidata(hObject, handles);
    updateImages(handles);
    if numel(handles.magDatasets)==0
        set(handles.MAGRadio,'Enable','off')
    end 

% --- Outputs from this function are returned to the command line.
function varargout = AnalyzePWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;

    
% --- ZOOM BUTTON - CALLBACK
function ZoomButton_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    disp('Draw an ROI to zoom in');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4));
    spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));
    zoomed = 1;
    updateImages(handles);

    
% --- UNZOOM BUTTON - CALLBACK
function UnzoomButton_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    zoomed = 0; spanUD = 0; spanLR = 0;
    updateImages(handles);

    
% --- DRAWROI BUTTON - CALLBACK
function DrawROIbutton_Callback(hObject, eventdata, handles)
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
                              disp('User pressed the escape key. Deleting ROI.')
                              DeleteROIbutton_Callback(hObject, eventdata, handles)
                              break % break out of the while loop
                          case 13 % 13 is the return key 
                              disp('ROI selected');
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
    
    if get(handles.PCRadio,'Value')
        set(handles.MAGRadio,'Enable','off');
    else 
        set(handles.PCRadio,'Enable','off');
    end 
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.PlanePopup,'Enable','off');
    set(handles.ZoomButton,'Enable','off');
    set(handles.UnzoomButton,'Enable','off');
    
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('An ROI has already been placed!\n');
    end 

    
% --- DELETE ROI BUTTON - CALLBACK
function DeleteROIbutton_Callback(hObject, eventdata, handles)
    hold off
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        handles.pcDatasets(planeNum).ROI = [];
    else 
        handles.magDatasets(planeNum).ROI = [];
    end 
    guidata(hObject,handles);
    updateImages(handles);
    set(handles.DrawROIbutton,'Enable','on');
    set(handles.PlanePopup,'Enable','on');
        if get(handles.PCRadio,'Value')
            set(handles.MAGRadio,'Enable','on');
        else 
            set(handles.PCRadio,'Enable','on');
        end 
    set(handles.ZoomButton,'Enable','on');
    set(handles.UnzoomButton,'Enable','on');

    
% --- ANALYZE ROI BUTTON - CALLBACK
function AnalyzeROIbutton_Callback(hObject, eventdata, handles)
    global spanUD spanLR zoomed
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        v = handles.pcDatasets(planeNum).Data.v;
        if zoomed
            v = v(spanUD,spanLR,:);
        end 
        circle = handles.pcDatasets(planeNum).ROI;
        radius = circle.Radius; 
        center = round(circle.Center);

        [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
        X = X-center(2); %shift coordinate grid
        Y = Y-center(1);
        roiMask = sqrt(X.^2+Y.^2)<=radius;

        for i=1:size(v,3)
            vTemp = v(:,:,i);
            roiDataRaw(:,i) = vTemp(roiMask);
            meanROI(i) = mean(vTemp(roiMask));
            medianROI(i) = median(vTemp(roiMask));
            maxROI(i) = max(vTemp(roiMask));
            minROI(i) = min(vTemp(roiMask));
            stdROI(i) = std(vTemp(roiMask));
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
            area = sum(roiMask(:))*(xres)^2;
            flowROI(i) = area.*meanROI(i);
        end         
        times = timeres.*(1:length(flowROI));
        dFlow = derivative(abs(flowROI));
        dMean = derivative(abs(meanROI));
        roiStatistics.radius = radius;
        roiStatistics.center = center;
        roiStatistics.roiMask = roiMask;
        roiStatistics.times = times;
        roiStatistics.roiDataRaw = roiDataRaw;
        roiStatistics.meanROI = meanROI;
        roiStatistics.dMean = dMean;
        roiStatistics.medianROI = medianROI;
        roiStatistics.maxROI = maxROI;
        roiStatistics.minROI = minROI;
        roiStatistics.stdROI = stdROI;
        roiStatistics.flowROI = flowROI;
        roiStatistics.dFlow = dFlow;
        if isstruct(handles.pcDatasets(planeNum).ROIdata)
            handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
        else 
            handles.pcDatasets(planeNum).ROIdata = roiStatistics;
        end 
        clear times dFlow dMean radius center roiMask
        clear roiDataRaw meanROI medianROI maxROI minROI stdROI flowROI
    end 
    
    % now start plotting flow curves
    axes(handles.VelocityPlot)
    handles.flow = organizeFlowInfo(handles);
    legendSet = {'Baseline',flow.Name};
    count = length(legendSet)-1;
    if count==1
        hold on; plot( zeros(1,length(roiStatistics.flowROI)) ,'Color','black','LineWidth',2);
    end 
 
    if mean(roiStatistics.flowROI)<0
        hold on; plot(-1*roiStatistics.flowROI);
    else 
        hold on; plot(roiStatistics.flowROI);
    end 
    legend(legendSet);
    
    clear roiStatistics
    guidata(hObject,handles);
    DeleteROIbutton_Callback(hObject, eventdata, handles);
    updateImages(handles);
    axes(handles.PlanePlot);

    
% --- PLANE DROPDOWN - CALLBACK
function PlanePopup_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    zoomed = 0; spanUD = 0; spanLR = 0;
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(planeNum).Data));
    end 
    set(handles.DatasetPopup,'Value',1);
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

    
    
% --- PLANE SLIDER - CALLBACK
function PlaneSlider_Callback(hObject, eventdata, handles)
    updateImages(handles);

% --- PLANE SLIDER - CREATE FUNCTION
function PlaneSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --- PC RADIO CLICK - CALLBACK
function PCRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'Enable','on');
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.MAGRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    updateImages(handles);


% --- MAG RADIO CLICK - CALLBACK
function MAGRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.magDatasets.Names});
    set(handles.DatasetPopup,'Enable','off');
    set(handles.PCRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    set(handles.DatasetPopup,'String',' ');
    updateImages(handles);
    

% --- ANATOMICAL SLIDER - CALLBACK
function SliderAnatomical_Callback(hObject, eventdata, handles)


% --- ANATOMICAL SLIDER - CREATE FUNCTION
function SliderAnatomical_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --- ANATOMICAL LISTBOX - CALLBACK
function AnatListbox_Callback(hObject, eventdata, handles)


% --- ANATOMICAL LISTBOX - CREATE FUNCTION
function AnatListbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- ZOOM ANATOMICAL - CALLBACK
function ZoomAnatomicalButton_Callback(hObject, eventdata, handles)


% --- UNZOOM ANATOMICAL - CALLBACK
function UnzoomAnatomicalButton_Callback(hObject, eventdata, handles)


% --- DRAW CENTERLINE - CALLBACK
function DrawCenterlineButton_Callback(hObject, eventdata, handles)


% --- DELETE CENTERLINE - CALLBACK
function DeleteCenterlineButton_Callback(hObject, eventdata, handles)


% --- COMPUTE CENTERLINE - CALLBACK
function ComputeCenterlineButton_Callback(hObject, eventdata, handles)

% --- SHOW PLANE RADIO - CALLBACK
function ShowPlanesRadio_Callback(hObject, eventdata, handles)

    
% --- PLANE PLOT - CREATE FUNCTION
function PlanePlot_CreateFcn(hObject, eventdata, handles)

% --- VELOCITY PLOT - CREATE FUNCTION
function VelocityPlot_CreateFcn(hObject, eventdata, handles)

% --- ANATOMICAL PLOT - CREATE FUNCTION
function AnatomicalPlot_CreateFcn(hObject, eventdata, handles)


%%%%%%%  MY FUNCTIONS  %%%%%%%%%

% --- Update Images in PLANE PLOT
function varargout = updateImages(varargin)
    global spanUD spanLR zoomed 
    handles = varargin{1};
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
        steps = [1 1];
        set(handles.PlaneSlider, 'SliderStep', steps);
        slice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.PlaneSlider, 'SliderStep', steps);
        sliceNum = 1+round( get(handles.PlaneSlider,'Value').*(dim3size-1) );
        slice = images(:,:,sliceNum);
    end 

    if zoomed
        slice = slice(spanUD,spanLR);
    end 

    varargout{1} = slice;   
    axes(handles.PlanePlot);
    if ~isempty(handles.pcDatasets(planeNum).ROI)
        hold on
        imshow(slice,[]);
    else 
        cla(handles.PlanePlot,'reset')
        imshow(slice,[])
    end 
    
% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    datasetNum = get(handles.DatasetListbox,'Value');
    imageSet = handles.anatDatasets(datasetNum).Data;

    dim3size = size(imageSet,3);
    steps = [1/(dim3size-1) 10/(dim3size-1)];
    set(handles.ImageSlider, 'SliderStep', steps);
    sliceNum = 1+round( get(handles.SliderAnatomical,'Value').*(dim3size-1) );
    anatSlice = imageSet(:,:,sliceNum);

    axes(handles.AnatomicalPlot); %show image
    imshow(anatSlice,[]);

 
% --- Calculate derivative   
function fPrime = derivative(f)
    f(end+1) = f(1);
    fPrime = diff(f);
    
    
% --- TTP - time to peak calculation
function ttp = computeTTP(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        [~,idx] = max(flow(i).Data.flowROI);
        ttp(i).maxFlowIdx = idx;     
    end 
    
    timeres = flow(i).Data.times(1);
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        comparer = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            comparer(iterator) = timeres.*( ttp(iterator).maxFlowIdx - ttp(i).maxFlowIdx );
        end 
        ttp(i).deltaT = comparer;
    end 

    
% --- TTF - time to foot calculation
function ttf = computeTTF(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        
    end 

    
% --- TTU - time to upstroke calculation    
function ttu = computeTTU(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        
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
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdata)
                    name = handles.pcDatasets(i).Names;
                    planeName = [name ' Plane ' num2str(j)];
                    flow(count).Name = planeName;
                    flow(count).Data = handles.pcDatasets(i).ROIdata(j);
                    flow(count).Info = handles.pcDatasets(i).Info;
                    count = count+1;
                end 
            end 
        end 
    end 

    
% --- Executes on selection change in InterpolatePopup.
function InterpolatePopup_Callback(hObject, eventdata, handles)
    global interp
    interpType = get(handles.InterpolatePopup,'Value');
    switch interpType
        case 2 
            interp = 2;
        case 3 
            interp = 3;
        case 4
            interp = 'Spline';
        otherwise
    end 


function InterpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(handles.InterpolatePopup,'String',{'1x','2x','3x','Spline'});



