function varargout = LoadPWV(varargin)
    % LOADPWV MATLAB code for LoadPWV.fig
    %      LOADPWV, by itself, creates a new LOADPWV or raises the existing
    %      singleton*.
    %
    %      H = LOADPWV returns the handle to a new LOADPWV or the handle to
    %      the existing singleton*.
    %
    %      LOADPWV('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in LOADPWV.M with the given input arguments.
    %
    %      LOADPWV('Property','Value',...) creates a new LOADPWV or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before LoadPWV_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to LoadPWV_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help LoadPWV

    % Last Modified by GUIDE v2.5 22-Jul-2020 09:22:22

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @LoadPWV_OpeningFcn, ...
                       'gui_OutputFcn',  @LoadPWV_OutputFcn, ...
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


% --- Executes just before LoadPWV is made visible.
function LoadPWV_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to LoadPWV (see VARARGIN)

    % Choose default command line output for LoadPWV
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    
    
    % Initialize 'global' variables into handles struct
    handles.anatDatasets.Names = ''; %list of anatomical dataset names
    handles.anatDatasets.Data = 0; %space for image data
    handles.global.anatDataIter = 1; %index of dataset (in handles) 
    
    handles.pcDatasets.Names = ''; %list of PC velocity dataset names
    handles.pcDatasets.Data = 0; %space for image data
    handles.global.pcDataIter = 1; %index of dataset (in handles) 
    
    handles.magDatasets.Names = ''; %list of PC magnitude dataset names
    handles.magDatasets.Data = 0; %space for image data
    handles.global.magDataIter = 1; %index of dataset (in handles) 

    guidata(hObject, handles);
    
% --- Outputs from this function are returned to the command line.
function varargout = LoadPWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;

    
    
%%%%%%%%%%%%%% PRELOAD CASE PANEL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PRELOAD CASE CALLBACK
function PreloadCaseButton_Callback(hObject, eventdata, handles)
    [casefile, casefileDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the case file (NAME_pwvcasefile_DATE.mat)');
    load([casefileDir '\' casefile]); %load casefile
    cd(casefileDir);
    
    % Place loaded datasets back into current handles
    handles.anatDatasets = anatDatasets; 
    handles.pcDatasets = pcDatasets;
    handles.magDatasets = magDatasets;
    
    % Set selectable listbox strings to dataset names
    set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names});
    set(handles.Listbox2DPC,'String',{handles.pcDatasets.Names});
    set(handles.Listbox2DMAG,'String',{handles.magDatasets.Names});
    set(handles.MessageBar,'String','Case file loaded successfully')

    handles.global.anatDataIter = numel(anatDatasets) + 1;
    % If only 1 entry (=0), no data was loaded. Subtract this from the iter
    if numel(anatDatasets(1).Data)==1 
        handles.global.anatDataIter = handles.global.anatDataIter - 1;
    end 
    
    handles.global.pcDataIter = numel(pcDatasets) + 1;
    if ~isstruct(pcDatasets(1).Data(1))
        handles.global.pcDataIter = handles.global.pcDataIter - 1;
    end 
    
    handles.global.magDataIter = numel(magDatasets) + 1;
    if ~isstruct(magDatasets(1).Data(1))
        handles.global.magDataIter = handles.global.magDataIter - 1;
    end 
      
    guidata(hObject, handles);


    
%%%%%%%%%%%%%% ANATOMICAL PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME ANATOMICAL CALLBACK
function DataNameAnatomical_Callback(hObject, eventdata, handles)

% --- DATA NAME ANATOMICAL CREATE FUNCTION
function DataNameAnatomical_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- BUTTON LOAD ANATOMICAL CALLBACK
function ButtonLoadAnatomical_Callback(hObject, eventdata, handles)
    anatDataIter = handles.global.anatDataIter; %grab number for listbox
    set(handles.MessageBar,'String',''); %erase message bar
    if isempty(get(handles.DataNameAnatomical,'String')) %if no name typed
        set(handles.MessageBar,'String', ... %type name in text field
            'Name is required to load data! Type a name into the text box and press the Load Data button.')
    else 
        name = get(handles.DataNameAnatomical,'String'); %grab name
        % Select matlab or dicom file of anatomical image(s)
        [anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE anatomical file in the dataset');
        [~,~,extension] = fileparts(anatomicalFile); %get file extension
        dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
        if isequal(extension,'.dcm') %if our extension is a dicom file
            handles.anatDatasets(anatDataIter).Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
            for i=1:length(dirInfo) %read all dcm files
                data(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
            end        
        else %if a single matlab file (with all images)
            temp = load(fullfile(anatomicalDir,dirInfo(1).name));
            data = single(struct2array(temp)); %typecast to single
        end
        
        % Check if proposed dataset name is already used with strcmp
        if sum(strcmp({handles.anatDatasets.Names},name))>0 %if used
            set(handles.MessageBar,'String', ...
                'This name has already been used, try loading again with a different name')
            set(handles.DataNameAnatomical,'String','') %erase name
        else %else, we're good to go
            set(handles.MessageBar,'String', ...
                ['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            handles.anatDatasets(anatDataIter).Names = name; %add name to list
            handles.anatDatasets(anatDataIter).Data = data; %add data to handles
            handles.anatDatasets(anatDataIter).RootDir = anatomicalDir;
            handles.global.anatDataIter = anatDataIter + 1; %move iter up 1
            set(handles.DataNameAnatomical,'String','') %erase name
            set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names}); %update listbox names
        end 
    end 

    guidata(hObject, handles);

% --- BUTTON LOAD ANATOMICAL CREATE FUNCTION
function ButtonLoadAnatomical_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX ANATOMICAL CALLBACK
function ListboxAnatomical_Callback(hObject, eventdata, handles)

% --- LISTBOX ANATOMICAL CREATE FUNCTION
function ListboxAnatomical_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE ANATOMICAL CALLBACK
function RemoveAnatomical_Callback(hObject, eventdata, handles)
    index = get(handles.ListboxAnatomical, 'value'); %get current index to remove
    handles.anatDatasets(index) = []; %remove everything at index
    set(handles.ListboxAnatomical,'Value',1);
    set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names}); %reset listbox names
    handles.global.anatDataIter = handles.global.anatDataIter - 1; %bring back iter by 1
    
    guidata(hObject, handles);

% --- REMOVE ANATOMICAL CREATE FUNCTION
function RemoveAnatomical_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%%%% 2DPC PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 2D PC DATA NAME CALLBACK
function DataName2DPC_Callback(hObject, eventdata, handles)

% --- 2D PC DATA NAME CREATE FUNCTION
function DataName2DPC_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- BUTTON LOAD 2DPC CALLBACK
function ButtonLoad2DPC_Callback(hObject, eventdata, handles)
    pcDataIter = handles.global.pcDataIter; %index in handles
    set(handles.MessageBar,'String','')
    if isempty(get(handles.DataName2DPC,'String'))
        set(handles.MessageBar,'String', ...
            'Name is required to load data! Type a name into the text box and press the Load Data button.')
    elseif sum(strcmp({handles.pcDatasets.Names},get(handles.DataName2DPC,'String')))>0
        set(handles.MessageBar,'String', ...
            'This name has already been used, try loading again with a different name')
        set(handles.DataName2DPC,'String','')
    else 
        name = get(handles.DataName2DPC,'String');
        [pcFile, pcDir] = uigetfile({'*.dcm;*.dat;','Useable Files (*.dcm,*.dat)';
       '*.dcm', 'DICOM files (*.dcm)'; ...
       '*.dat', 'Data files (*.dat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2D phase contrast file in the dataset');
        [~,~,extension] = fileparts(pcFile);
        dirInfo = dir(fullfile(pcDir,['*' extension]));
        if isequal(extension,'.dcm')
            handles.pcDatasets(pcDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name)); %get dicom metadata (from 1st dicom)
            for i=1:length(dirInfo)
                temp(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name))); %read dicoms and cast to single
            end  
            mag = temp(:,:,floor(length(dirInfo)/2)+1:end); %magnitude is last half
            v = temp(:,:,1:floor(length(dirInfo)/2)); %velocity is first half of images
            data.MAG = mean(mag,3); %time-averaged magnitude
            %data.VMEAN = BGPhaseCorrect(data.v);
            data.VMEAN = mean(v,3); %time-averaged velocity
            data.mag = mag;
            data.v = v;
        else
            fid = fopen([pcDir '\pcvipr_header.txt'], 'r'); %open header
            dataArray = textscan(fid,'%s%s%[^\n\r]','Delimiter',' ', ...
                'MultipleDelimsAsOne',true,'ReturnOnError',false); %parse header info
            fclose(fid);
            dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:),'UniformOutput',false);
            pcviprHeader = cell2struct(dataArray{1,2}(:),dataArray{1,1}(:),1); %turn to structure
            handles.pcDatasets(pcDataIter).Info = pcviprHeader; %add pcvipr header to handles
            resx = pcviprHeader.matrixx; %resolution in x
            resy = pcviprHeader.matrixy; %resolution in y
            nframes = pcviprHeader.frames; %number of cardiac frames
            MAG = load_dat(fullfile(pcDir,'MAG.dat'),[resx resy]); %Average magnitude
            CD = load_dat(fullfile(pcDir,'CD.dat'),[resx resy]); %Average complex difference
            VMEAN = load_dat(fullfile(pcDir,'comp_vd_3.dat'),[resx resy]); %Average velocity
            
            % Initialize data time-resolved data arrays
            mag = zeros(resx,resy,nframes); %Time-resolved magnitude
            cd = zeros(resx,resy,nframes); %Time-resolved complex difference
            v = zeros(resx,resy,nframes); %Time-resolved velocity 
            for j = 1:nframes  %velocity is placed in v3 for 2D (through-plane)
                mag(:,:,j) = load_dat(fullfile(pcDir,['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
                cd(:,:,j) = load_dat(fullfile(pcDir,['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
                v(:,:,j) = load_dat(fullfile(pcDir,['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
            end
            data.MAG = flipud(MAG);
            data.CD = flipud(CD);
            data.VMEAN = flipud(VMEAN);
            data.mag = flipud(mag);
            data.cd = flipud(cd);
            data.v = flipud(v);
        end
        set(handles.MessageBar,'String', ...
            ['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
        handles.pcDatasets(pcDataIter).Names = name; %add to name to namelist in handles
        handles.pcDatasets(pcDataIter).Data = data; %add data to handles
        handles.pcDatasets(pcDataIter).RootDir = pcDir;
        handles.global.pcDataIter = pcDataIter + 1; %move up index for next dataset
        set(handles.DataName2DPC,'String','');
        set(handles.Listbox2DPC,'String',{handles.pcDatasets.Names}); %update listbox names
    end 
    
    guidata(hObject, handles);

% --- BUTTON LOAD 2DPC CREATE FUNCTION
function ButtonLoad2DPC_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX 2DPC CALLBACK
function Listbox2DPC_Callback(hObject, eventdata, handles)

% --- LISTBOX 2DPC CREATE FUNCTION
function Listbox2DPC_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE 2DPC CALLBACK
function Remove2DPC_Callback(hObject, eventdata, handles)
    index = get(handles.Listbox2DPC, 'value'); %get current index to remove
    handles.pcDatasets(index) = []; %remove everything at index
    set(handles.Listbox2DPC,'Value',1);
    set(handles.Listbox2DPC,'String',{handles.pcDatasets.Names}); %reset listbox
    handles.global.pcDataIter = handles.global.pcDataIter - 1; %bring index back 1
    
    guidata(hObject, handles);
        
% --- REMOVE 2DPC CREATE FUNCTION
function Remove2DPC_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%%%% 2DMAG PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME 2DMAG CALLBACK
function DataName2DMAG_Callback(hObject, eventdata, handles)

% --- DATA NAME 2DMAG CREATE FUNCTION
function DataName2DMAG_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- BUTTON LOAD 2DMAG CALLBACK
function ButtonLoad2DMAG_Callback(hObject, eventdata, handles)
    magDataIter = handles.global.magDataIter;
    set(handles.MessageBar,'String','')

    if isempty(get(handles.DataName2DMAG,'String'))
        set(handles.MessageBar,'String', ...
            'Name is required to load data! Type a name into the text box and press the Load Data button.')
    else 
        name = get(handles.DataName2DMAG,'String');
        [magFile, magDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2D magnitude file in the dataset');
        [~,~,extension] = fileparts(magFile);
        dirInfo = dir(fullfile(magDir,['*' extension]));
        
        if isequal(extension,'.dcm')
            handles.magDatasets(magDataIter).Info = dicominfo(fullfile(magDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(magDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(magDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end
        
        if sum(strcmp({handles.magDatasets.Names},name))>0
            set(handles.MessageBar,'String', ...
                'This name has already been used, try loading again with a different name')
            set(handles.DataName2DMAG,'String','')
        else 
            set(handles.MessageBar,'String', ...
                ['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            handles.magDatasets(magDataIter).Names = name;
            handles.magDatasets(magDataIter).Data = data;
            handles.magDatasets(magDataIter).RootDir = magDir;
            handles.global.magDataIter = magDataIter + 1;
            set(handles.DataName2DMAG,'String','')
            set(handles.Listbox2DMAG,'String',{handles.magDatasets.Names});
        end 
    end 
    
    guidata(hObject, handles);

% --- BUTTON LOAD 2DMAG CREATE FUNCTION
function ButtonLoad2DMAG_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX 2DMAG CALLBACK
function Listbox2DMAG_Callback(hObject, eventdata, handles)

% --- LISTBOX 2DMAG CREATE FUNCTION
function Listbox2DMAG_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE 2DMAG CALLBACK
function Remove2DMAG_Callback(hObject, eventdata, handles)
    index = get(handles.Listbox2DMAG, 'value');
    handles.global.magDatasets(index) = [];
    set(handles.Listbox2DMAG,'Value',1);
    set(handles.Listbox2DMAG,'String',{handles.global.magDatasets.Names});
    handles.global.magDataIter = handles.global.magDataIter - 1;
    
    guidata(hObject, handles);

% --- REMOVE 2DMAG CREATE FUNCTION
function Remove2DMAG_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%%%% POST LOADING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- CASENAME CALLBACK
function CaseName_Callback(hObject, eventdata, handles)

% --- CASENAME CREATE FUNCTION
function CaseName_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SAVE CASE BUTTON CALLBACK
function SaveCaseButton_Callback(hObject, eventdata, handles)
    anatDatasets = handles.anatDatasets; %make as variable (for saving) 
    pcDatasets = handles.pcDatasets;
    magDatasets = handles.magDatasets;
    
    if isempty(get(handles.CaseName,'String')) %need a name
        set(handles.MessageBar,'String', ...
            'Name is required to save data! Type a name into the text box and press the Save Case File button.')
    else
        set(handles.MessageBar,'String', ...
            'Select the directory where you will save the Case File.');
        casefileDir = uigetdir(pwd,'Select the directory where you will save the Case File.'); %get save directory
        casename = get(handles.CaseName,'String'); %pull name from textfield
        date = datestr(now); %get current date/time
        chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)]; %chop date up
        filename = [casename '_pwvcasefile_' chopDate '.mat']; %combine name and date
        set(handles.MessageBar,'String',['Saving Data as "' filename '"']);
        filenameWithPath = [casefileDir '\' filename];
        save(filenameWithPath,'anatDatasets','pcDatasets','magDatasets'); %save anat, pc, and mag datasets
        set(handles.MessageBar,'String','Data saved successfully');
        set(handles.CaseName,'String','');
    end 
    
    guidata(hObject, handles);

    
% --- MESSAGE BAR CALLBACK
function MessageBar_Callback(hObject, eventdata, handles)

% --- MESSAGE BAR CREATE FUNCTION
function MessageBar_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- COMPLETE LOADING BUTTON CALLBACK
function ButtonCompleteLoading_Callback(hObject, eventdata, handles)
    if handles.global.anatDataIter==1 %if we haven't loaded anatomical data (iter still at 1)
        set(handles.MessageBar,'String', ...
            'Anatomical datasets required to complete loading process.'); %don't move on
    elseif handles.global.pcDataIter==1 % or if we haven't loaded pc data (iter still at 1)
        set(handles.MessageBar,'String', ...
            '2D PC datasets required to complete loading process.'); %don't move on
    else
        set(handles.MessageBar,'String','Loading Process Complete!');
        AnalyzePWV(handles.anatDatasets,handles.pcDatasets,handles.magDatasets) %MOVE TO NEXT GUI
    end 

% --- COMPLETE LOADING BUTTON CREATE FUNCTION
function ButtonCompleteLoading_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Dat files
function v = load_dat(name, res)
    [fid,errmsg]= fopen(name,'r');
    if fid < 0  %if name does not exist in directory
        set(handles.MessageBar,'String',['Error Opening Data : ',errmsg]);
    end

    % Reads in as short, reshapes by image res.
    v = reshape(fread(fid,'short=>single'),res);
    fclose(fid);

    
% Perform 2D background phase correction
function vmean = BGPhaseCorrect(v)
vmean = mean(abs(v),3);
vmax = max(abs(v),[],3);

cdThresh = 0.08*max(vmax(:));
mask = vmax<cdThresh;
[rows,cols] = find(mask);
vmean = double(vmean.*mask);
vmeanVector = nonzeros(vmean(:));

f = fit( [rows, cols], vmeanVector, 'poly33' );
c = coeffvalues(f);

x0 = 1:size(mask,1);
y0 = 1:size(mask,2);
[y,x] = meshgrid(x0,y0);

fits = c(1) + c(2).*x + c(3).*y + c(4).*(x.^2) + c(5).*x.*y + c(6).*(y.^2) ...
     + c(7).*(x.^3) + c(8).*(x.^2).*y + c(9).*x.*(y.^2) + c(10).*(y.^3); 
vmean = vmean-fits;





%%%%%%%%%%%%%% UNUSED %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- MAIN PANEL CREATE FUNCTION
function MainPanel_CreateFcn(hObject, eventdata, handles)
% --- MAIN PANEL DELETE FUNCTION
function MainPanel_DeleteFcn(hObject, eventdata, handles)
% --- HEADER MAIN CREATE FUNCTION
function HeaderMain_CreateFcn(hObject, eventdata, handles)
% --- HEADER ANATOMICAL CREATEFUNCTION
function HeaderAnatomical_CreateFcn(hObject, eventdata, handles)
% --- HEADER ANATOMICAL LISTBOX CREATEFUNCTION
function HeaderAnatomicalListing_CreateFcn(hObject, eventdata, handles)
% --- LOAD ANATOMICAL TEXT CREATEFUNCTION
function LoadAnatomicalText_CreateFcn(hObject, eventdata, handles)
% --- 2D PC HEADER CREATE FUNCTION
function Header2DPC_CreateFcn(hObject, eventdata, handles)
% --- LOAD 2DPC TEXT CREATE FUNCTION
function Load2DPCText_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DPC LISTBOX CREATE FUNCTION
function Header2DPCListing_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DMAG CREATE FUNCTION
function Header2DMAG_CreateFcn(hObject, eventdata, handles)
% --- LOAD 2DMAG TEXT CREATE FUNCTION
function Load2DMAGText_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DMAG LISTING CREATE FUNCTION
function Header2DMAGListing_CreateFcn(hObject, eventdata, handles)
