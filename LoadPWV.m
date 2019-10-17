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

% Last Modified by GUIDE v2.5 10-Oct-2019 15:08:11

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

% UIWAIT makes LoadPWV wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadPWV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%% PRELOAD CASE PANEL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in PreloadCaseButton.
function PreloadCaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreloadCaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anatDataIter pcDataIter magDataIter pcviprDataIter preloaded
[casefile, casefileDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, 'Select the case file (NAME_pwvcasefile_DATE.mat)');
load([casefileDir '\' casefile]);
set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
set(handles.Listbox2DPC,'String',{pcDatasets.Names});
set(handles.Listbox2DMAG,'String',{magDatasets.Names});
set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});

disp('Case file loaded successfully');
set(handles.ErrorMessageBar,'String','Case file loaded successfully')

anatDataIter = numel(anatDatasets);
pcDataIter = numel(pcDatasets);
magDataIter = numel(magDatasets);
pcviprDataIter = numel(pcviprDatasets);
preloaded=1;


%%%%%%%%%%%%%% ANATOMICAL PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function HeaderAnatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HeaderAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function LoadAnatomicalText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadAnatomicalText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function DataNameAnatomical_Callback(hObject, eventdata, handles)
% hObject    handle to DataNameAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataNameAnatomical as text
%        str2double(get(hObject,'String')) returns contents of DataNameAnatomical as a double

% --- Executes during object creation, after setting all properties.
function DataNameAnatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataNameAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global anatDatasets anatDataIter preloaded
if isempty(preloaded)
    anatDatasets.Names = '';
    anatDatasets.Data = 0;
    anatDataIter = 1;
end 

% --- Executes on button press in ButtonLoadAnatomical.
function ButtonLoadAnatomical_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anatDatasets anatDataIter
set(handles.ErrorMessageBar,'String','')
if isempty(get(handles.DataNameAnatomical,'String'))
    disp('Name is required to load data!');
    disp('Type a name into the text box and press the Load Data button.');
    set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
else 
    name = get(handles.DataNameAnatomical,'String');
    [anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE anatomical file in the dataset');
    [~,~,extension] = fileparts(anatomicalFile);
    dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
    if isequal(extension,'.dcm')
        anatDatasets(anatDataIter).Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
        for i=1:length(dirInfo)
            data(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
        end        
    else
        temp = load(fullfile(anatomicalDir,dirInfo(1).name));
        data = single(struct2array(temp));
    end
    
    if sum(strcmp({anatDatasets.Names},name))>0
        set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
        set(handles.DataNameAnatomical,'String','')
    else 
        set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
        anatDatasets(anatDataIter).Names = name;
        anatDatasets(anatDataIter).Data = data;
        anatDataIter = anatDataIter+1;
        set(handles.DataNameAnatomical,'String','')
        set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
    end 
end 

% --- Executes during object creation, after setting all properties.
function ButtonLoadAnatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonLoadAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function HeaderAnatomicalListing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HeaderAnatomicalListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in ListboxAnatomical.
function ListboxAnatomical_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListboxAnatomical contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListboxAnatomical

% --- Executes during object creation, after setting all properties.
function ListboxAnatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListboxAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RemoveAnatomical.
function RemoveAnatomical_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global anatDatasets anatDataIter
index = get(handles.ListboxAnatomical, 'value');
anatDatasets(index) = [];
set(handles.ListboxAnatomical,'Value',1);
set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
anatDataIter = anatDataIter-1;

% --- Executes during object creation, after setting all properties.
function RemoveAnatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RemoveAnatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



%%%%%%%%%%%%%% 2DPC PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function Header2DPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Load2DPCText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Load2DPCText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function DataName2DPC_Callback(hObject, eventdata, handles)
% hObject    handle to DataName2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataName2DPC as text
%        str2double(get(hObject,'String')) returns contents of DataName2DPC as a double

% --- Executes during object creation, after setting all properties.
function DataName2DPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataName2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global pcDatasets pcDataIter preloaded
if isempty(preloaded)
    pcDatasets.Names = '';
    pcDatasets.Data = 0;
    pcDataIter = 1;
end

% --- Executes on button press in ButtonLoad2DPC.
function ButtonLoad2DPC_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoad2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pcDatasets pcDataIter
set(handles.ErrorMessageBar,'String','')
if isempty(get(handles.DataName2DPC,'String'))
    disp('Name is required to load data!');
    disp('Type a name into the text box and press the Load Data button.');
    set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
elseif sum(strcmp({pcDatasets.Names},get(handles.DataName2DPC,'String')))>0
    set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
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
        pcDatasets(pcDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
        for i=1:length(dirInfo)
            temp(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name)));
        end   
        data.mag = temp(:,:,1:floor(length(dirInfo)/2)); %%maybe change to velocity
        data.cd = temp(:,:,floor(length(dirInfo)/2)+1:end);
    else
        fid = fopen([pcDir '\pcvipr_header.txt'], 'r');
        dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
        fclose(fid);
        dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
        pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
        pcDatasets(pcDataIter).Info = pcviprHeader;
        resx = pcviprHeader.matrixx;  
        resy = pcviprHeader.matrixy;  
        nframes = pcviprHeader.frames;      
        MAG = load_dat(fullfile(pcDir,'MAG.dat'),[resx resy]);
        CD = load_dat(fullfile(pcDir,'CD.dat'),[resx resy]);
        VMEAN = load_dat(fullfile(pcDir,'comp_vd_3.dat'),[resx resy]);
        v = zeros(resx,resy,nframes);
        mag = zeros(resx,resy,nframes);
        cd = zeros(resx,resy,nframes);
        %%% Find the non-zero velocity direction
        for j = 1:nframes   
            v(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
            mag(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
            cd(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
        end 
        data.cd = cd;
        data.mag = mag;
        data.v = v;
        data.CD = CD;
        data.MAG = MAG;
        data.VMEAN = VMEAN;
    end
set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
pcDatasets(pcDataIter).Names = name;
pcDatasets(pcDataIter).Data = data;
pcDataIter = pcDataIter+1;
set(handles.DataName2DPC,'String','')
set(handles.Listbox2DPC,'String',{pcDatasets.Names});
end 

% --- Executes during object creation, after setting all properties.
function ButtonLoad2DPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonLoad2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Header2DPCListing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header2DPCListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in Listbox2DPC.
function Listbox2DPC_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox2DPC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox2DPC

% --- Executes during object creation, after setting all properties.
function Listbox2DPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Remove2DPC.
function Remove2DPC_Callback(hObject, eventdata, handles)
% hObject    handle to Remove2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pcDatasets pcDataIter
index = get(handles.Listbox2DPC, 'value');
pcDatasets(index) = [];
set(handles.Listbox2DPC,'Value',1);
set(handles.Listbox2DPC,'String',{pcDatasets.Names});
pcDataIter = pcDataIter-1;

% --- Executes during object creation, after setting all properties.
function Remove2DPC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Remove2DPC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



%%%%%%%%%%%%%% 2DMAG PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function Header2DMAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Load2DMAGText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Load2DMAGText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function DataName2DMAG_Callback(hObject, eventdata, handles)
% hObject    handle to DataName2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataName2DMAG as text
%        str2double(get(hObject,'String')) returns contents of DataName2DMAG as a double

% --- Executes during object creation, after setting all properties.
function DataName2DMAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataName2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global magDatasets magDataIter preloaded
if isempty(preloaded)
    magDatasets.Names = '';
    magDatasets.Data = 0;
    magDataIter = 1;   
end 

% --- Executes on button press in ButtonLoad2DMAG.
function ButtonLoad2DMAG_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoad2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global magDatasets magDataIter
set(handles.ErrorMessageBar,'String','')

if isempty(get(handles.DataName2DMAG,'String'))
    disp('Name is required to load data!');
    disp('Type a name into the text box and press the Load Data button.');
    set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
else 
    name = get(handles.DataName2DMAG,'String');
    [magFile, magDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE 2D magnitude file in the dataset');
    [~,~,extension] = fileparts(magFile);
    dirInfo = dir(fullfile(magDir,['*' extension]));
    if isequal(extension,'.dcm')
        magDatasets(magDataIter).Info = dicominfo(fullfile(magDir,dirInfo(1).name));
        for i=1:length(dirInfo)
            data(:,:,i) = single(dicomread(fullfile(magDir,dirInfo(i).name)));
        end        
    else
        temp = load(fullfile(magDir,dirInfo(1).name));
        data = single(struct2array(temp));
    end
if sum(strcmp({magDatasets.Names},name))>0
    set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
    set(handles.DataName2DMAG,'String','')
else 
    set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
    magDatasets(magDataIter).Names = name;
    magDatasets(magDataIter).Data = data;
    magDataIter = magDataIter+1;
    set(handles.DataName2DMAG,'String','')
    set(handles.Listbox2DMAG,'String',{magDatasets.Names});
end 
end 

% --- Executes during object creation, after setting all properties.
function ButtonLoad2DMAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonLoad2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Header2DMAGListing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header2DMAGListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in Listbox2DMAG.
function Listbox2DMAG_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox2DMAG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox2DMAG

% --- Executes during object creation, after setting all properties.
function Listbox2DMAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Remove2DMAG.
function Remove2DMAG_Callback(hObject, eventdata, handles)
% hObject    handle to Remove2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global magDatasets magDataIter
index = get(handles.Listbox2DMAG, 'value');
magDatasets(index) = [];
set(handles.Listbox2DMAG,'Value',1);
set(handles.Listbox2DMAG,'String',{magDatasets.Names});
magDataIter = magDataIter-1;

% --- Executes during object creation, after setting all properties.
function Remove2DMAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Remove2DMAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



%%%%%%%%%%%%%% 4D FLOW PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function Header4DFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Load4DFlowText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Load4DFlowText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function DataName4DFlow_Callback(hObject, eventdata, handles)
% hObject    handle to DataName4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataName4DFlow as text
%        str2double(get(hObject,'String')) returns contents of DataName4DFlow as a double

% --- Executes during object creation, after setting all properties.
function DataName4DFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataName4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global pcviprDatasets pcviprDataIter preloaded
if isempty(preloaded)
    pcviprDatasets.Names = '';
    pcviprDatasets.Data = 0;
    pcviprDataIter = 1;
end

% --- Executes on button press in ButtonLoad4DFlow.
function ButtonLoad4DFlow_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoad4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global pcviprDatasets pcviprDataIter
set(handles.ErrorMessageBar,'String','')

if isempty(get(handles.DataName4DFlow,'String'))
    disp('Name is required to load data!');
    disp('Type a name into the text box and press the Load Data button.');
    set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
else 
    name = get(handles.DataName4DFlow,'String');
    [pcviprFile, pcviprDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
   '*.dcm',  'DICOM files (*.dcm)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, 'Select ONE 4D flow file in the dataset');
    [~,~,extension] = fileparts(pcviprFile);
    dirInfo = dir(fullfile(pcviprDir,['*' extension]));
    if isequal(extension,'.dcm')
        pcviprDatasets(pcviprDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
        for i=1:length(dirInfo)
            data(:,:,i) = single(dicomread(fullfile(pcviprDir,dirInfo(i).name)));
        end        
    else
        temp = load(fullfile(pcviprDir,dirInfo(1).name));
        data = single(struct2array(temp));
    end
if sum(strcmp({pcviprDatasets.Names},name))>0
    set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
    set(handles.DataName4DFlow,'String','')
else 
    set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
    pcviprDatasets(pcviprDataIter).Names = name;
    pcviprDatasets(pcviprDataIter).Data = data;
    pcviprDataIter = pcviprDataIter+1;
    set(handles.DataName4DFlow,'String','')
    set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});
end 
end 

% --- Executes during object creation, after setting all properties.
function ButtonLoad4DFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonLoad4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Header4DFlowListing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Header4DFlowListing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in Listbox4DFlow.
function Listbox4DFlow_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Listbox4DFlow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox4DFlow

% --- Executes during object creation, after setting all properties.
function Listbox4DFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Remove4DFlow.
function Remove4DFlow_Callback(hObject, eventdata, handles)
% hObject    handle to Remove4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pcviprDatasets pcviprDataIter
index = get(handles.Listbox4DFlow, 'value');
pcviprDatasets(index) = [];
set(handles.Listbox4DFlow,'Value',1);
set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});
pcviprDataIter = pcviprDataIter-1;

% --- Executes during object creation, after setting all properties.
function Remove4DFlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Remove4DFlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



%%%%%%%%%%%%%% POST LOADING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in SaveCaseButton.
function CaseName_Callback(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CaseName as text
%        str2double(get(hObject,'String')) returns contents of CaseName as a double

% --- Executes during object creation, after setting all properties.
function CaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SaveCaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global anatDatasets pcDatasets magDatasets pcviprDatasets
if isempty(get(handles.CaseName,'String'))
    disp('Name is required to save data!');
    disp('Type a name into the text box and press the Save Case File button.');
    set(handles.ErrorMessageBar,'String','Name is required to save data! Type a name into the text box and press the Save Case File button.')
else
    disp('Select the directory where you will save the Case File.')
    set(handles.ErrorMessageBar, 'String', 'Select the directory where you will save the Case File.');
    casefileDir = uigetdir('C:\','Select the directory where you will save the Case File.');
    casename = get(handles.CaseName,'String');
    date = datestr(now);
    chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)];
    filename = [casename '_pwvcasefile_' chopDate '.mat'];
    set(handles.ErrorMessageBar,'String',['Saving Data as "' filename '"']);
    filenameWithPath = [casefileDir '\' filename];
    save(filenameWithPath,'anatDatasets','pcDatasets','magDatasets','pcviprDatasets');
    set(handles.ErrorMessageBar,'String','Data saved successfully');
end 


function ErrorMessageBar_Callback(hObject, eventdata, handles)
% hObject    handle to ErrorMessageBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ErrorMessageBar as text
%        str2double(get(hObject,'String')) returns contents of ErrorMessageBar as a double

% --- Executes during object creation, after setting all properties.
function ErrorMessageBar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ErrorMessageBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in ButtonCompleteLoading.
function ButtonCompleteLoading_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonCompleteLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global anatDataIter pcDataIter
global anatDatasets magDatasets pcDatasets
if anatDataIter==1
    disp('Anatomical datasets required to complete loading process.');
    set(handles.ErrorMessageBar,'String','Anatomical datasets required to complete loading process.');
elseif pcDataIter==1
    disp('2D PC datasets required to complete loading process.');
    set(handles.ErrorMessageBar,'String','2D PC datasets required to complete loading process.');
else
    disp('Loading Process Complete...');
    set(handles.ErrorMessageBar,'String','Loading Process Complete...');
    Anatomical_Images(anatDatasets);
    AnalyzePWV(pcDatasets,magDatasets)
end 

% --- Executes during object creation, after setting all properties.
function ButtonCompleteLoading_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonCompleteLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%%%%% UNUSED %%%%%%
% --- Executes during object creation, after setting all properties.
function MainPanel_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function MainPanel_DeleteFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function HeaderMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HeaderMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

