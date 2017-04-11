
%% project directory
% PROJpath = '/home/ahn/metastasis';
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% location of excel file
LOCpath  = fullfile(PROJpath,'coordinates');

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170106.xlsx');
[a,b,xlsData] = xlsread(fn_list);
subjlist = cell2mat(xlsData(2:end,3));




%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------
outpath = uigetdir([],'Select directory to save');
fn_out = fullfile(outpath,'count_metastasis.csv');
fid = fopen(fn_out,'w+');
fprintf(fid,'subjname,nBM\n');
%%
nsubj = length(subjlist);
for c=1:nsubj,
    
    % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    
    
    % Load center coordinates of metastasis
    %----------------------------------------------------------------------
    fn_xls = fullfile(LOCpath,sprintf('%s.xls',subjname));
    try
        [a,b,xlsData] = xlsread(fn_xls);
        Vxyz = round(cell2mat(xlsData(2:end,3:5)));
        nBM = size(Vxyz,1);
    catch
        fprintf('no data for %s\n',subjname);
        nBM = nan;
    end
    fprintf(fid,'%s,%d\n',subjname,nBM);
end
fclose(fid);