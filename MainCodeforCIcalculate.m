clear
CropIntensity = zeros(4800,4800);
load('Threshold.mat')
InPath = '...\MODIS\MOD13Q1\';
tile = dir(fullfile(Path,'h*probability.tif'));
for HV_i = 1:length(tile)
    HV = tile(HV_i).name;
    HVstr = HV(1:6);
    CropMaskFile = strcat(Path,'\Crop\',lower(HVstr),'_probability.tif');
    CropMask = imread(CropMaskFile);
    [num_x,num_y] = find(CropMask >= 10);
    ProjFile =  strcat(Path,'250m_',lower(HVstr),'.tif');
    [~,~,R] = geotiffread(ProjFile);
    proj=geotiffinfo(ProjFile);
    for year = 2001:2019
        WorkSpace0 = strcat(InPath,'\',num2str(year-1));
        WorkSpace1 = strcat(InPath,'\',num2str(year));
        WorkSpace2 = strcat(InPath,'\',num2str(year+1));
        searchKeyword = strcat('\*\*',HV_group{HV_i},'*');
        folder0 = dir(fullfile(WorkSpace0,searchKeyword));
        folder1 = dir(fullfile(WorkSpace1,searchKeyword));
        folder2 = dir(fullfile(WorkSpace2,searchKeyword));
        folderList = [folder0(end-5:end,1);folder1;folder2(1:6,1)];
        EVI =  zeros(4800,4800,length(folderList));
        QA =  zeros(4800,4800,length(folderList));
        DQA = zeros(4800,4800,length(folderList));
        for f = 1:length(folderList)
            Filefolder = folderList(f).folder;
            Filename = folderList(f).name;
            FilenameSource = strcat(Filefolder,'\',Filename);
            EVI(:,:,f) = hdfread(FilenameSource,'250m 16 days EVI');
            QA(:,:,f) = hdfread(FilenameSource,'250m 16 days pixel reliability');
            DQA(:,:,f) = hdfread(FilenameSource,'250m 16 days VI Quality');
        end
        fprintf('Processing: done %s year perparision ...\r',num2str(year));
        tic
        bias = 100;
        ts = 5;
        for n = 1:length(num_x)
            i = num_x(n);
            j = num_y(n);
            ndvi = EVI(i,j,:);
            ndvi = ndvi(:)';
            qa = QA(i,j,:);
            qa = qa(:)';
            dqa = DQA(i,j,:);
            dqa = dqa(:)';
            CropIntensity(i,j) = polyfitCI(QAD,evi,qa,dqa,peak,bias,ts);
        end 
        T_time = toc;
        fprintf('Processing: year [%s/2019] Tile [%s] in %s time ...\r',num2str(year),HVstr,num2str(T_time));
        OutFilename = strcat(Path,'\CI_',lower(HVstr),'_',num2str(year),'.tif');
        geotiffwrite(OutFilename, CropIntensity, R, 'GeoKeyDirectoryTag', proj.GeoTIFFTags.GeoKeyDirectoryTag);
        clear CropIntensity
    end
end