function dcm2nii_3DVANE(indir,outdir)
    % A custom DICOM to NIFTI converter to handle 3D VANE (Philips) images.
    % The code will create a NIFTI file for each image type that the scan
    % produce.
    %
    % Usage:
    %   dcm2nii_3DVANE(indir,outdir)
    %       indir: the file path containing all the DICOM files for the 3D
    %       VANE acquisition.
    %       outdir: the file path where the output images are written. File
    %       name are 3D_VANE_*.nii where * is F, W, IP, and OP for Fat,
    %       Water, In-Phase, and Out-of-Phase, respectively.

    % Make sure target directory exists
    if ~isfolder(outdir)
        mkdir(outdir)
    end

    % Get list of dicom files
    dcms = dir(fullfile(indir,'*.dcm'));
    % Parse the header for position information
    heads = dicominfo(fullfile(dcms(1).folder,dcms(1).name));
    parfor n=2:length(dcms)
        heads(n) = dicominfo(fullfile(dcms(n).folder,dcms(n).name));
    end
    heads = struct2table(heads);

    % Get Philips private head info specifying which image type each file
    % belongs to (i.e., W, F, IP, OP)
    try
        imtype = categorical(heads.Private_2005_1011);
    catch
        try
        temprep = [heads.Private_2005_1011{:}]';
        imtype = cell(length(temprep),1);
        for n=1:length(temprep)
            imtype{n} = strtrim(char(temprep(n,:)));
        end
        imtype = categorical(imtype);
        catch
            cellfun(@(f) strrep(f,'\JP2K LOSSY',''),heads.ImageType,'UniformOutput',false)
            imtype = categorical(heads.ImageType,...
                {'DERIVED\PRIMARY\F\F\DERIVED','DERIVED\PRIMARY\W\W\DERIVED','DERIVED\PRIMARY\IP\IP\DERIVED','DERIVED\PRIMARY\OP\OP\DERIVED'},...
                {'F','W','IP','OP'});
        end
    end
    types = unique(imtype);

    % For each image type
    for n=1:length(types)
        % get headers belonging to that image type
        theads = heads(imtype==types(n),:);
        % find the unique slice locations and their order (ascending)
        [~,ia,~] = unique(theads.SliceLocation);
        theads = theads(ia,:);
        % Instantiate image array
        I = zeros(theads.Rows(1),theads.Columns(1),height(theads));
        parfor s=1:height(theads)
            % Read the image, applying the rescale parameters
            I(:,:,s) = double(dicomread(theads.Filename{s})).*...
                theads.RescaleSlope(s)+theads.RescaleIntercept(s);
        end
        % Flip along z-axis
        I = flip(I,3);
        % Create a nifti header based on dicom metadata
        fname = fullfile(outdir,['3D_VANE_' char(types(n)) '.nii']);
        mn = min(I(:));
        mx = max(I(:));
        I = gray2ind(mat2gray(I,[mn mx]),4^8);
        info.AdditiveOffset = mn;
        info.MultiplicativeScaling = mx/(4^8);
        info.DisplayIntensityRange = [mn mx];
        info.Datatype = 'uint16';
        info.Filename = fname;
        info.Filemoddate = char(datetime);
        info.Version = 'NIfTI1';
        info.Description = '';
        info.ImageSize = size(I);
        info.PixelDimensions = [theads.PixelSpacing{1}' theads.SpacingBetweenSlices(1)];
        info.BitsPerPixel = 16;
        info.SpaceUnits = 'Millimeter';
        info.TimeUnits = 'Second';
        info.TimeOffset = 0;
        info.SliceCode = 'Unknown';
        info.FrequencyDimension = 1;
        info.PhaseDimension = 2;
        info.SpatialDimension = 3;
        info.TransformName = 'Sform';
        info.Transform.Dimensionality = 3;
        info.Transform.T = [info.PixelDimensions 1].*[-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
        info.Transform.T(4,1:3) = theads.ImagePositionPatient{1};
        info.Qfactor = 1;
        % Save the nifti file
        niftiwrite(flip(rot90(I)),fname,info);
    end
end