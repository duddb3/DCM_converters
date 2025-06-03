function dcm2nii_mDIXON(indir,outdir)
    % A custom DICOM to NIFTI converter to handle mDIXON Quant 
    % (from Philips) images. The code will create a NIFTI file for each 
    % image type that the scan produces.
    %
    % Usage:
    %   dcm2nii_mDIXON(indir,outdir)
    %       indir: the file path containing all the DICOM files for the 
    %           mDIXON Quant acquisition.
    %       outdir: the file path where the output images are written. File
    %       name are mDIXON_*.nii where * is W, F, FF, and T2_STAR for
    %       Water, Fat, Fat Fraction, and T2*, respectively.

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
    % belongs to (i.e., W, F, FF, T2_STAR)
    try
        imtype = categorical(cellfun(@(f) strtrim(char(f')),heads.Private_2005_1011,'Uni',0));
    catch
        utypes = unique(heads.ImageType,'stable');
        wi = find(contains(utypes,'\W\W\'));
        fi = find(contains(utypes,'\F\F\'));
        ffi = find(contains(utypes,'\FF\FF\'));
        t2i = find(contains(utypes,'\T2_STAR'));
        imtype = categorical(heads.ImageType,...
            utypes([wi fi ffi t2i]),...
            {'W','F','FF','T2_STAR'});
    end
    types = unique(imtype);

    % For each image type
    for n=1:length(types)
        % get headers belonging to that image type
        theads = heads(imtype==types(n),:);
        % find the unique slice locations and their order (ascending)
        [~,ia,~] = unique(theads.SliceLocation,'sorted');
        theads = theads(ia,:);
        % Instantiate image array
        tmpI = dicomread(theads.Filename{1});
        switch class(tmpI)
            case 'int16'
                I = int16(zeros(theads.Rows(1),theads.Columns(1),height(theads)));
                info.Datatype = 'int16';
            case 'uint16'
                I = uint16(zeros(theads.Rows(1),theads.Columns(1),height(theads)));
                info.Datatype = 'uint16';
            otherwise
                continue
        end
        parfor s=1:height(theads)
            % % Read the image
            I(:,:,s) = dicomread(theads.Filename{end-s+1});
        end
        % Create a nifti header based on dicom metadata
        fname = fullfile(outdir,['mDIXON_' char(types(n)) '.nii']);
        info.Filename = fname;
        info.AdditiveOffset = theads.RescaleIntercept(1);
        info.MultiplicativeScaling = theads.RescaleSlope(1);
        info.DisplayIntensityRange = [
            theads.WindowCenter(1)-theads.WindowWidth(1)/2
            theads.WindowCenter(1)+theads.WindowWidth(1)/2];
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
        ImOrPat = theads.ImageOrientationPatient{1};
        RotMat = [ImOrPat(1:3) ImOrPat(4:6) cross(ImOrPat(1:3),ImOrPat(4:6))];
        info.Transform.T = info.PixelDimensions.*RotMat;
        info.Transform.T(4,4) = 1;
        info.Transform.T(4,1:3) = theads.ImagePositionPatient{1};
        info.Qfactor = 1;
        % Save the nifti file
        niftiwrite(rot90(flip(I)),fname,info);

        % Get the quarternions, add to header
        q = dcm2quat(info.Transform.T(1:3,1:3)');
        fid = fopen(fname,'r+');
        fseek(fid,252,'bof');
        fwrite(fid,1,'short');
        fseek(fid,256,'bof');
        fwrite(fid,[q(2:4)' info.Transform.T(4,1:3)],'float');
        fclose(fid);
        gzip(fname)
        delete(fname)


    end


    % Retrun quaternion abcd from normalized matrix R (3x3)
    function [q, proper] = dcm2quat(R)
        % Simplied from Quaternions by Przemyslaw Baranski 
        proper = sign(det(R));
        if proper<0
            R(:,3) = -R(:,3);
        end
        
        q = sqrt([1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1] * diag(R) + 1) / 2;
        if ~isreal(q(1))
            % if trace(R)+1<0, zero it
            q(1) = 0;
        end 
        [mx, ind] = max(q);
        mx = mx * 4;
        
        if ind == 1
            q(2) = (R(3,2) - R(2,3)) /mx;
            q(3) = (R(1,3) - R(3,1)) /mx;
            q(4) = (R(2,1) - R(1,2)) /mx;
        elseif ind ==  2
            q(1) = (R(3,2) - R(2,3)) /mx;
            q(3) = (R(1,2) + R(2,1)) /mx;
            q(4) = (R(3,1) + R(1,3)) /mx;
        elseif ind == 3
            q(1) = (R(1,3) - R(3,1)) /mx;
            q(2) = (R(1,2) + R(2,1)) /mx;
            q(4) = (R(2,3) + R(3,2)) /mx;
        elseif ind == 4
            q(1) = (R(2,1) - R(1,2)) /mx;
            q(2) = (R(3,1) + R(1,3)) /mx;
            q(3) = (R(2,3) + R(3,2)) /mx;
        end
        if q(1)<0
            q = -q;
        end
    end
end