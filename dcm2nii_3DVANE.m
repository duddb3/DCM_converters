function dcm2nii_3DVANE(indir,outdir,rootname)
    % A custom DICOM to NIFTI converter to handle 3D VANE (Philips) images.
    % The code will create a NIFTI file for each image type that the scan
    % produce.
    %
    % Usage:
    %   dcm2nii_3DVANE(indir,outdir,rootname)
    %       indir: the file path containing all the DICOM files for the 3D
    %       VANE acquisition.
    %       outdir: the file path where the output images are written. File
    %       name are {rootname}_*.nii where * is F, W, IP, and OP for Fat,
    %       Water, In-Phase, and Out-of-Phase, respectively.

    if ~exist('rootname','var')
        [~,rootname] = fileparts(indir);
    elseif ~strcmp(rootname(1),'_')
        rootname = ['_' rootname];
    end

    % Get list of dicom files
    dcms = dir(fullfile(indir,'*.dcm'));
    % Parse the header for position information
    heads = dicominfo(fullfile(dcms(1).folder,dcms(1).name));
    if isscalar(dcms) && isfield(heads,'PerFrameFunctionalGroupsSequence')
        % enhanced DICOM
        A = squeeze(dicomread(heads));
        try
            % get the private tag info from the nested objects
            tags = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).Private_2005_140f.Item_1.Private_2005_1011',1:size(A,3),'Uni',0);
            imtype = categorical(cellfun(@(f) strtrim(char(f)),tags,'Uni',0));
        catch
            fprintf(2,'You need to code this Jon\n')
            return
        end
        types = unique(imtype);
        for n=1:length(types)
            idx = find(imtype==types(n));


            % populate metadata table with requisite information
            theads = table();
            % Rescale slope and intercept
            theads.RescaleSlope = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PixelValueTransformationSequence.Item_1.RescaleSlope,idx)';
            theads.RescaleIntercept = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PixelValueTransformationSequence.Item_1.RescaleIntercept,idx)';
            % Window level
            theads.WindowCenter = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).Private_2005_140f.Item_1.WindowCenter,idx)';
            theads.WindowWidth = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).Private_2005_140f.Item_1.WindowWidth,idx)';
            % Pixel info
            theads.PixelSpacing = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PixelMeasuresSequence.Item_1.PixelSpacing,idx,'Uni',0)';
            theads.SpacingBetweenSlices = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PixelMeasuresSequence.Item_1.SpacingBetweenSlices,idx)';
            % Orientation & Position
            theads.ImageOrientationPatient = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PlaneOrientationSequence.Item_1.ImageOrientationPatient,idx,'Uni',0)';
            theads.ImagePositionPatient = arrayfun(@(f) heads.PerFrameFunctionalGroupsSequence.(sprintf('Item_%i',f)).PlanePositionSequence.Item_1.ImagePositionPatient,idx,'Uni',0)';

            % Get the image in slice-ascending order
            [~,order] = sort(cellfun(@(f) f(3),theads.ImagePositionPatient));
            theads = theads(order,:);
            I = A(:,:,idx(order));

            % Save the nifti file
            fname = fullfile(outdir,[rootname '_' char(types(n)) '.nii']);
            if ~isfolder(outdir)
                mkdir(outdir)
            end
            writenii(I,fname,theads);
        end
    else
        parfor n=2:length(dcms)
            heads(n) = dicominfo(fullfile(dcms(n).folder,dcms(n).name));
        end
        heads = struct2table(heads);
        % Get Philips private head info specifying which image type each file
        % belongs to (i.e., W, F, IP, OP)
        try
            imtype = categorical(heads.Private_2005_1011);
            if ~all(ismember(imtype,{'W','F','IP','OP'}))
                fprintf(2,'Unexpected image types in input directory. Aborting\n')
                return
            end
        catch
            try
                utypes = unique(heads.ImageType,'stable');
                wi = find(contains(utypes,'\W\W\'));
                fi = find(contains(utypes,'\F\F\'));
                ipi = find(contains(utypes,'\IP\IP\'));
                opi = find(contains(utypes,'\OP\OP\'));
                if length(utypes)~=4 ||  ~all(ismember([wi fi ipi opi],1:4))
                    fprintf(2,'Unexpected image types in input directory. Aborting\n')
                    return
                end
                imtype = categorical(heads.ImageType,...
                    utypes([wi fi ipi opi]),...
                    {'W','F','IP','OP'});
            catch
                fprintf(2,'Input directory does not contain imaging data. Aborting\n')
                return
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
            tmpI = dicomread(theads.Filename{1});
            switch class(tmpI)
                case 'int16'
                    I = int16(zeros(theads.Rows(1),theads.Columns(1),height(theads)));
                case 'uint16'
                    I = uint16(zeros(theads.Rows(1),theads.Columns(1),height(theads)));
                otherwise
                    fprintf(sprintf('Image type was %s (expected int16 or uint16), skipping\n',class(tmpI)))
                    continue
            end
            parfor s=1:height(theads)
                % % Read the image
                I(:,:,s) = dicomread(theads.Filename{end-s+1});
            end

            % Save the nifti file
            fname = fullfile(outdir,['mDIXON_' char(types(n)) suffix '.nii']);
            if ~isfolder(outdir)
                mkdir(outdir)
            end
            writenii(I,fname,theads);
        end
    end

    function writenii(I,fname,theads)
        % Create a nifti header based on dicom metadata
        info.Filename = fname;
        info.Datatype = class(I);
        info.AdditiveOffset = theads.RescaleIntercept(1);
        info.MultiplicativeScaling = theads.RescaleSlope(1);
        info.DisplayIntensityRange = round(mean(cat(2,...
            theads.WindowCenter-theads.WindowWidth./2,...
            theads.WindowCenter+theads.WindowWidth./2)))';
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
        niftiwrite(rot90(flip(I)),fname,info);

        % Get the quarternions, add to header
        [u,s,v] = svd(info.Transform.T(1:3,1:3)');
        R = u*v';
        q = dcm2quat(R);
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