classdef MHAReader < handle
    % MHAReader  A class for reading non-compressed Volume Sequence Image (.mha) files
    %
    % This MATLAB class is adapted from the provided C++ MHAReader code. 
    % It parses the header fields of an .mha file, then reads the raw 
    % (uncompressed) volume data. For best compatibility, ensure your .mha 
    % has 'ElementDataFile = LOCAL'.
    %
    % EXAMPLE:
    %    reader = MHAReader('path_to_your_volume.mha');
    %    success = reader.readVolumeImage();
    %    if success
    %        hdr = reader.getMHAHeader();
    %        vol = reader.getMHAVolume();
    %        % Use 'vol' as needed...
    %    end
    
    properties (Access = private)
        volumeimage   % Stores the raw binary data of the voxels (uint8 array)
        filename      % Full path to the .mha file
        header        % Struct storing .mha header fields
    end
    
    methods
        function obj = MHAReader(filename)
            % Constructor that initializes the object with the .mha filename.
            if ~isfile(filename)
                error('File does not exist or is not accessible: %s', filename);
            end
            
            obj.filename = filename;
            
            % Initialize the header struct fields akin to the C++ MHAHeader.
            % You can expand these if you discover your .mha file uses more fields.
            obj.header = struct( ...
                'ObjectType',                 "", ...
                'NDims',                      0, ...
                'BinaryData',                false, ...
                'BinaryDataByteOrderMSB',    false, ...
                'CompressedData',            false, ...
                'TransformMatrix',           [], ...
                'DimSize',                   [], ...
                'Offset',                    [], ...
                'CenterOfRotation',          [], ...
                'AnatomicalOrientation',     "", ...
                'ElementSpacing',            [], ...
                'ElementType',               "", ...
                'UltrasoundImageOrientation', "", ...
                'UltrasoundImageType',       "", ...
                'ElementDataFile',           "" ...
            );
            
            obj.volumeimage = [];
        end
        
        function success = readVolumeImage(obj)
            % READVOLUMEIMAGE  Reads the entire .mha file, parses the header,
            %                  and extracts the raw binary volume data.
            %
            % RETURNS:
            %   success (boolean) indicating success/failure.
            
            fid = fopen(obj.filename, 'rb');
            if fid == -1
                warning('Unable to open file: %s', obj.filename);
                success = false;
                return;
            end
            
            % Read the entire file into a uint8 array (binary-safe).
            rawData = fread(fid, Inf, '*uint8');
            fclose(fid);
            
            % Convert that to a char array for searching/parsing text.
            % (In modern MATLAB, char(0) is allowed—so embedded null bytes 
            %  won't necessarily terminate the string, but can complicate searching.)
            fileContents = char(rawData.');
            
            % Look for the line that indicates the start of LOCAL data.
            % Avoid including newlines, as different OS may use different EOL.
            separator = 'ElementDataFile = LOCAL';
            idx = strfind(fileContents, separator);
            if isempty(idx)
                % If we didn't find the "ElementDataFile = LOCAL" marker,
                % we can't properly split text from binary as in the original logic.
                warning('Could not find "%s" in file. Aborting.', separator);
                success = false;
                return;
            end
            
            % The header text ends right after the 'ElementDataFile = LOCAL' substring.
            % We also skip over any trailing spaces/newlines to find the first 
            % binary byte:
            binaryStartPos = idx(1) + length(separator);
            
            % Possibly skip over trailing space, \r, or \n:
            while binaryStartPos <= length(fileContents) ...
                    && ismember(fileContents(binaryStartPos), [9 10 13 32]) % tab/nl/cr/space
                binaryStartPos = binaryStartPos + 1;
            end
            
            % The text portion is everything from the beginning up to 
            % (but not including) the binary payload start.
            normalText = fileContents(1 : binaryStartPos - 1);
            
            % The remaining bytes in rawData from 'binaryStartPos' onward:
            % Remember, 'binaryStartPos' is an index in fileContents (char array).
            % We need the same index in the rawData (uint8). 
            % Because we built fileContents = char(rawData.').  
            % => If everything is 1:1, that means rawData(binaryStartPos) 
            %    corresponds to fileContents(binaryStartPos). 
            % That should be correct as long as we haven't lost any data 
            % (no early null termination).
            binaryData = rawData(binaryStartPos : end);
            
            % Parse the header lines (everything before the binary).
            successReadHeader = obj.readHeader(normalText);
            if ~successReadHeader
                warning('Error parsing header text section.');
                success = false;
                return;
            end
            
            % Transfer binary data to volumeimage. (Uncompressed only)
            try
                if ~obj.header.CompressedData
                    if ~strcmpi(obj.header.ElementDataFile, 'LOCAL')
                        warning(['ElementDataFile is not "LOCAL". ' ...
                                 'External file reading not implemented.']);
                        success = false;
                        return;
                    end
                    obj.volumeimage = uint8(binaryData);
                else
                    % The original C++ code complains about compressed data. 
                    % We do the same here.
                    warning('Compressed data is not handled in this implementation.');
                    success = false;
                    return;
                end
            catch ME
                warning(sprintf('Exception while reading volume data: %s', ME.message));
                success = false;
                return;
            end
            
            success = true;
        end
        
        function hdr = getMHAHeader(obj)
            % GETMHAHEADER  Returns the parsed .mha header as a struct.
            hdr = obj.header;
        end
        
        function volume = getMHAVolume(obj)
            % GETMHAVOLUME  Returns the volume data (uint8) 
            %               that was read from the file.
            volume = obj.volumeimage;
        end
    end
    
    methods (Access = private)
        function success = readHeader(obj, textData)
            % READHEADER  Parses key-value lines from the .mha header text 
            %             and populates obj.header accordingly.
            %
            % RETURNS:
            %   success (bool) indicating if it parsed without issues.
            
            lines = strsplit(textData, {'\r', '\n'});
            
            for iLine = 1:numel(lines)
                line = strtrim(lines{iLine});
                if isempty(line)
                    continue;  % skip empty lines
                end
                
                % Find the '=' delimiter.
                equalPos = strfind(line, '=');
                if isempty(equalPos)
                    continue;  % skip lines that don't have a key=value pair
                end
                
                key   = strtrim(line(1:equalPos(1)-1));
                value = strtrim(line(equalPos(1)+1:end));
                
                switch key
                    case 'ObjectType'
                        obj.header.ObjectType = value;
                        
                    case 'NDims'
                        obj.header.NDims = str2double(value);
                        
                    case 'BinaryData'
                        obj.header.BinaryData = strcmpi(value, 'True');
                        
                    case 'BinaryDataByteOrderMSB'
                        obj.header.BinaryDataByteOrderMSB = strcmpi(value, 'True');
                        
                    case 'CompressedData'
                        obj.header.CompressedData = strcmpi(value, 'True');
                        
                    case 'TransformMatrix'
                        vals = str2double(strsplit(value));
                        obj.header.TransformMatrix = vals;
                        
                    case 'DimSize'
                        vals = str2double(strsplit(value));
                        obj.header.DimSize = vals;
                        
                    case 'Offset'
                        vals = str2double(strsplit(value));
                        obj.header.Offset = vals;
                        
                    case 'CenterOfRotation'
                        vals = str2double(strsplit(value));
                        obj.header.CenterOfRotation = vals;
                        
                    case 'AnatomicalOrientation'
                        obj.header.AnatomicalOrientation = value;
                        
                    case 'ElementSpacing'
                        vals = str2double(strsplit(value));
                        obj.header.ElementSpacing = vals;
                        
                    case 'ElementType'
                        obj.header.ElementType = value;
                        
                    case 'UltrasoundImageOrientation'
                        obj.header.UltrasoundImageOrientation = value;
                        
                    case 'UltrasoundImageType'
                        obj.header.UltrasoundImageType = value;
                        
                    case 'ElementDataFile'
                        obj.header.ElementDataFile = value;
                        
                    otherwise
                        % Unrecognized header key, ignoring...
                end
            end
            
            % If we get this far, we consider it a success. 
            % You might add more checks here if you want to verify 
            % certain fields are actually populated.
            success = true;
        end
    end
end
