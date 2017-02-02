%=========================================================================%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%=========================================================================%
% Copyright (C) 2013 - Cyclotron Research Centre
% University of Liege, Belgium

function hdr = eb_read_protocol(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:
% this script reads almost everything contained in the header (meas*.hdr)
% of raw data acquired under VB15/VB17, including ASCII parameters and
% XProtocol parameters of all kind. The result is stored in a structure.
% In general, most of the useful parameters can be found in
% hdr.asc and/or hdr.Meas.YAPS.
%
% USAGE:
% hdr = eb_read_protocol(filename)
% where filename is the header file name (format meas*.hdr) and hdr is a
% structure containing all the retrieved parameters.
%
% WARNING AND DISCLAIMER: 
% This software is for research use only. Do not use it for clinical or
% diagnostic purposes. 
%--------------------------------------------------------------------------
% History
% - first implemented in February 2013 for VB15/17 
% - modified to read VD13 headers from Prisma in October 2013
%       the structure repeats twice and only the second one contains
%       actually useful parameters (the first structure seems to contain
%       limits, parameter attributes and the like). The current code has
%       been slightly modified to take this into account, using cell arrays
%       to differenciate iterations. The ASCII portion of the hdr appears
%       actually 4 times and is also changed into a cell array to avoid
%       overwriting values. 
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENABLE_DEBUG = false;

fid = fopen(filename);
if fid == -1
    disp('File does not exist or cannot be opened');
    return;
end

disp(['Reading headers from file ' filename ]);

% initialise a few variables
hdr = [];
cdepth = 1;
cbranch{cdepth} = 'hdr';
cstruct = 'hdr';
xprot = {'Config','Dicom','Meas','Phoenix','Spice'};  
hdrtype = 'xprot';
new_field_name = '';
new_field_value = '';
new_field_type = '';

% starts reading...
%tic
cline = fgets(fid); %returns the next line of the file, if end-of-file, returns -1
oline = '';
clinenum = 1;
expecting_value = 1;
while (cline ~= -1 & cdepth>0)
    if ENABLE_DEBUG
        disp(['Line #' num2str(clinenum) ' - depth = ' num2str(cdepth) ' - cbranch = ' cbranch{cdepth} ' - ' cstruct ' - [' cline(1:end-1) ']']);
    end
    
    idx_asconv_begin = findstr(cline,'### ASCCONV BEGIN'); % for Prisma data
    idx_asconv_end = findstr(cline,'### ASCCONV END ###');
    
    if ~isempty(idx_asconv_begin)
        hdrtype = 'ascii';
        disp('Reading ASCII headers');

        new_field_name = 'asc';
        if isfield(hdr, new_field_name)
            rep = eval(['length(hdr.' new_field_name ')'])+1;
        else
            rep = 1;
        end
        new_field_name = [new_field_name '{' num2str(rep) '}'];
        
        oline = cline;
        cline = fgets(fid); % go to next line where ascii headers start
    end
    if ~isempty(idx_asconv_end)
        hdrtype = 'xprot';
        %break;
    end

    switch (hdrtype)
        case 'ascii'
            %tmp = textscan(cline,'%s = %s'); % problematic if any space in the value... 
            %hdrnam = tmp{1}{1};
            %hdrval = tmp{2}{1};
            
            [hdrnam, hdrval] = strtok(cline,'=');
            hdrval = strtok(hdrval,'='); 
            hdrval = hdrval(1:end-1); % to remove the end-of-line character...
            hdrval = strtrim(hdrval); % to remove leading and trailing white space from string 
            
            % first process hdrnam
            % skip if contains '_' characters (usually parameter attribute
            % in Prisma data) 
            if isempty(strfind(hdrnam,'_'))
                % convert indexes if any: C++ [0],[1],... -> matlab (1),(2),...
                idx = [strfind(hdrnam,'[');strfind(hdrnam,']')];
                if ~isempty(idx)
                    % remplace [] by ()
                    hdrnam(idx(1,:)) = '('; hdrnam(idx(2,:)) = ')';
                    % increment indexes
                    for in = 1:(size(idx,2))
                        oldidx = hdrnam(idx(1,in)+1:idx(2,in)-1);
                        newidx = num2str(str2num(oldidx)+1);
                        hdrnam = [hdrnam(1:idx(1,in)) newidx hdrnam(idx(2,in):end)];
                        if (length(newidx)>length(oldidx) && in<size(idx,2))
                            idx(:,in+1:end) = idx(:,in+1:end)+1;
                        end
                    end
                end
                % now process hdrval
                idx = strfind(hdrval, '"');
                if ~isempty(idx)
                    eval(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval(idx(1)+1:idx(2)-1) '''};']);
                elseif ~isempty(strfind(hdrval, 'x'))
                    eval(['hdr.' new_field_name '.' hdrnam ' = {''' hdrval '''};']);
                else
%                     ['hdr.asc.' hdrnam ' = str2num(hdrval);']
%                     hdrval
                    eval(['hdr.' new_field_name '.' hdrnam ' = str2num(hdrval);']);
                end
            end

        case 'xprot'
            idx = find(cline=='<'|cline=='{'|cline=='}');
            if ~isempty(idx)
                i = 1;
                while (i<length(idx)+1)
                    switch cline(idx(i))
                        case '<'
                            % truc = '     <ParamBool."ucDixonSaveOriginal"> '
                            % kk = textscan(truc,'%s %q','delimiter','<>."','MultipleDelimsAsOne',1)
                            % kk{1} = 'ParamBool' & kk{2} = 'ucDixonSaveOriginal'
                            old_field_type = new_field_type;
                            new_field_name = cline(idx(i):end);
                            new_field_name = new_field_name(find(new_field_name=='<')+1:find(new_field_name=='>')-1);
                            [new_field_type, new_field_name] = strtok(new_field_name,'.');
                            idxq = find(new_field_name == '"');
                            if ~isempty(idxq);new_field_name = new_field_name(idxq(1)+1:idxq(2)-1);end
                            if isempty(new_field_name)
                                switch new_field_type
                                    case 'XProtocol'
                                        % Scan the current and previous
                                        % lines to determine which
                                        % XProtocol it is:
                                        new_field_name = 'unknown';
                                        for cxprot=1:length(xprot)
                                            if (~isempty(findstr(cline, xprot{cxprot})) | ~isempty(findstr(oline, xprot{cxprot})) )
                                                new_field_name = xprot{cxprot};
                                                disp(['Reading XProtocol ' xprot{cxprot}]);
                                            end
                                        end
                                        if isfield(hdr, new_field_name)
                                            rep = eval(['length(hdr.' new_field_name ')'])+1;
                                        else
                                            rep = 1;
                                        end
                                        new_field_name = [new_field_name '{' num2str(rep) '}'];
                                    case 'ParamMap'
                                        cbranch{cdepth+1} = '';
                                    case {'Precision','MaxSize','MinSize','Default'}
                                        new_field_type = old_field_type;
                                    otherwise
                                end
                            end
                            if ~isempty(new_field_name)
                                % some field names start with a number 
                                % -> add 'x' before in that case:
                                if ~isletter(new_field_name(1))
                                    new_field_name = ['x' new_field_name];
                                end
                                cstruct = [];
                                cbranch{cdepth+1} = ['.' new_field_name];
                                for cd=1:cdepth+1;cstruct=[cstruct cbranch{cd}];end
                            end

                        case '{'
                            cdepth = cdepth+1;
                            expecting_value = true;
                            new_field_value = cline(idx(i):end);
                            % NB: we ONLY deal with values given between
                            % a single pair of curly brackets. If
                            % <Precision> given, it must be removed before
                            % further processing...
                            
                            % retrieve what is between the first and last
                            % curly brackets if a pair is present (returns
                            % empty new_field_value otherwise):
                            new_field_value = new_field_value(find(new_field_value=='{',1,'first')+1:find(new_field_value=='}',1,'last')-1);
                            
                            % if not empty...
                            if ~isempty(new_field_value)
                                str2eval = retrieve_value(new_field_value, new_field_type, cstruct);
                                eval(str2eval);
                                i = length(idx);
                                cdepth = cdepth-1;
                            end

                        case '}'
                            cdepth = cdepth-1;
                            expecting_value = false;

                        otherwise
                    end
                    i = i+1;
                end
            else
                if expecting_value
                    new_field_value = cline;
                    str2eval = retrieve_value(new_field_value, new_field_type, cstruct);
                    eval(str2eval);
                end
            end
    end
    oline = cline;
    cline = fgets(fid);
    clinenum = clinenum+1;
end
%toc

function str2eval = retrieve_value(new_field_value, new_field_type, cstruct)
% subfunction returning a string to be evaluated in order to assign a value
% in the hdr structure (made substructure since it can be called from 2
% distinct places in the code)

str2eval = '';

% check whether there are no other pair of {} within the string
if isempty(find(new_field_value=='{'))
    % check whether <Precision> is given
    if findstr(new_field_value,'<Precision>')
        tmp = textscan(new_field_value,'<Precision> %s %s');
        if ~isempty(tmp{2})
            new_field_value = tmp{2}{1};
        else
            new_field_value = '  ';
        end
    end
    % check whether other params (<MinSize>, <MaxSize>, etc) are given
    % (skip in this case)... 
    if isempty(strfind(new_field_value,'<'))
        if isempty(find(~isspace(new_field_value)))
            str2eval = [cstruct ' = [];'];
        else
            switch new_field_type
                case 'ParamString'
                    idxq = find(new_field_value=='"');
                    if ~isempty(idxq)
                        if length(idxq)==1;idxq(2)=length(new_field_value);end
                        new_field_value = new_field_value(idxq(1)+1:idxq(2)-1);
                        str2eval = [cstruct ' = ''' new_field_value ''';'];
                    else
                        str2eval = [cstruct ' = '''';'];
                    end
                case 'ParamLong'
                    str2eval = [cstruct ' = [' new_field_value '];'];
                case 'ParamDouble'
                    str2eval = [cstruct ' = [' new_field_value '];'];
                case 'ParamBool'
                    idxq = find(new_field_value=='"');
                    new_field_value = new_field_value(idxq(1)+1:idxq(2)-1);
                    str2eval = [cstruct ' = ' new_field_value ';'];
                otherwise
            end
        end
    end
    % For Prisma...
    idx = findstr(str2eval, '@');
    for i=length(idx):-1:1
        str2eval = [str2eval(1:idx(i)-1) 'AT' str2eval(idx(i)+1:end)];
    end
end
