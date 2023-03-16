function vtkwrite_pathlines( filename,dataType,varargin )
% VTKWRITE Writes 3D Matlab pathlines into VTK file format.

%  vtkwrite(filename, 'polydata', 'lines', x, y, z) exports a 3D line where
%  x,y,z are coordinates of the points that make the line. x, y, z are
%  vectors containing the coordinates of points of the line, where point(n)
%  is specified by x(n), y(n) and z(n).
%  
%  Version 2.3
%  Copyright, Chaoyuan Yeh, 2016
%  Codes are modified by Xinan Chen on 05/29/20

if strcmpi(filename,'execute'), filename = 'matlab_export.vtk'; end
fid = fopen(filename, 'w'); 
% VTK files contain five major parts
% 1. VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 2.0\n');
% 2. Title
fprintf(fid, 'VTK from Matlab\n');


binaryflag = any(strcmpi(varargin, 'BINARY'));
if any(strcmpi(varargin, 'PRECISION'))
    precision = num2str(varargin{find(strcmpi(vin, 'PRECISION'))+1});
else
    precision = '2';
end

switch upper(dataType)        
    case {'STRUCTURED_POINTS','STRUCTURED_GRID','UNSTRUCTURED_GRID'}
        warning('Only POLYDATA in vtkwrite_pathlines!')
        return
    case 'POLYDATA'

        fprintf(fid, 'ASCII\n');
        if numel(varargin)<2, error('Not enough input arguments'); end
        fprintf(fid, 'DATASET POLYDATA\n');
        
        countpts = zeros(1,size(varargin{2},2));
        x = []; y = x; z = x;
        for plN = 1:size(varargin{2},2)
            xtmp = varargin{2}{plN}(:,1);%varargin{2}{plN}(:,2);
            ytmp = varargin{2}{plN}(:,2);%varargin{2}{plN}(:,1);
            ztmp = varargin{2}{plN}(:,3);
            if sum(size(xtmp)==size(ytmp) & size(ytmp)==size(ztmp))~= length(size(xtmp))
                error('Input dimesions do not match')
            end
            x = [x;xtmp]; y = [y;ytmp]; z = [z;ztmp];
            countpts(plN) = numel(xtmp);            
        end
        n_elements = numel(x);            
        if mod(n_elements,3)==1
            x(n_elements+1:n_elements+2,1)=[0;0];
            y(n_elements+1:n_elements+2,1)=[0;0];
            z(n_elements+1:n_elements+2,1)=[0;0];
        elseif mod(n_elements,3)==2
            x(n_elements+1,1)=0;
            y(n_elements+1,1)=0;
            z(n_elements+1,1)=0;
        end
        output = [x(1:3:end-2), y(1:3:end-2), z(1:3:end-2),...
                      x(2:3:end-1), y(2:3:end-1), z(2:3:end-1),...
                      x(3:3:end), y(3:3:end), z(3:3:end)]';
        %nbpoint = sum(countpts);
        nbpoint = numel(x);
        spec = [repmat(['%0.', precision, 'f '], 1, 9), '\n'];
        fprintf(fid, ['POINTS ' num2str(nbpoint) ' float\n']);
        fprintf(fid, spec, output);
        
        switch upper(varargin{1})
            case 'LINES'
                %nbLine = 2*(n_elements-1 - size(varargin{2},2)+1);
                nbLine = 2*(n_elements-1);
                % create conns to be deleted
                countsum = cumsum(countpts); countsum = countsum(1:end-1);
                tmp1 = 1:nbLine/2;
                [~,ia,~] = intersect(tmp1,countsum);
                tmp1(ia)=[];

                conn21 = tmp1;%1:nbLine/2;
                conn11 = conn21-1;
                conn12 = tmp1;%1:nbLine/2;
                conn22 = conn12-1;
                conn1 = [conn11,conn12];
                conn2 = [conn21,conn22];
                %s = [conn1',conn2'];
                fprintf(fid,'\nLINES %d %d\n',length(conn1),3*length(conn1));
                fprintf(fid,'2 %d %d\n',[conn1;conn2]);
                
                %% write start and end points/point values
                fprintf(fid,'\nPOINT_DATA %d\n',nbpoint);
                fprintf(fid,'SCALARS PathPoint float\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                ptsv = zeros(1,nbpoint);
                
                
                % ==== type 1 ====: set start as -1, end as 1 and intermediate as 0
                %{
                countsum2 = cumsum(countpts);
                ptsv([1,countsum+1]) = -1; % start points
                ptsv(countsum2) = 1; % end points
                %}
                % ==== type 2 ====: set start as -1, end as 1 and intermediate as
                % equally spaced
                rr = [1,countsum+1];
                for plN = 1:size(varargin{2},2)
                    ptsv(rr(plN):rr(plN)+countpts(plN)-1) = linspace(-1,1,countpts(plN));
                end
                
                
                %spec2 = [repmat(['%0.', precision, 'f '], 1, 9), '\n'];
                %fprintf(fid, spec2, ptsv);
                fprintf(fid, '%d\n', ptsv);
                
            case {'TRIANGLE','TETRAHEDRON'}
                warning('Only LINES in vtkwrite_pathlines!')
                return

        end     
end
fclose(fid);
if strcmpi(filename,'matlab_export.vtk')
    switch computer
        case {'PCWIN','PCWIN64'}
            !paraview.exe --data='matlab_export.vtk' &
            % Exclamation point character is a shell escape, the rest of the
            % input line will be sent to operating system. It can not take
            % variables, though. The & at the end of line will return control to 
            % Matlab even when the outside process is still running. 
        case {'GLNXA64','MACI64'}
            !paraview --data='matlab_export.vtk' &
    end
end
end

function setdataformat(fid, binaryflag)

if ~binaryflag
    fprintf(fid, 'ASCII\n');
else
    fprintf(fid, 'BINARY\n');
end
end
    