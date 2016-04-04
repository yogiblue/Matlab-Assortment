function [ ] = BatchRoundness( )
%BATCHROUNDNESS Processes a directory of images for roundness
%   this is what it does
AllFiles = dir('*.jpg');
 
% Computes the number of files in the current directory
NumPhotos = length(AllFiles);

display('Creating roundness_out spreadsheet');
display('Processing ');
display(num2str(NumPhotos));
filename = 'roundness_out.xlsx';

results = cell(NumPhotos+1,5); % creating a cell array to store the names of the images

results(1,1) = cellstr('Name');
results(1,2) = cellstr('Roundness');
results(1,3) = cellstr('Area');
results(1,4) = cellstr('Xpos');
results(1,5) = cellstr('Ypos');

for i = 1:NumPhotos
   results(i+1,1)=cellstr(AllFiles(i).name);
   [Roundness, Area, Centroid] = mCircle(AllFiles(i).name);
   results(i+1,2)=cellstr(num2str(Roundness));
   results(i+1,3)=cellstr(num2str(Area));
   results(i+1,4)=cellstr(num2str(Centroid(1)));
   results(i+1,5)=cellstr(num2str(Centroid(2)));
end

xlswrite(filename,results)
display('Finished');


end

