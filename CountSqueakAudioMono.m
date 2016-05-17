function [  ] = CountSqueakAudioMonoFlac( )
%COUNTSQUEAKAUDIOMONOFLAC Look for a loud, wide squeak
%   Read a directory of flac (or wav) audio files
%   Apply the Matlab spectrogram function
%   Sum the PSD, and then look for a certain type of squeak

AllFiles = dir('201*.flac');
% Computes the number of files in the current directory
NumFiles = length(AllFiles);

if NumFiles==0
    AllFiles = dir('201*.wav');
    % Computes the number of files in the current directory
    NumFiles = length(AllFiles);
end


display('Creating count_out spreadsheet');
display('Processing ');
display(num2str(NumFiles));
filename = 'count_out.xlsx';

results = cell(NumFiles+1,4); % creating a cell array to store the data

results(1,1) = cellstr('Name');
results(1,2) = cellstr('Count');

% empirical investigations show this to be a good threshold value for the 
% squeak of interest
threshold  = -12250;

for i = 1:NumFiles
    % we split the file into two because of memory constraints
    try 
        display(AllFiles(i).name)
        display('Part1');
        samples = [1,1800*44100];

        % read in the audio signal
        [y2,Fs2]=audioread(AllFiles(i).name,samples);
        % use the spectrogram function to get the PSD estimate
        [~,~,~,P] = spectrogram(y2,256,50,256,Fs2);
        % sum the PSD estimate over time
        sumPSDOverTime=sum(10*log10(P));
        % find which areas are over the threshold we set
        hits = sumPSDOverTime>threshold;
        len = length(sumPSDOverTime);
        % now we process the summed PSD
        % a counter to see how wide the squeak is
        count=0;
        % the number of squeaks
        numSqueaks = 0;
        % loop over the whole summed PSD
        for j = 1:len
            if hits(j)==1
                count=count+1;
            else
                % 10 gives us a squeak width of > 0.0421s 
                if count>10
                    numSqueaks = numSqueaks+1;
                end
                count=0;
            end
        end

        % repeat for the second half
        display('Part2');
        samples = [1800*44100,3594*44100];
        % read in the audio signal
        [y2,Fs2]=audioread(AllFiles(i).name,samples);
        % use the spectrogram function to get the PSD estimate
        [~,~,~,P] = spectrogram(y2,256,50,256,Fs2);
        % sum the PSD estimate over time
        sumPSDOverTime=sum(10*log10(P));
        % find which areas are over the threshold we set
        hits = sumPSDOverTime>threshold;
        len = length(sumPSDOverTime);
        % now we process the summed PSD
        % a counter to see how wide the squeak is
        count=0;
        % loop over the whole summed PSD
        for j = 1:len
            if hits(j)==1
                count=count+1;
            else
                % 10 gives us a squeak width of > 0.0421s 
                if count>10
                    numSqueaks = numSqueaks+1;
                end
                count=0;
            end
        end


       results(i+1,1)=cellstr(AllFiles(i).name);
       results(i+1,2)=cellstr(num2str(numSqueaks));

       display(num2str(numSqueaks));
    catch
       display('OOps, an error occurred'); 
    end
end

xlswrite(filename,results)
display('Finished');

end

