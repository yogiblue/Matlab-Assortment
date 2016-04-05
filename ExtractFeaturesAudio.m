function [  ] = ExtractFeaturesAudio(outputfile)
%EXTRACTFEATURESAUDIO Extract spectral and temporal features from an audio
%sample
%   Takes an excel filename (outputfile) as an argument
%   Read the current directory of flac (or wav) audio files
%   Filter and threshold the signal to remove noise
%   Calculate fft and then use this for peak frequency
%   Calculate psd and work out the mean power of signal
%   Calculate ar coefficients if required
%   Calculate zero crossing estimate
%   Calculate mean frequency
%   Calculate flatness
%   Calculate 85 percent rolloff point
%   Attempt cepstrum estimate - doesn't seem to work
%   Output the results into the outputfile

AllFiles = dir('201*.flac');
% Computes the number of files in the current directory
NumFiles = length(AllFiles);

if NumFiles==0
    AllFiles = dir('201*.wav');
    % Computes the number of files in the current directory
    NumFiles = length(AllFiles);
end


display('Creating features_out spreadsheet');
display('Processing ');
display(num2str(NumFiles));
filename = outputfile;

%thefilter = load('Beetle_HighPass_1kHz');
%thefilter = load('Beetle_HighPass_3kHz');

% thefilter = load('Beetle_BandPass_1_5kHz');
 thefilter = load('Beetle_BandPass_5_10kHz');
%thefilter = load('Beetle_BandPass_10_15kHz');
%thefilter = load('Beetle_BandPass_15_20kHz');

% use modelSize>0 if you want to calculate ar coefficients
modelSize = 0;
results = cell(NumFiles+1,1+modelSize+8); % creating a cell array to store the data

results(1,1) = cellstr('Name');

factor = 1;
noise_thres = 0.015;
noise_thres = 0;

resCount=0;
results(1,1) = cellstr('Name');
results(1,2+modelSize) = cellstr('Max Freq');
results(1,3+modelSize) = cellstr('Time');
results(1,4+modelSize) = cellstr('Avg Power');
results(1,5+modelSize) = cellstr('Mean Freq');
results(1,6+modelSize) = cellstr('Fund Freq (crossings)');
results(1,7+modelSize) = cellstr('Fund Freq (cepstrum)');
results(1,8+modelSize) = cellstr('Spectral roll off');
results(1,9+modelSize) = cellstr('Spectral flatness');
resCount = resCount+1;

for i = 1:NumFiles
    try 
        display('Processing')
        display(AllFiles(i).name)
        % read in the audio signal
        [x,Fs]=audioread(AllFiles(i).name);

        %meanSignal = mean(abs(x));
        %x = x * 0.03 / meanSignal;

        %meanSignal = mean(abs(x));
               
         %data = abs(x(:,1));
         data = x;
         amp_data = data.*factor;
         
         amp_data(amp_data>1)=1;
         amp_data(amp_data<-1)=-1;
         %amp_data = x;
         filt_data = filter(thefilter.myfilt,amp_data);
         filt_data (filt_data<=noise_thres & filt_data>0) = 0;
         filt_data (filt_data>=-noise_thres & filt_data<0) = 0;
                
         x_raw = x;
         x = filt_data;
         
      
        fftLength = length(x);  % Always make sure to be at least as long as your data.
        xdft = fft(x,fftLength);
        maxAmp = max(abs(xdft));
        freq = [0:fftLength-1].*(Fs/fftLength); % This is your total freq-axis
        freqsYouCareAbout = freq(freq < Fs/2);  % You only care about either the pos or neg
        % frequencies, since they are redundant for
        % a real signal.
        xdftYouCareAbout = abs(xdft(1:round(fftLength/2))); % Take the absolute magnitude.
        [maxVal, index] = max(xdftYouCareAbout); % maxVal is your (un-normalized) maximum amplitude
        %skip low frequency (105 means 1khz)
        [maxVal, index] = max(xdftYouCareAbout(105:length(xdftYouCareAbout),:));
        maxFreq = freqsYouCareAbout(index+104); % This is the frequency of your dominant signal. Skip low frequency

        
        %calculate the power in db/Hz
        N = length(x_raw);
        xdft = xdft(1:round(N/2)+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x):Fs/2;
        meanSignal = mean(10*log10(psdx));
         
         if modelSize>0
            arcoeffs = arburg(x,modelSize);
            [Pxx,F] = pburg(x,modelSize,[],Fs);
             % plot(F,10*log10(Pxx))
            % plot(arcoeffs(2:modelSize+1))
            % arcoeffs = arburg(filt_data,modelSize);
         end
         
        %zero crossing estimate of fundamental frequency
        xin = abs(x_raw);
        xin = xin-mean(xin);
        x2 = zeros(length(xin),1);
        x2(1:length(x_raw)-1) = xin(2:length(x_raw));
        zc=length(find((xin>0 & x2<0) | (xin<0 & x2>0)));
        F0=0.5*Fs*zc/length(x_raw);
        fprintf('Zero-crossing F0 estimate is %3.2f Hz.\n',F0)
        
        % calculate center of gravity - mean frequency
        sumy = 0;
        total_sum = 0;
        for ii=1:length(xdftYouCareAbout)
            sumy = ii * abs(xdftYouCareAbout(ii)) + sumy;
            total_sum = total_sum + abs(xdftYouCareAbout(ii));
        end

        mean_freq_index = round(sumy/total_sum);
        mean_freq = freqsYouCareAbout(mean_freq_index);
        
        %flatness
        flatness = geomean(xdftYouCareAbout) / mean(xdftYouCareAbout);
        
        %rolloff
        eightyfive  = 0.85*sumy;
        sumy=0;
        for ii=1:length(xdftYouCareAbout)
            sumy = ii * abs(xdftYouCareAbout(ii)) + sumy;
            if sumy>=eightyfive
                rolloffpoint = ii;
                break;
            end
        end
        rolloff = freqsYouCareAbout(rolloffpoint);
        
        %cepstrum estimate of fundamental frequenct
        c = cceps(x_raw);        
        dt = 1/Fs;
        t = 0:dt:length(x)*dt-dt;
        [~,I] = max(c(12:500));
        fprintf('Complex cepstrum F0 estimate is %3.2f Hz. mx(%3.2f) min(%3.2f)\n', 1/(t(I+12)), 1/(t(12)), 1/(t(500)))
        cepfreq = 1/(t(I+12));
        
        if maxFreq>00 && mean(x(1000:1500))~=0
            % plot(F,10*log10(Pxx))
            results(resCount+1,1)=cellstr(AllFiles(i).name);
       
            if modelSize>0
                for j=1:modelSize
                    results(resCount+1,1+j)=cellstr(num2str(arcoeffs(1+j)));
                end
            end
        
            results(resCount+1,1+modelSize+1)=cellstr(num2str(maxFreq));
            results(resCount+1,1+modelSize+2)=cellstr(num2str(length(x)/Fs));
            results(resCount+1,1+modelSize+3)=cellstr(num2str(meanSignal));
            results(resCount+1,1+modelSize+4)=cellstr(num2str(mean_freq));
            results(resCount+1,1+modelSize+5)=cellstr(num2str(F0));
            results(resCount+1,1+modelSize+6)=cellstr(num2str(cepfreq));
            results(resCount+1,1+modelSize+7)=cellstr(num2str(rolloff));
            results(resCount+1,1+modelSize+8)=cellstr(num2str(flatness));
            resCount=resCount+1;
        end
       
    catch
       display('OOps, an error occurred'); 
    end
end

xlswrite(filename,results)
display('Finished');



end


