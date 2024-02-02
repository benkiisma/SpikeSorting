%% File reading and amplification

clear
clc
read_Intan_RHD2000_file()

data=amplifier_data;

dataset_path = './data';

warning('off', 'all');

clearvars path filename ddata

%% Set the curation parameters

%Curating using the amplitude 
% Set your threshold value and set amp_cur to 1 if you want to use it for curation
amplitude_thr = -15; 
amp_cur = 1;

%Curating using the peak to valley distance
% Set your threshold value and set peak_cur to 1 if you want to use it for curation
peak_to_valley_thr = 40; 
peak_cur = 1;

%Curating using the repolarization time 
% Set your thresholds values and set rep_cur to 1 if you want to use it for curation
min_repolarization_dur = 0.6;%[ms]
max_repolarization_dur = 1; %[ms]
rep_cur = 1;

%Curating using the number of spikes
% Set your threshold value and set nb_spikes_cur to 1 if you want to use it for curation
nb_spikes_thr = 400; 
nb_spikes_cur = 1;

%Curating using the SNR
% Set your threshold value and set SNR_cur to 1 if you want to use it for curation
SNR_thr = 4;
SNR_cur = 1;

%%Plot units 
plot_single_unit = 1;
plot_multi_unit = 1;

%% Parameters settings

userInputs.nChan = length(amplifier_channels);
intan.raw = data;
samplingFrequencies= frequency_parameters.amplifier_sample_rate  ;
disp(['Duration of cut neuro data: ' num2str(size(intan.raw,2)/samplingFrequencies) ' sec'])
time_intan = (1:length(intan.raw)) / samplingFrequencies;
normData=intan.raw-median(intan.raw,1);
normDataZ = zscore(normData,0,2); %Zscore along the lines using the sample standard deviation.

%% Spike Extraction

('EXTRACTING SPIKES');
%warning('off','MATLAB:colon:nonIntegerIndex')  

userInputs.stdev = median(normData/0.6745);
userInputs.thr = 5.0 * userInputs.stdev; %considers the threshold 5 times the stdev
userInputs.nUnits = 3; 
userInputs.twindow = 1e-3; %time window for spikes in s ! (Refractory period)
userInputs.wp = userInputs.twindow * samplingFrequencies; % window - corresponding number of samples
userInputs.tstep = 10e-3 * samplingFrequencies; %timestep in points

Spikes = cell(1, userInputs.nChan); %Creation of a cell that will contain the spikes

for i = 1:userInputs.nChan
    [Spikes{i}.peak , Spikes{i}.loc] = findpeaks(-normDataZ(i,:), samplingFrequencies, 'MinPeakHeight',4,'MinPeakDistance',userInputs.twindow);

    Spikes{i}.peak = - Spikes{i}.peak;
    if  isempty(Spikes{i}.loc)
        Spikes{i}.loc=1;
        Spikes{i}.peak=0;
    end
    if ~mod(i,10)
        disp(['EXTRACTED SPIKES IN ' num2str(i) ' Channels']);
    end
end
clearvars i ans

%% Cleaning spikes

% Remove spikes happening in more than n channels at the same time (noise)
spikeTimes = cell(userInputs.nChan,1); 
idx = cell(userInputs.nChan,1);
countRemoved = 0;
for jj = 1:10000
    for i = 1:userInputs.nChan
        if isstruct(Spikes{1,i})
            spikeTimes{i,:} = round((Spikes{1,i}.loc*samplingFrequencies+2 ... 
                            -mod(Spikes{1,i}.loc*samplingFrequencies,4))/4);
        end

    end
    A = cat(2,spikeTimes{1:end,1}); % Get all spike times together
    [tcommPeak,timesPeak] = mode(A); % Find spike time that occurs the most

    % if it's a spike in more than n channels at the same time, remove
    if timesPeak >= 5 && tcommPeak>0

    disp("Common spike detected ...Removing...")

    % Find the times at which each peak occurred
   
        for k = 1:userInputs.nChan
            if ~isempty(Spikes{1,k})
                [temp,~] = find(round((Spikes{1,k}.loc(:)*samplingFrequencies ...
                           +2 -mod(Spikes{1,k}.loc(:)*samplingFrequencies,4))/4) == tcommPeak);
                if ~isempty(temp)
                   idx{k,:} = temp;
                else
                    idx{k,:} = 0;
                end
            else
            idx{k,:} = 0;
            end
        end

         % Remove them 
        for l = 1:length(idx)
            if sum(idx{l,:})>0
                Spikes{1,l}.peak(idx{l,:}) = [];
                Spikes{1,l}.loc(idx{l,:}) = [];
            end
        end
        countRemoved = countRemoved+1;
    else % If it's less than 40 channels then stop finding them
        disp([mat2str(countRemoved) ' Spikes removed'])
        break
    end
end    

clearvars jj k i l timesPeak tcommPeak A idx 

disp("---------------------------------------")

%% Waveform computation

time=1:1:71;
time=(time./30000).*1000; %% About 2.4ms windows
Waveform_inter = cell(userInputs.nChan,1);

for r=1:userInputs.nChan
    temp=zeros(length(Spikes{1,r}.loc),71);
    for h=2:length(Spikes{1,r}.loc)-1
        temp(h-1,:)=normData(r, (round(Spikes{1,r}.loc(1,h)*samplingFrequencies-30)):1:(round(Spikes{r}.loc(1,h)*samplingFrequencies+40)));
    end
    Waveform_inter{r,1}=temp;
end

clearvars r h temp

%% SNR Computation

clearvars SNR
Noise_sig=1:1:46;
Noise_sig=(Noise_sig./30000).*1000;
SNR=zeros(32,1);

for r=1:32
    noise=zeros(length(Spikes{1,r}.loc),46);
    for h=2:length(Spikes{1,r}.loc)-1
        noise(h-1,:)=normData(r, (round(Spikes{1,r}.loc(1,h)*samplingFrequencies+45)):1:(round(Spikes{r}.loc(1,h)*samplingFrequencies+90)));
    end
    Noise{r,1}=noise;
    noise_std=mean(std(Noise{r,1}));
    signal=min(mean(Waveform_inter{r,1}));
    SNR(r,1)=-(signal./noise_std);
end

clearvars r h Noise_sig noise Noise noise_std signal

%% SNR threshold 
%Take out units whose SNR is below the threhold

if SNR_cur == 1 % If the SNR curation is selected
    count = 0;
    remove = 0;
    for channel=1:length(Waveform_inter)
        if SNR(channel,1) > SNR_thr
            count = count + 1; %Count the number of good units
        end
    end

    Waveform = cell(count,1);
    interTime = cell(count,1);
    
    for channel=1:length(Waveform_inter)
        if SNR(channel,1) > SNR_thr
            Waveform{channel-remove,1} = Waveform_inter{channel,1};%Selecting only the good units
            interTime{channel-remove,1} = spikeTimes{channel,1};
        else
            remove = remove + 1;
        end
    end
else
    Waveform = Waveform_inter;
    interTime = spikeTimes;
end

disp("---------------------------------------")
disp([mat2str(remove) ' Spikes removed'])
disp("---------------------------------------")

clearvars channel count remove
%% adjacency matrix
% Reordering of the channels based on the electrodes adjacency

adj_mat = cell(userInputs.nChan,1);
adj_mat{1,1} = cat(1,Waveform{16,1},Waveform{1,1});
adj_mat{2,1} = cat(1,Waveform{16,1},Waveform{1,1},Waveform{17,1});
adj_mat{3,1} = cat(1,Waveform{1,1},Waveform{17,1},Waveform{32,1});
adj_mat{4,1} = cat(1,Waveform{17,1},Waveform{32,1},Waveform{18,1});
adj_mat{5,1} = cat(1,Waveform{32,1},Waveform{18,1},Waveform{31,1});
adj_mat{6,1} = cat(1,Waveform{18,1},Waveform{31,1},Waveform{3,1});
adj_mat{7,1} = cat(1,Waveform{31,1},Waveform{3,1},Waveform{19,1},Waveform{13,1},Waveform{28,1});
adj_mat{8,1} = cat(1,Waveform{3,1},Waveform{19,1},Waveform{28,1},Waveform{13,1});
adj_mat{9,1} = cat(1,Waveform{3,1},Waveform{13,1},Waveform{28,1},Waveform{19,1});
adj_mat{10,1} = cat(1,Waveform{19,1},Waveform{13,1},Waveform{3,1},Waveform{28,1},Waveform{15,1});
adj_mat{11,1} = cat(1,Waveform{28,1},Waveform{15,1},Waveform{30,1});
adj_mat{12,1} = cat(1,Waveform{15,1},Waveform{30,1},Waveform{2,1});
adj_mat{13,1} = cat(1,Waveform{30,1},Waveform{2,1},Waveform{20,1});
adj_mat{14,1} = cat(1,Waveform{2,1},Waveform{20,1},Waveform{6,1});
adj_mat{15,1} = cat(1,Waveform{20,1},Waveform{6,1},Waveform{22,1},Waveform{10,1},Waveform{26,1});
adj_mat{16,1} = cat(1,Waveform{6,1},Waveform{22,1},Waveform{26,1},Waveform{10,1});
adj_mat{17,1} = cat(1,Waveform{6,1},Waveform{10,1},Waveform{26,1},Waveform{22,1});
adj_mat{18,1} = cat(1,Waveform{22,1},Waveform{10,1},Waveform{6,1},Waveform{26,1},Waveform{14,1});
adj_mat{19,1} = cat(1,Waveform{26,1},Waveform{14,1},Waveform{4,1});
adj_mat{20,1} = cat(1,Waveform{14,1},Waveform{4,1},Waveform{28,1});
adj_mat{21,1} = cat(1,Waveform{4,1},Waveform{28,1},Waveform{21,1});
adj_mat{22,1} = cat(1,Waveform{28,1},Waveform{21,1},Waveform{25,1},Waveform{12,1},Waveform{8,1});
adj_mat{23,1} = cat(1,Waveform{21,1},Waveform{25,1},Waveform{8,1},Waveform{12,1});
adj_mat{24,1} = cat(1,Waveform{21,1},Waveform{12,1},Waveform{8,1},Waveform{25,1});
adj_mat{25,1} = cat(1,Waveform{25,1},Waveform{12,1},Waveform{21,1},Waveform{8,1},Waveform{5,1});
adj_mat{26,1} = cat(1,Waveform{8,1},Waveform{5,1},Waveform{24,1});
adj_mat{27,1} = cat(1,Waveform{5,1},Waveform{24,1},Waveform{23,1});
adj_mat{28,1} = cat(1,Waveform{24,1},Waveform{23,1},Waveform{27,1},Waveform{11,1},Waveform{7,1});
adj_mat{29,1} = cat(1,Waveform{23,1},Waveform{27,1},Waveform{7,1},Waveform{11,1});
adj_mat{30,1} = cat(1,Waveform{23,1},Waveform{11,1},Waveform{7,1},Waveform{27,1});
adj_mat{31,1} = cat(1,Waveform{27,1},Waveform{11,1},Waveform{23,1},Waveform{7,1},Waveform{9,1});
adj_mat{32,1} = cat(1,Waveform{7,1},Waveform{9,1});

adj_time = cell(userInputs.nChan,1);
adj_time{1,1} = cat(2,interTime{16,1},interTime{1,1});
adj_time{2,1} = cat(2,interTime{16,1},interTime{1,1},interTime{17,1});
adj_time{3,1} = cat(2,interTime{1,1},interTime{17,1},interTime{32,1});
adj_time{4,1} = cat(2,interTime{17,1},interTime{32,1},interTime{18,1});
adj_time{5,1} = cat(2,interTime{32,1},interTime{18,1},interTime{31,1});
adj_time{6,1} = cat(2,interTime{18,1},interTime{31,1},interTime{3,1});
adj_time{7,1} = cat(2,interTime{31,1},interTime{3,1},interTime{19,1},interTime{13,1},interTime{28,1});
adj_time{8,1} = cat(2,interTime{3,1},interTime{19,1},interTime{28,1},interTime{13,1});
adj_time{9,1} = cat(2,interTime{3,1},interTime{13,1},interTime{28,1},interTime{19,1});
adj_time{10,1} = cat(2,interTime{19,1},interTime{13,1},interTime{3,1},interTime{28,1},interTime{15,1});
adj_time{11,1} = cat(2,interTime{28,1},interTime{15,1},interTime{30,1});
adj_time{12,1} = cat(2,interTime{15,1},interTime{30,1},interTime{2,1});
adj_time{13,1} = cat(2,interTime{30,1},interTime{2,1},interTime{20,1});
adj_time{14,1} = cat(2,interTime{2,1},interTime{20,1},interTime{6,1});
adj_time{15,1} = cat(2,interTime{20,1},interTime{6,1},interTime{22,1},interTime{10,1},interTime{26,1});
adj_time{16,1} = cat(2,interTime{6,1},interTime{22,1},interTime{26,1},interTime{10,1});
adj_time{17,1} = cat(2,interTime{6,1},interTime{10,1},interTime{26,1},interTime{22,1});
adj_time{18,1} = cat(2,interTime{22,1},interTime{10,1},interTime{6,1},interTime{26,1},interTime{14,1});
adj_time{19,1} = cat(2,interTime{26,1},interTime{14,1},interTime{4,1});
adj_time{20,1} = cat(2,interTime{14,1},interTime{4,1},interTime{28,1});
adj_time{21,1} = cat(2,interTime{4,1},interTime{28,1},interTime{21,1});
adj_time{22,1} = cat(2,interTime{28,1},interTime{21,1},interTime{25,1},interTime{12,1},interTime{8,1});
adj_time{23,1} = cat(2,interTime{21,1},interTime{25,1},interTime{8,1},interTime{12,1});
adj_time{24,1} = cat(2,interTime{21,1},interTime{12,1},interTime{8,1},interTime{25,1});
adj_time{25,1} = cat(2,interTime{25,1},interTime{12,1},interTime{21,1},interTime{8,1},interTime{5,1});
adj_time{26,1} = cat(2,interTime{8,1},interTime{5,1},interTime{24,1});
adj_time{27,1} = cat(2,interTime{5,1},interTime{24,1},interTime{23,1});
adj_time{28,1} = cat(2,interTime{24,1},interTime{23,1},interTime{27,1},interTime{11,1},interTime{7,1});
adj_time{29,1} = cat(2,interTime{23,1},interTime{27,1},interTime{7,1},interTime{11,1});
adj_time{30,1} = cat(2,interTime{23,1},interTime{11,1},interTime{7,1},interTime{27,1});
adj_time{31,1} = cat(2,interTime{27,1},interTime{11,1},interTime{23,1},interTime{7,1},interTime{9,1});
adj_time{32,1} = cat(2,interTime{7,1},interTime{9,1});

clearvars Waveform interTime

%% Clustering Running K-means

nb_iter = 1; %Setting the number of clustering iteration
k_values = zeros(length(adj_mat),nb_iter);
labels_value = cell(length(adj_mat),1);

for iter=1:nb_iter
    disp(iter);
    for channel=1:length(adj_mat)
        X = min(adj_mat{channel,1},1); % Set the dataset to cluster
        eva = evalclusters(X,'kmeans','silhouette','klist',1:2); % Assess the optimal k number of clusters
        K = eva.OptimalK;
        k_values(channel,iter) = K;

        labels = kmeans(X,K); % Applying k-means
        
        labels_inter = labels; 
        labels_value{channel,1}(:,iter) = labels_inter; %Storing the labels into a matrix
    end
end

clearvars iter  channel X K labels_inter

%% Mean over the iterations
%Taking the mean over the nb_iter iterations to have the most representative result.

labels_final = cell(length(adj_mat),1);
for channel=1:length(adj_mat)
    labels_final{channel,1} = round(mean(labels_value{channel,1}(:,:),2));
end
k_final = round(mean(k_values,2));

clearvars channel labels_value
%% Form the units by taking out unit with less than 20% representation of the data
% If there is a cluster that doesn't represent the data correctly (Too few datapoints in it), we take it out.

Units = cell(length(labels_final),max(k_final));
Units_time = cell(length(labels_final),max(k_final));
count=0;

%We link the spikes to their respective cluster

for channel=1:userInputs.nChan
    for signal=1:length(labels_final{channel,1})
        for k=1:k_final(channel)
            if labels_final{channel,1}(signal)==k
                Units{channel,k} = [Units{channel,k};adj_mat{channel, 1}(signal,:)];
                Units_time{channel,k} = [Units_time{channel,k};adj_time{channel, 1}(:,signal)];
                
                break;
            end
        end 
    end
    for k=1:max(k_final)
        if size(Units{channel,k},1) > 0.2*length(labels_final{channel,1})
            count = count+1;     
        end
    end
end

Units_inter = cell(count,1);
Units_time_inter = cell(count,1);
remove = 0;

%We remove the non-representative clusters and align the others

for k=1:(max(k_final))
    for unit=1:userInputs.nChan
        if size(Units{unit,k},1) > 0.2*length(labels_final{unit,1})
            Units_inter{unit-remove+((k-1)*userInputs.nChan),1} = Units{unit,k};
            Units_time_inter{unit-remove+((k-1)*userInputs.nChan),1} = Units_time{unit,k};
        else 
            remove = remove + 1;
        end
    end
end

clearvars unit remove k channel signal unit count

%% Curation
%Taking out spikes that do not meet the choosen thresholds

%Amplitude
if amp_cur ==1 %If the amplitude curation is selected
    for channel=1:length(Units_inter)
        inter = [];
        inter_time = [];
        for spike=1:size(Units_inter{channel,1},1)
            if min(Units_inter{channel,1}(spike,:)) < amplitude_thr
                inter = [inter;Units_inter{channel,1}(spike,:)]; %Selecting good spikes
                inter_time = [inter_time;Units_time_inter{channel,1}(spike,1)]; %Selecting good spikes

            end
        end
        Units_inter{channel,1} = inter;
        Units_time_inter{channel,1} = inter_time;
    end
end
clearvars inter inter_time

%Peak to valley
if peak_cur == 1 %If the peak curation is selected
    for channel=1:length(Units_inter)
        inter = [];
        inter_time = [];
        for spike=1:size(Units_inter{channel,1},1)
            [min_value,index_min] = min(Units_inter{channel,1}(spike,:));
            [max_value,index_max] = max(Units_inter{channel,1}(spike,:));
            if (max_value-min_value) > peak_to_valley_thr
                inter = [inter;Units_inter{channel,1}(spike,:)];
                inter_time = [inter_time;Units_time_inter{channel,1}(spike,1)];
            end
        end
        Units_inter{channel,1} = inter;
        Units_time_inter{channel,1} = inter_time;
    end
end

%Repolarization slope
if rep_cur == 1 %If the repolarization curation is selected
    for channel=1:length(Units_inter)
        inter = [];
        inter_time = [];
        for spike=1:size(Units_inter{channel,1},1)
            [min_value,index_min] = min(Units_inter{channel,1}(spike,:));
            [max_value,index_max] = max(Units_inter{channel,1}(spike,:));
            if min_repolarization_dur < (((index_max-index_min)*2.4)/71) < max_repolarization_dur
                inter = [inter;Units_inter{channel,1}(spike,:)];
                inter_time = [inter_time;Units_time_inter{channel,1}(spike,1)];
            end
        end
        Units_inter{channel,1} = inter;
        Units_time_inter{channel,1} = inter_time;
    end
end

% Number of spikes
if nb_spikes_cur == 1 %If the curation based on the nb od spikes is selected
    count = 0;
    remove = 0;
    for channel=1:length(Units_inter)
        if size(Units_inter{channel,1},1) > nb_spikes_thr
            count = count + 1;
        end
    end
    inter = cell(count,1);
    inter_time = cell(count,1);
    for channel=1:length(Units_inter)
        if size(Units_inter{channel,1},1) > nb_spikes_thr
            inter{channel-remove,1} = Units_inter{channel,1};
            inter_time{channel-remove,1} = Units_time_inter{channel,1};
        else
            remove = remove + 1;
        end
    end
    Units_inter = inter;
    Units_time_inter = inter_time;
end

clearvars channel spike inter min_value index_min max_value index_max

%% plot clusters
%Choose the units that you want to see
% Be aware that sometimes only 1 cluster will be resulting from the
% clustering

%Unit = 1;
%spike1 = 1;
%spike2 = 1;

%figure;
%plot(time, Units{Unit,1}(spike1,:), '.-'); hold on
%plot(time, Units{Unit,2}(spike2,:),'r.-'); hold on

%clearvars Unit spike1 spike2
%% Plot
%We plot the resulting units to visualize them.

% figure;
% for i=1:length(Units_inter)
%     subplot(8,9,i)
%     clearvars d s t c1 c2
%     
%     d=mean(Units_inter{i,1});
%     s=std(Units_inter{i,1});
%     c1=d+s;
%     c2=d-s;
%     t=[time,fliplr(time)];
%     inBetween = [c1, fliplr(c2)];
%     fill(t, inBetween, [0.7 0.8 0.2], 'faceAlpha', 0.7, 'LineStyle', 'none');
%     hold on;
%     plot(time, d, 'k', 'LineWidth', 2);
%     % if i~=1
%     %yticks([-50 0 50 100])
%     ylim([-50 50])
%     title(['E',num2str(i)])
% end
% 
% clearvars i d s c1 c2 t 

%% Merging clusters 
%We merge clusters that have a high similarity ratio (Using the max and min values)
eps = 1;
Merged = [0];
Units_inter_copy = Units_inter;
Units_time_inter_copy = Units_time_inter;

while isempty(Merged)== 0
    Merged = [];
    for i=1:length(Units_inter_copy)
        if ismember(i,Merged) == 0
            for j=(i+1):length(Units_inter_copy)
                if ismember(j,Merged) == 0
                    if abs(min(mean(Units_inter_copy{i,1},1))-min(mean(Units_inter_copy{j,1},1))) < eps  %If the difference between the minimums are within epsilon
                        if abs(max(mean(Units_inter_copy{i,1},1))-max(mean(Units_inter_copy{j,1},1))) < eps % If the difference between the maximums are within epsilon
                            disp(['Merging cluster ' num2str(i) ' with cluster ' num2str(j)]);
                            Units_inter_copy{i,1}=cat(1,Units_inter_copy{i,1}(:,:),Units_inter_copy{j,1}(:,:));
                            Units_time_inter_copy{i,1}=cat(1,Units_time_inter_copy{i,1}(:,:),Units_time_inter_copy{j,1}(:,:));
                            Merged = [Merged;j];
                            break;    
                        end
                    end
                end
            end
        end
    end

    count = length(Units_inter_copy) - length(Merged);
    Units_final = cell(count,1);
    Units_time_final = cell(count,1);
    off = 0;

    for i=1:length(Units_inter_copy)
        if ismember(i,Merged) == 0
            Units_final{i-off,1} = Units_inter_copy{i,1};
            Units_time_final{i-off,1} = Units_time_inter_copy{i,1};
        elseif ismember(i,Merged) == 1
            off = off + 1;
        end
    end
    Units_inter_copy = Units_final;
    Units_time_inter_copy = Units_time_final;
end       

clearvars count i Units_inter_copy Units_time_inter_copy off j Merged 
clearvars Units_time_inter Units_time Units_inter 

%% Plot
% We plot again to visualize the results of merging.

figure;
single_unit_inter = [];
multi_unit_inter = [];
for i=1:length(Units_final)
    subplot(8,8,i)
    clearvars d s t c1 c2
    
    d=mean(Units_final{i,1});
    s=std(Units_final{i,1});
    c1=d+s;
    c2=d-s;
    t=[time,fliplr(time)];
    inBetween = [c1, fliplr(c2)];
    fill(t, inBetween, [0.7 0.8 0.2], 'faceAlpha', 0.7, 'LineStyle', 'none');
    hold on;
    plot(time, d, 'k', 'LineWidth', 2);
    %if i~=1
    %yticks([-50 0 50 100])
    ylim([-50 50])
    title(['E',num2str(i)])
    
    % We determine which unit correspond to MUA and which to SUA based on std
    if mean(c1-c2) > 25
        multi_unit_inter = [multi_unit_inter; i];
    else
        single_unit_inter = [single_unit_inter; i];
    end
end

clearvars i d s c1 c2 t 

%% Determination of multiunit activity
%Putting the single-units in a matrix 
    
single_unit = cell(length(single_unit_inter),1);
multi_unit = cell (length(multi_unit_inter),1);

for i=1:length(multi_unit_inter)
    multi_unit{i,1} = Units_final{multi_unit_inter(i),1};
end

for i=1:length(single_unit_inter)
    single_unit{i,1} = Units_final{single_unit_inter(i),1};
end

clearvars i
%% Display results

disp('----------------------------');
disp([num2str(length(Units_final)) ' Units were found']);
disp("From which : ");
disp([num2str(length(multi_unit)) ' were multi-units']);
disp([num2str(length(single_unit)) ' were single-units']);
disp('----------------------------');   

%% Plot fig

% Plot Single unit activity and multi-unit activity

% %Plot the single_units

if plot_single_unit == 1
    figure;
    for i=1:length(single_unit)
        
        subplot(6,6,i)
        plot(time,single_unit{i,1},'LineWidth' ,0.5, 'Color', [0.7 0 0 0.075])
        ylim([-40 40])
        hold on
        title(['E',num2str(i) '- Single unit Activity' ])
        plot(time,mean(single_unit{i,1}),'k', 'LineWidth',2)
    end
end

% Plot multi-unit

if plot_multi_unit == 1
    figure;
    for i=1:length(multi_unit)
        
        subplot(6,6 ,i)
        plot(time,multi_unit{i,1},'LineWidth' ,0.5, 'Color', [0.7 0 0 0.075])
        ylim([-40 40])
        hold on
        title(['E',num2str(i) '- Multi-unit Activity'])
        plot(time,mean(multi_unit{i,1}),'k', 'LineWidth',2)
        
    end
end
 
 clearvars i
