%{
%% Initialization
files = ['../dataset/Data_Eval_E_1.mat'; '../dataset/Data_Eval_E_2.mat'; '../dataset/Data_Eval_E_3.mat'; '../dataset/Data_Eval_E_4.mat'];

allData=zeros(4,1440000);
spike_Times = cell(4,1);
spike_Class = cell(4,1);

for i=1:1:4
    load(files(i,:));
    allData(i,:) = data;
    spike_Times{i} = spikeTimes;
    spike_Class{i} = spikeClass;
end
%change name for more intuition & clear temp variables
spikeTimes = spike_Times;   
spikeClass = spike_Class;
clear spike_Times spike_Class;

k=@(sigma)(1.861248757651653+0.250156158913673./sigma-0.008006893531367./(sigma.*sigma)-2.410871628915119e-05/(sigma.*sigma.*sigma));
disp("Q2.1")
%% Q2.1
sigmas=zeros(4,1);
measuredNumSpikes=zeros(4,1);
for i=1:1:4
    sigmas(i)=median(abs(allData(i,:)))/0.6745;
    T=k(sigmas(i))*sigmas(i);
    previousMeasuredSpike=0;
    for m=1:1:1440000
       if ((allData(i,m))>T && previousMeasuredSpike==0)
           measuredNumSpikes(i)=measuredNumSpikes(i)+1;
           previousMeasuredSpike=m;
       elseif (previousMeasuredSpike~=0 && (allData(i,m))<T)
           previousMeasuredSpike=0;
       end
    end
end

%spikeTimesEst {estimation of when has a spike appeared} initialization
spikeTimesEst = cell(4,1);
for i=1:1:4
    spikeTimesEst{i} = zeros(measuredNumSpikes(i),1);
end

for i=1:4
    T=k(sigmas(i))*sigmas(i);
    previousMeasuredSpike=0;
    count=0;
    for m=1:1:1440000
       if ((allData(i,m))>=T && previousMeasuredSpike==0)
           count=count+1;
           spikeTimesEst{i}(count)=m;
           previousMeasuredSpike=m;
       elseif (previousMeasuredSpike~=0 && (allData(i,m))<T)
           previousMeasuredSpike=0;
       end
    end
end
disp("Q2.2")
%% Q2.2
%spikeEst {4-cell matrix, containing arrays that display the waveforms of all the measured spikes}
spikesEst = cell(4,1);
for j=1:1:4
    spikesEst{j} = zeros(length(spikeTimesEst{j}),64);
    centers=zeros(length(spikeTimesEst{j}),1);
    for i=1:1:length(spikeTimesEst{j})
       k=spikeTimesEst{j}(i);
       l=0;
       while (allData(1, k-l)<=allData(1, k-l+1))
            l=l+1;
       end
       t=0;
       while ((allData(1, k+t)>=allData(1, k) && allData(1, k+t)>=allData(1, k+t-1)))
            t=t+1;
       end
       if (allData(1, k-l)<=-allData(1, k))
           centers(i)=k-l;
       else
           centers(i)=k+t;
       end
       spikesEst{j}(i,:)=allData(1, centers(i)-31:centers(i)+32);
       %spikeTimesNewEst_1(i) = centers(i);
    end
    figure()
    plot(1:1:64, spikesEst{j}(:,:));
end

disp("Q2.3")
%% Q2.3
correctSpikes = zeros(4,1);
spikesCounted = cell(4,1);

for m=1:1:4
    %initialize spikesCounted {spikes correlated to the real ones} 
    if (length(spikeTimes{m})>length(spikeTimesEst{m}))
        spikesCounted{m}=zeros(length(spikeTimesEst{m}),1);
    else
        spikesCounted{m}=zeros(length(spikeTimes{m}),1);
    end 
    %for every real spike, find one of the measured ones to correlate to
    for i=1:1:length(spikeTimes{m})
       d=inf;
       j=1;
       while(j<=i && j<length(spikeTimesEst{m}))
          %if the spike is the closest one to the currently real examined
          %and it wasn't chosen before, correlate it now
          if(abs(spikeTimesEst{m}(j)-spikeTimes{m}(i))<d && ~ismember(spikeTimesEst{m}(j), spikesCounted{m}))
              d=abs(spikeTimesEst{m}(j)-spikeTimes{m}(i));
              spikesCounted{m}(i) = j;
              correctSpikes(m) = correctSpikes(m)+1;
          end
          j = j+1;
       end
    end
end
%}
load('eval-data.mat');
disp("Q2.4")
%% Q2.4
attr = cell(4,1);
for i=1:1:4
    %let's begin with two attributes
    attr{i} = zeros(length(spikesEst{i}(:,64)),2);
    for j=1:1:length(spikesEst{i}(:,64))
        attr{i}(j,1) = max(spikesEst{i}(j,:));
        for k=1:1:63
            if (spikesEst{i}(j,k)*spikesEst{i}(j,k+1) < 0)
                attr{i}(j,2) = attr{i}(j,2) + 1;
            end
        end
    end
    figure()
    plot(attr{i}(:,1), attr{i}(:,2), 'o');
end
%% Q2.5
acc = zeros(4,1);
for i=1:1:4
   data = zeros(length(spikesCounted{i}),2);
   for j=1:1:length(spikesCounted{i})
       data(j,:) = attr{i}(spikesCounted{i}(j),:);
   end
   acc(i) = MyClassify(data,spikeClass{i}(:));
end

disp(acc);