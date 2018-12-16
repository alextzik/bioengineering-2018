%%Script for Question 2
allData=zeros(4,1440000);

load('../dataset/Data_Eval_E_1.mat');
allData(1,:)=data;
spikeTimes1=spikeTimes;
spikeClass1=spikeClass;

load('../dataset/Data_Eval_E_2.mat');
allData(2,:)=data;
spikeTimes2=spikeTimes;
spikeClass2=spikeClass;

load('../dataset/Data_Eval_E_3.mat');
allData(3,:)=data;
spikeTimes3=spikeTimes;
spikeClass3=spikeClass;

load('../dataset/Data_Eval_E_4.mat');
allData(4,:)=data;
spikeTimes4=spikeTimes;
spikeClass4=spikeClass;

k=@(sigma)(1.861248757651653+0.250156158913673./sigma-0.008006893531367./(sigma.*sigma)-2.410871628915119e-05/(sigma.*sigma.*sigma));

%% Q2.1
sigmas=zeros(4,1);
measuredNumSpikes=zeros(4,1);
for i=1:4
    sigmas(i)=median(abs(allData(i,:)))/0.6745;
    %sigmas(i)=std(allData(i,:));
    T=k(sigmas(i))*sigmas(i);
    previousMeasuredSpike=0;
    for m=1:1:1440000
       if (m>1 && (allData(i,m))>T && (allData(i,m-1))<=(allData(i,m)) && (allData(i,m))>=(allData(i,m)) && previousMeasuredSpike==0)
           measuredNumSpikes(i)=measuredNumSpikes(i)+1;
           previousMeasuredSpike=m;
       elseif (m==1 && (allData(i,m))>T && (allData(i,m))>=(allData(i,m)) && previousMeasuredSpike==0)
           measuredNumSpikes(i)=measuredNumSpikes(i)+1;
           previousMeasuredSpike=m;    
       elseif (m==1440000 && (allData(i,m))>T && (allData(i,m-1))<=(allData(i,m)) && previousMeasuredSpike==0)
           measuredNumSpikes(i)=measuredNumSpikes(i)+1;
           previousMeasuredSpike=m;
       elseif (previousMeasuredSpike~=0 && (allData(i,m))<T)
           previousMeasuredSpike=0;
       end
    end
end
%%
spikeTimesEst_1=zeros(measuredNumSpikes(1),1);
spikeTimesEst_2=zeros(measuredNumSpikes(2),1);
spikeTimesEst_3=zeros(measuredNumSpikes(3),1);
spikeTimesEst_4=zeros(measuredNumSpikes(4),1);
for i=1:4
    T=k(sigmas(i))*sigmas(i);
    previousMeasuredSpike=0;
    count=0;
    for m=1:1:1440000
       if ((allData(i,m))>=T && previousMeasuredSpike==0)
           count=count+1;
           if(i==1)
               spikeTimesEst_1(count)=m;
           elseif(i==2)
               spikeTimesEst_2(count)=m;
           elseif(i==3)
               spikeTimesEst_3(count)=m;
           elseif(i==4)
               spikeTimesEst_4(count)=m;
           end
           previousMeasuredSpike=m;
       elseif (previousMeasuredSpike~=0 && (allData(i,m))<T)
           previousMeasuredSpike=0;
       end
    end
end


%% Q2.2

SpikesEst_1=zeros(length(spikeTimesEst_1), 64);
SpikesEst_2=zeros(length(spikeTimesEst_2), 64);
SpikesEst_3=zeros(length(spikeTimesEst_3), 64);
SpikesEst_4=zeros(length(spikeTimesEst_4), 64);

centers1=zeros(length(spikeTimesEst_1),1);
for i=1:1:length(spikeTimesEst_1)
   k=spikeTimesEst_1(i);
   l=0;
   while (allData(1, k-l)<=allData(1, k-l+1))
        l=l+1;
   end
   t=0;
   while ((allData(1, k+t)>=allData(1, k) && allData(1, k+t)>=allData(1, k+t-1)))
        t=t+1;
   end
   if (allData(1, k-l)<=-allData(1, k))
       centers1(i)=k-l;
   else
       centers1(i)=k+t;
   end
   SpikesEst_1(i,:)=allData(1, centers1(i)-31:centers1(i)+32);
end
plot(1:1:64, SpikesEst_1(:,:));


%% Q2.3
correctlyFoundspikes=zeros(4,1);
if (length(spikeTimes1)>length(spikeTimesEst_1))
    timesofSpikesCounted1=zeros(length(spikeTimesEst_1),1);
    realTimesofSpikesCounted=zeros(length(spikeTimesEst_1),1);
else
    timesofSpikesCounted1=zeros(length(spikeTimes1),1);
    realTimesofSpikesCounted=zeros(length(spikeTimes1),1);
end 
k=1;
for i=1:1:length(spikeTimes1)
   d=inf;
   for j=1:1:length(spikeTimesEst_1)
      if(abs(spikeTimesEst_1(j)-spikeTimes1(i))<d && ~ismember(spikeTimesEst_1(j), timesofSpikesCounted1))
          d=abs(spikeTimesEst_1(j)-spikeTimes1(i));
          timesofSpikesCounted1(k)=spikeTimesEst_1(j);
          correctlyFoundspikes(1)=correctlyFoundspikes(1)+1;
          realTimesofSpikesCounted(k)=spikeTimes1(i);
          k=k+1;
      end
   end
end
%%Needs to be done for the other 3 files

%% Q2.4
maxVals1=zeros(length(spikeTimesEst_1),1);
minVals1=zeros(length(spikeTimesEst_1),1);
crossingTimes=zeros(length(spikeTimesEst_1),1);
timeofChange=zeros(length(spikeTimesEst_1),1);
for i=1:1:length(spikeTimesEst_1)
    maxVals1(i)=max(SpikesEst_1(i,:));
    minVals1(i)=min(SpikesEst_1(i,:));
    for j=1:1:size(SpikesEst_1,2)-1
       if(SpikesEst_1(i,j)>=0 && SpikesEst_1(i, j+1)<=0)
          crossingTimes(i)=crossingTimes(i)+1;
       end
    end
    for j=32:1:63
       if((SpikesEst_1(i,j)>=0 && SpikesEst_1(i, j+1)<=0) || SpikesEst_1(i,j)<=0 && SpikesEst_1(i, j+1)>=0)
          timeofChange(i)=j-32;
          break;
       end
    end
end

%Needs to be done for the other 3 files as well


%% Q2.5
for i=1:1:length(spikeTimes1)
   if(~ismember(spikeTimes1(1,i),realTimesofSpikesCounted))
       spikeClass1(i)=0;
   end
end

perc1=MyClassify(dataset1, spikeClass1(find(spikeClass1)))