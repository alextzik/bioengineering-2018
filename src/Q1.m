%% Initialization
files = ['../dataset/Data_Test_1.mat';'../dataset/Data_Test_2.mat';'../dataset/Data_Test_3.mat';'../dataset/Data_Test_4.mat'; ...
        '../dataset/Data_Test_5.mat'; '../dataset/Data_Test_6.mat'; '../dataset/Data_Test_7.mat'; '../dataset/Data_Test_8.mat'];

spikeNumbers=zeros(8,1);
allData=zeros(8,1440000);

%load the samples and num of spikes from all the data files
for i=1:1:8
    load(files(i,:));
    allData(i,:) = data;
    spikeNumbers(i) = spikeNum;
end

%% Q1.1
for i=1:1:8
   figure();
   plot(linspace(1,10000, 10000), allData(i, 1:10000));
end
   
%% Q1.2
s=@(x)(median(abs(x))/0.6745);  %lamda function for computing variance 
numOfks=2001;                   %number of different values of k to test
sigmas=zeros(8,1);
measuredNumSpikes=zeros(8,numOfks);

for i=1:1:8
   sigmas(i)=s(allData(i,:));
   count=1;
   for k=2:15/2000:17
       sameSpike=0;
       T=k*sigmas(i);
       for m=1:1:1440000    %for all the samples
           %every distinct spike may cover more than one samples, so we
           %avoid counting more spikes due to this fact
           if ((allData(i,m))>T && sameSpike==0)  
               %we count a new spike only if we have left behind the
               %previous one {samespike==0}
               measuredNumSpikes(i,count)=measuredNumSpikes(i,count)+1;
               sameSpike=m;
           elseif (sameSpike~=0 && (allData(i,m))<T) 
               %search for the next spike in the next loop iteration
               sameSpike=0;
           end
       end
       count=count+1;
   end
   ks=2:15/2000:17;
   %figure();
   %plot(ks, measuredNumSpikes(i,:));
    
end

%find the best estimation and save the value of k that gives it
closestK=zeros(8,1);
ds=ones(8,1);
for i=1:1:8
    ds(i)=inf;
    for k=1:1:2001
       if (abs(measuredNumSpikes(i, k)-spikeNumbers(i))<ds(i))
            ds(i)=abs(measuredNumSpikes(i, k)-spikeNumbers(i));
            closestK(i)=2+(k-1)*15/2000;
       end
    end
end

figure();
plot(sigmas, closestK,'o');
hold on;

%% Q1.3
%regression for specifying a function estimation for the s-k relationship
fun=@(b,x)(b(1)+b(2)./x+b(3)./(x.*x)+b(4)./(x.*x.*x));
beta0=[0 0 0 0];
beta=nlinfit(sigmas, closestK, fun, beta0);

plot(sigmas, fun(beta,sigmas));
