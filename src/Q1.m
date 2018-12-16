%%Answers to the questions of the first Topic
spikeNumbers=zeros(8,1);
allData=zeros(8,1440000);

load('../Data/Data_Test_1.mat');
allData(1,:)=data;
spikeNumbers(1)=spikeNum;

load('../Data/Data_Test_2.mat');
allData(2,:)=data;
spikeNumbers(2)=spikeNum;

load('../Data/Data_Test_3.mat');
allData(3,:)=data;
spikeNumbers(3)=spikeNum;

load('../Data/Data_Test_4.mat');
allData(4,:)=data;
spikeNumbers(4)=spikeNum;

load('../Data/Data_Test_5.mat');
allData(5,:)=data;
spikeNumbers(5)=spikeNum;

load('../Data/Data_Test_6.mat');
allData(6,:)=data;
spikeNumbers(6)=spikeNum;

load('../Data/Data_Test_7.mat');
allData(7,:)=data;
spikeNumbers(7)=spikeNum;

load('../Data/Data_Test_8.mat');
allData(8,:)=data;
spikeNumbers(8)=spikeNum;

%%Q1.1
for i=1:1:8
   figure();
   plot(linspace(1,10000, 10000), allData(i, 1:10000));
end

   
%%Q1.2
s=@(x)(median(abs(x))/0.6745);
numOfks=2001;
sigmas=zeros(8,1);
measuredNumSpikes=zeros(8,numOfks);

%{
measuredNumofSpikes=zeros(8,1);
bestKs=zeros(8,1);
for i=1:1:8
    k=0;
    low=2;
    high=17;
    while(abs(measuredNumofSpikes(i, 1)-spikeNumbers(i))>3)
          measuredNumofSpikes(i, 1)
          k=(low+high)/2;
          sameSpike=0;
          T=k*sigmas(i);
          measuredNumofSpikes(i,1)=0;
          for m=1:1:1440000
                if (abs(allData(i,m))>T && sameSpike==0)
                        measuredNumofSpikes(i,1)=measuredNumofSpikes(i,1)+1;
                        sameSpike=m;
                elseif (sameSpike~=0 && abs(allData(i,m))<T)
                        sameSpike=0;
                end
          end
          if(measuredNumofSpikes(i, 1)>spikeNumbers(i))
              high=k;
          end
          if(measuredNumofSpikes(i, 1)<spikeNumbers(i))
              low=k;
          end
    end
    bestKs(i)=k;
end
%}


for i=1:1:8
   sigmas(i)=s(allData(i,:));
   count=1;
   for k=2:15/2000:17
       sameSpike=0;
       T=k*sigmas(i);
       for m=1:1:1440000
           if ((allData(i,m))>T && sameSpike==0)
               measuredNumSpikes(i,count)=measuredNumSpikes(i,count)+1;
               sameSpike=m;
           elseif (sameSpike~=0 && (allData(i,m))<T)
               sameSpike=0;
           end
       end
       count=count+1;
   end
   ks=2:15/2000:17;
   %figure();
   %plot(ks, measuredNumSpikes(i,:));
    
end

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


fun=@(b,x)(b(1)+b(2)./x+b(3)./(x.*x)+b(4)./(x.*x.*x));
beta0=[0 0 0 0];
beta=nlinfit(sigmas, closestK, fun, beta0);

plot(sigmas, fun(beta,sigmas));



