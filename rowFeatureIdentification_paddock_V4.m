function [angleMain,period,tType,bench] = rowFeatureIdentification_paddock_V4(im,resolution,mask)
%This estimate the periodical features in agricultural areas
%   INPUTS:
%       -im:     the high resolution optical image of the reasearch area
%       -mask:   a mask image with the selected paddock being 1 and the
%                background being 0. The mask should have the same size with im
%       -resolution: the spatial resolution of im
%   OUTPUTs:
%       -angleMain: the azimuth anlge of rows from North
%       -period: the periods of row features, with the first being that of
%       the domaint component 
%       -tType: the type of tillage for bare soil
%           1: Sinusoidal
%           2: Sinusoidal Bench
%           3: Bench
%       -bench: the length of bench
%Ref: Zhu, L., Walker, J.P., Rüdiger, C., & Xiao, P. (2020). Identification of agricultural row features using optical data for scattering and reflectance modelling over periodic surfaces. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 13, 1729-1739. 
%Liujun.zhu@monash.edu Dec. 23 2020


% Data preprocessing
if nargin==3
[ im,~] = paddockExtract( mask,im ); %get the maximum inscribed recontagle in the selected paddock
end
im=rgb2gray(uint8(im));%RGB 2 gray
im=im-mean(im(:));

% Fourier transformation 
F=fft2(double(im));%FFT
F=fftshift(F);%shift the 0 frequency to the center
Famp=abs(F);% Amplitude in frequency domain

% find the main orientation
[EnergyInAzimuth,maskMain,angleMain,flag,Fprofile,maxAmp]=lineScan(Famp);% scan in azimuth from 0 to 180

if flag==1% for paddocks with row features
    [damLoc,damH]=profileSegmentation(Fprofile,maxAmp);% scan in the main direction
    len=length(damLoc);
    [fre,angleMain,period,tType,bench]=PeriodicFeaturesExtraction(damH,F);
    period=period*resolution;
    bench=bench*resolution;
    angleMain=angleMain(1);
else
    angleMain=nan;
    period=nan;
    tType=nan;
    bench=nan;
end

end

function [damLoc,damH]=profileSegmentation(Fprofile,maxAmp)

tempF=Fprofile;
tempF(end-3:end)=0;%remove components influenced by the DC component
[pks,locs,dams1]=findMainPeaks(tempF,3);% merge peaks with a distance less than 5
% spring method
[maxloc]=find(Fprofile==maxAmp);
maxID=find(locs==maxloc);
dams2=zeros(length(locs),1);
for i=1:length(pks)
    if i<=maxID 
        dams2(i)=(max(pks(1:i))==pks(i));% for the left part
    else
        dams2(i)=(max(pks(i:end))==pks(i));% for the right part
    end
end
damLoc=locs((dams1+dams2)==2);
damH=tempF(damLoc);
[damH,Sid]=sort(damH,'descend');
if length(damH)>=5
damH=damH(1:5);%only the three largest components were considered for signal reconstruction
damLoc=damLoc(Sid(1:5));
end
end
function [EnergyInAzimuth,mask,angle,flag,Fprofile,maxAmp]=lineScan(Famp)
FampO=Famp;
[m,n]=size(Famp);
%Get a band pass filter to scan the circle area
%parts
yloc=repmat(1:m,n,1)-0.5;
xloc=repmat(1:n,m,1)-0.5;
yloc=yloc';
cX=n/2+0.5;% the location of 0 frequency
cY=m/2+0.5;
circleMask1=sqrt((xloc-cX).^2+(yloc-cY).^2)<(min(m,n)/2);
circleMask2=sqrt((xloc-cX).^2+(yloc-cY).^2)>0.01*(min(m,n)/2);
circleMask=(circleMask1+circleMask2)==2;
Famp(~circleMask)=nan;

% scan in the azimuth from 0(north) to 180 deg
EnergyInAzimuth=nan(360,1);

for i=1:360
    if i~=360
        k=tan(deg2rad(i/2+90));
        d=abs(k*xloc-yloc-(k*cX-cY))./sqrt(1+k^2);
        mask=d<1;% the scanning buffer 
        tempFamp=Famp(mask);
        EnergyInAzimuth(i+1)=nanmean(tempFamp);
    else
        mask=zeros(m,n);
        mask(:,floor(n/2)+1)=1;
        mask=mask==1;
        tempFamp=Famp(mask);
        EnergyInAzimuth(1)=nanmean(tempFamp);
    end
end

% Find the main angle and determine whether the isotropic features exist
angle=find(EnergyInAzimuth==max(EnergyInAzimuth))/2;%find the main direction
angle=angle(1)+0.00001;
[pks,locs,dams1]=findMainPeaks(EnergyInAzimuth,5);
flag=sum(EnergyInAzimuth(locs(dams1==1))>max(EnergyInAzimuth)*0.9048)<4;%-1 dB threshold C1

%scan the main direction and get the profile
k1=tan(deg2rad(angle+90));
d=abs(k1*xloc-yloc-(k1*cX-cY))./sqrt(1+k1^2);
mask=d<8;
FampO(~mask)=0;
flagD=((angle>=45)+(angle<=135))==2;
if ~flagD
    Fprofile=zeros(m,1);
    for i=1:m           
        Fprofile(i)=nanmax(FampO(i,:));
    end
else
    Fprofile=zeros(n,1);
    for i=1:n
        Fprofile(i)=nanmax(FampO(:,i));
    end
end
len=length(Fprofile);
Fprofile=Fprofile(1:floor(len/2)+1);
len=length(Fprofile);
tFprofile=Fprofile;
tFprofile(len-2:len)=0;
[maxAmp,maxIndex]=max(tFprofile);

flag=(flag+(maxIndex<=len-5))==2;
if flag==0
    angle=nan;
end
end
function[pks,locs,dams1]=findMainPeaks(profileRaw,T)

[pks,locs] = findpeaks(profileRaw);
% merge peaks with a distance less than T
temppks=pks;
temploc=locs;
dams1=zeros(length(locs),1);
while 1
    tempMaxLoc=find(temppks==nanmax(temppks));
    tempFlag=abs(locs-temploc(tempMaxLoc))<=T;
    temploc(tempFlag)=nan;
    temppks(tempFlag)=nan;
    dams1(tempMaxLoc)=1;
    if sum(~isnan(temploc))==0
        break;
    end
end
dams1=dams1==1;
end
function[fre,angleMain,period,tType,bench]=PeriodicFeaturesExtraction(damH,F)

tType=sum(damH/max(abs(F(:)))>0.0913);%C3
if tType~=0
    [m,n]=size(F);
    cm=floor(m/2)+1;
    cn=floor(n/2)+1;
    len=length(damH);
    fre=nan(len,1);
    angleMain=nan(len,1);
    phase=nan(len,1);
    period=nan(len,1);
    for i=1:len
        [pm,pn]=find(abs(F)==damH(i));
        fre(i)=sqrt((pm(1)-cm).^2+(pn(1)-cn).^2);
        angleMain(i)=90-rad2deg(atan((cn-pn(1))/(cm-pm(1))));% main orientation update
        phase(i)=angle(F(pm(1),pn(1)));
        if ((angleMain(i)<1)+(angleMain(i)>179))==1
           period(i)=abs(n/(cn-pn(1))*sin(deg2rad(90-angleMain(i)))); %Eq. 4
        else
           period(i)=abs(m/(cm-pm(1))*cos(deg2rad(90-angleMain(i)))); %Eq. 5
        end
    end
    
    % bench estimation
    x=0:0.1:1000;
    y=0;
    phase=phase-phase(1);
    for i=1:len
    y=y+damH(i)*sin(x*2*pi/period(i)+phase(i));
    end
    
    len=length(x(x<period(1)));
    [~,id]=min(y(1:len));
    x=x(id+10:id+len+20);
    y=y(id+10:id+len+20);
    y=y-min(y);
    [pk1,loc1]=findpeaks(y);
    index=pk1>mean(y);
    pk1=pk1(index);
    T=mean(pk1)/2;
    y1=(y-T)>0;
    y2=[y1(end),y1(1:end-1)];
    y1=abs(y1-y2);
    [loc]=find(y1==1);
    bench=0;
    for i=1:length(loc)-1
        bench=max(bench,loc(i+1)-loc(i));
    end
    if bench/len<0.5
        bench=(1-bench/len).*period(1);
    else
        bench=bench.*period(1)/len;
    end
else
    fre=nan;
    angleMain=nan;
    period=nan;
    tType=nan;
    bench=nan;
end

end
