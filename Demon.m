%script for the large data set
resolution=0.075;
Dirx='\dataExample\';
fileAll=dir(strcat(Dirx,'\*.jpg'));
len=length(fileAll);

angleMain=nan(len,1);
period=nan(len,1);
type=nan(len,1);
bench=nan(len,1);
flag=nan(len,1);
cN=nan(len,1);

for i=1:len
im=imread(strcat(Dirx,fileAll(i).name));
[angleMain(i),tperiod,type(i),bench(i)] = rowFeatureIdentification_paddock_V4(im,resolution);
period(i)=tperiod(1);
end
