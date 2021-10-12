%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to calculate the statistic results for wall-bounded   %
%   closed turbulent flows, calculated by Channel3d. Channel3d can be     %
%   freely downloaded from https://github.com/GongZheng-Justin/Channel3d  %
%                                                                         %
% There are 12 input parameters below                                     %
%   * xlx: Domain length in x-dir                                         %
%   * zlz: Domain lenght in z-dir                                         %
%   * nxc: Grid number in x-dir (nxc=nxp-1)                               %                                
%   * nzc: Gird number in z-dir (nzc=nzp-1)                               %
%   * xnu: Fluid kinematic viscosity                                      %
%   * iTimeSet: Starting time for statistics calculation                  %
%   * jSpecSet: The first y-index for energy spectra calculation          %     
%   * jSpecEnd: The last y-index for energy spectra calculation           %
%   * jSpecInc: The y-index interval for energy spectra                   %
%   * dir_statIn: The folder to store the original/raw statistic data     %
%   * dir_statOut:The folder to dump the final statistic results          %
%   * yMesh_str:  Ymesh file name                                         %
%                                                                         %
% The final statistic results will be dumped in the folder `dir_statOut`  %
%   commonly including the following files:                               %
%   * Profile.txt: Velocity/vorticity/pressure profiles etc.              % 
%   * kBudget.txt: Budget for TKE                                         %
%   * uuBudget.txt:Budgtt for <u'u'>                                      %
%   * vvBudget.txt:Budgtt for <v'v'>                                      %
%   * wwBudget.txt:Budgtt for <w'w'>                                      %
%   * uvBudget.txt:Budgtt for <u'v'>                                      %
%   * specx_u.txt: 1D u-velocity spectra in x-dir                         %
%   * specx_v.txt: 1D v-velocity spectra in x-dir                         %
%   * specx_w.txt: 1D w-velocity spectra in x-dir                         %
%   * specx_p.txt: 1D p-velocity spectra in x-dir                         %
%   * specz_u.txt: 1D u-velocity spectra in z-dir                         %
%   * specz_v.txt: 1D v-velocity spectra in z-dir                         %
%   * specz_w.txt: 1D w-velocity spectra in z-dir                         %
%   * specz_p.txt: 1D p-velocity spectra in z-dir                         %
%                                                                         %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-10-12                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

clear;clc;

xlx=4*pi;
zlz=2*pi;
nxc=384;
nzc=256;
xnu=2.3310E-4;
iTimeSet=6000;

IsUxConst=0;  % 1=True; 0=False
BodyForceX=0; % If IsUxConst=1, use BodyForceX to calculate u_tau, i.e. u_tau=sqrt(BodyForceX*height) 

jSpecSet = 5;
jSpecEnd = 95;
jSpecInc = 5;
dir_statIn = '../CFD/Results/';
dir_statOut='../StatOut/';
yMesh_str=[dir_statIn,'yMeshForCha180_2_.txt'];

%%========== Normally no need to change anything below ==========
%% read the mesh and calculate relative quantities
if ( exist(dir_statOut,'dir') ==false )
  mkdir(dir_statOut);
  fprintf( '%s\n\n',['crete directory: ',dir_statOut,'  sucessfully'] );
end

yp=load(yMesh_str);
nyp=size(yp,1); nyc=nyp-1; yp=yp(:,2);
yc=0.5*(yp(1:end-1)+yp(2:end));
height=0.5*yp(end);

dyp=zeros(nyc,1);
dyc=zeros(nyp,1);
for k=1:nyc
  dyp(k)=yp(k+1)-yp(k);
end
for k=2:nyc
  dyc(k)=yc(k)-yc(k-1);
end
dyc(1)=dyp(1);
dyc(nyp)=dyp(nyc);

am2c=zeros(nyc,1);
ap2c=zeros(nyc,1);
for k=1:nyc
  am2c(k)= 1/dyp(k)/dyc(k);
  ap2c(k)= 1/dyp(k)/dyc(k+1);
end
am2c(1)= 4.0/dyc(1)/( dyc(1)+2.0*dyc(2) );
ap2c(1)= 4.0/dyc(2)/( dyc(1)+2.0*dyc(2) );
am2c(nyc)= 4.0/dyc(nyc)/( dyc(nyp)+2.0*dyc(nyc) );
ap2c(nyc)= 4.0/dyc(nyp)/( dyc(nyp)+2.0*dyc(nyc) );
ac2c= -(am2c+ap2c);

am2p=zeros(nyp,1);
ap2p=zeros(nyp,1);
for k=2:nyc
  am2p(k)= 1.0/dyc(k)/dyp(k-1);        
  ap2p(k)= 1.0/dyc(k)/dyp(k);
end
am2p(1)=2.0/(dyp(1)*(dyp(1)+dyp(2)));
ap2p(1)=2.0/(dyp(2)*(dyp(1)+dyp(2)));
am2p(nyp)=2.0/(dyp(nyc-1)*(dyp(nyc)+dyp(nyc-1)));
ap2p(nyp)=2.0/(dyp(nyc)*(dyp(nyc)+dyp(nyc-1)));
ac2p= -(am2p+ap2p);

%% calculate averaged data
dir_output=dir(fullfile(dir_statIn, 'stats0*') );
file_names={dir_output.name};
file_num=length(file_names);

for k=1:file_num
  datapatht=cell2mat(file_names(k));
  if(str2double(datapatht(6:15))>iTimeSet); 
    datapath = [dir_statIn,datapatht];
    break;
  end
end

fid=fopen(datapath,'r');
line1 = 0;

while (feof(fid)==0)
  str=strtrim( fgets(fid) );
  dlt=sscanf(str,'%f');
  line1=line1 +1;
  if(isempty(dlt)==0 && length(dlt)>2)
    break;
  end
end
real_num=length(dlt);
frewind(fid);
line2 = 0;

while (feof(fid)==0)
  fgets(fid);
  line2= line2 +1;
end
fclose(fid);
file_len=line2-line1+1;
if(file_len ~=nyp) 
  error('nyp wrong');
end

file_ave=0;
for k=1:file_num
  datapath=cell2mat(file_names(k));
  if(str2double(datapath(6:15))>iTimeSet); 
    file_ave=file_ave+1;
  end
end

file_ave=0;
prgrad=zeros(file_ave,1);
dataE=zeros(nyp,real_num);
data_emerget=zeros(nyp,real_num);
for k=1:file_num
  datapath=cell2mat(file_names(k));
  if(str2double(datapath(6:15))<=iTimeSet);continue;end
  file_ave=file_ave+1;  
  datapath = [dir_statIn, cell2mat(file_names(k)) ];
  fid=fopen(datapath,'r');
  for kt=1:line1-1
    str=strtrim( fgets(fid) );
    prgradt=sscanf(str,'time averaged pressure gradient is:%f');
    if(isempty(prgradt)==0)
      prgrad(file_ave)=prgradt;
    end      
  end    
  
  for kt=line1:line2
    idl=kt-line1+1;
    str=strtrim( fgets(fid) );
    dlt=sscanf(str,'%f');
    data_emerget(idl,:)=dlt';
  end
  fclose(fid);
  dataE=dataE+ data_emerget;
  disp( ['read:   ',datapath,'  sucessfully'] );
end
dataE=dataE/file_ave;

% dudyp, dudyc
dudyp=zeros(nyp,1);
for k=2:nyc
  dudyp(k)=(dataE(k,1)-dataE(k-1,1))/dyc(k);
end
k=1;
dudyp(k)= 2.0*dataE(1,1)/dyc(k);
k=nyp;
dudyp(k)=-2.0*dataE(nyc,1)/dyc(k);
dudyc=0.5*(dudyp(1:end-1)+dudyp(2:end));

% dwdyp, dwdyc
dwdyp=zeros(nyp,1);
for k=2:nyc
  dwdyp(k)=(dataE(k,3)-dataE(k-1,3))/dyc(k);
end
k=1;
dwdyp(k)= 2.0*dataE(1,3)/dyc(k);
k=nyp;
dwdyp(k)=-2.0*dataE(nyc,3)/dyc(k);
dwdyc=0.5*(dwdyp(1:end-1)+dwdyp(2:end));

if(IsUxConst==0) 
  utau1=mean(sqrt(prgrad*height));
  utau2=0.5*(sqrt(xnu*dudyp(1))+sqrt(-xnu*dudyp(nyp)));
  utaufinal=0.5*(utau1+utau2);
else
  utaufinal=sqrt(BodyForceX*height);
end
disp(['utaufinal=',num2str(utaufinal)]);
vorMag=utaufinal^2/xnu;
budMag=utaufinal^4/xnu;
utauSqr=utaufinal^2;
nych=nyc/2;

%% profile first
dataPrf=zeros(nych,27);
for k=1:nych
  ksc=nyc+1-k; ksp=nyp+1-k;
  dataPrf(k,1)=yc(k)/height;
  dataPrf(k,2)=yc(k)*utaufinal/xnu;
  
  % ux, symmetry
  dataPrf(k,3)=0.5*(dataE(k,1)+dataE(ksc,1))/utaufinal;
  
  % uy, antisymmetry
  dataPrf(k,4)=0.25*(dataE(k,2)+dataE(k+1,2)-dataE(ksc,2)-dataE(ksp,2))/utaufinal;
  
  % uz, unknown
  dataPrf(k,5)=0.5*(dataE(k,3)+dataE(ksc,3))/utaufinal;  

  % pr, symmetry
  dataPrf(k,6)=0.5*(dataE(k,4)+dataE(ksc,4))/utauSqr; 
  
  % uRms
  dataPrf(k,7)=0.5*(sqrt(dataE(k,5)-dataE(k,1)^2)+sqrt(dataE(ksc,5)-dataE(ksc,1)^2));
  dataPrf(k,7)=dataPrf(k,7)/utaufinal;
  
  % vRms
  RealT1=sqrt(dataE(k,  6)-dataE(k,  2)^2);
  RealT2=sqrt(dataE(k+1,6)-dataE(k+1,2)^2);
  RealT3=sqrt(dataE(ksc,6)-dataE(ksc,2)^2);
  RealT4=sqrt(dataE(ksp,6)-dataE(ksp,2)^2);
  dataPrf(k,8)= 0.25*(RealT1+RealT2+RealT3+RealT4)/utaufinal;

  % wRms
  dataPrf(k,9)=0.5*(sqrt(dataE(k,7)-dataE(k,3)^2)+sqrt(dataE(ksc,7)-dataE(ksc,3)^2));
  dataPrf(k,9)=dataPrf(k,9)/utaufinal;

  % u'v'
  RealT1= dataE(k,  9)-dataE(k,  1)*(dataE(k,  2)+dataE(k+1,2))*0.5;
  RealT2= dataE(ksc,9)-dataE(ksc,1)*(dataE(ksc,2)+dataE(ksp,2))*0.5;
  dataPrf(k,10)= 0.5*(RealT1-RealT2)/utauSqr;

  % v'w'
  RealT1= dataE(k,  10)-dataE(k,  3)*(dataE(k,  2)+dataE(k+1,2))*0.5;
  RealT2= dataE(ksc,10)-dataE(ksc,3)*(dataE(ksc,2)+dataE(ksp,2))*0.5;
  dataPrf(k,11)= 0.5*(RealT1-RealT2)/utauSqr;

  % u'w'
  RealT1= dataE(k,  11)-dataE(k,  3)*dataE(k,  1);
  RealT2= dataE(ksc,11)-dataE(ksc,3)*dataE(ksc,1);
  dataPrf(k,12)= 0.5*(RealT1+RealT2)/utauSqr;

  % pRms
  dataPrf(k,13)=0.5*(sqrt(dataE(k,8)-dataE(k,4)^2)+sqrt(dataE(ksc,8)-dataE(ksc,4)^2));
  dataPrf(k,13)=dataPrf(k,13)/utaufinal^2; 

  % u'p'
  RealT1= dataE(k,  12)-dataE(k,  1)*dataE(k,  4);
  RealT2= dataE(ksc,12)-dataE(ksc,1)*dataE(ksc,4);
  dataPrf(k,14)=0.5*(RealT1+RealT2)/utaufinal^3;  

  % v'p'
  RealT1= dataE(k,  13)-0.5*(dataE(k,  2)+dataE(k+1,2))*dataE(k,  4);
  RealT2= dataE(ksc,13)-0.5*(dataE(ksc,2)+dataE(ksp,2))*dataE(ksc,4);
  dataPrf(k,15)=0.5*(RealT1-RealT2)/utaufinal^3;

  % w'p'
  RealT1= dataE(k,  14)-dataE(k,  3)*dataE(k,  4);
  RealT2= dataE(ksc,14)-dataE(ksc,3)*dataE(ksc,4);
  dataPrf(k,16)=0.5*(RealT1+RealT2)/utaufinal^3;

  % S(u)
  RealT1= dataE(k,  19)-3*dataE(k,  1)*dataE(k,  5)+2.0*dataE(k,  1)^3;
  RealT2= dataE(ksc,19)-3*dataE(ksc,1)*dataE(ksc,5)+2.0*dataE(ksc,1)^3;
  RealT1= RealT1/(dataE(k,  5)-dataE(k,  1)^2)^(3/2);
  RealT2= RealT2/(dataE(ksc,5)-dataE(ksc,1)^2)^(3/2);
  dataPrf(k,17)=0.5*(RealT1+RealT2);

  % S(v)
  RealT1= 0.5*(dataE(k,  16)-3*dataE(k,  6)*dataE(k,  2)+2.0*dataE(k,  2)^3+ dataE(k+1,16)-3*dataE(k+1,6)*dataE(k+1,2)+2.0*dataE(k+1,2)^3);
  RealT2= 0.5*(dataE(ksc,16)-3*dataE(ksc,6)*dataE(ksc,2)+2.0*dataE(ksc,2)^3+ dataE(ksp,16)-3*dataE(ksp,6)*dataE(ksp,2)+2.0*dataE(ksp,2)^3);
  RealT3= 0.5*((dataE(k,  6)-dataE(k,  2)^2)^(3/2) + (dataE(k+1,6)-dataE(k+1,2)^2)^(3/2));
  RealT4= 0.5*((dataE(ksc,6)-dataE(ksc,2)^2)^(3/2) + (dataE(ksp,6)-dataE(ksp,2)^2)^(3/2));
  dataPrf(k,18)=0.5*(RealT1/RealT3-RealT2/RealT4);

  % S(w)
  RealT1= dataE(k,  20)-3*dataE(k,  3)*dataE(k,  7)+2.0*dataE(k,  3)^3;
  RealT2= dataE(ksc,20)-3*dataE(ksc,3)*dataE(ksc,7)+2.0*dataE(ksc,3)^3;
  RealT1= RealT1/(dataE(k,  7)-dataE(k,  3)^2)^(3/2);
  RealT2= RealT2/(dataE(ksc,7)-dataE(ksc,3)^2)^(3/2);
  dataPrf(k,19)=0.5*(RealT1+RealT2);

  % F(u)
  RealT1= dataE(k,  21)-4*dataE(k,  1)*dataE(k,  19)+6*dataE(k,  5)*dataE(k,  1)^2-3*dataE(k,  1)^4;
  RealT2= dataE(ksc,21)-4*dataE(ksc,1)*dataE(ksc,19)+6*dataE(ksc,5)*dataE(ksc,1)^2-3*dataE(ksc,1)^4;
  RealT1= RealT1/(dataE(k,  5)-dataE(k,  1)^2)^2;
  RealT2= RealT2/(dataE(ksc,5)-dataE(ksc,1)^2)^2;
  dataPrf(k,20)=0.5*(RealT1+RealT2);

  % F(v)
  RealT1= 0.5*(dataE(k,  22)-4*dataE(k,  2)*dataE(k,  16)+6*dataE(k,  6)*dataE(k,  2)^2-3*dataE(k,  2)^4 +...
               dataE(k+1,22)-4*dataE(k+1,2)*dataE(k+1,16)+6*dataE(k+1,6)*dataE(k+1,2)^2-3*dataE(k+1,2)^4);
  RealT2= 0.5*(dataE(ksc,22)-4*dataE(ksc,2)*dataE(ksc,16)+6*dataE(ksc,6)*dataE(ksc,2)^2-3*dataE(ksc,2)^4 +...
               dataE(ksp,22)-4*dataE(ksp,2)*dataE(ksp,16)+6*dataE(ksp,6)*dataE(ksp,2)^2-3*dataE(ksp,2)^4);
  RealT3= 0.5*((dataE(k,  6)-dataE(k,  2)^2)^2 + (dataE(k+1,6)-dataE(k+1,2)^2)^2);
  RealT4= 0.5*((dataE(ksc,6)-dataE(ksc,2)^2)^2 + (dataE(ksp,6)-dataE(ksp,2)^2)^2);
  dataPrf(k,21)=0.5*(RealT1/RealT3+RealT2/RealT4);

  % F(w)
  RealT1= dataE(k,  23)-4*dataE(k,  3)*dataE(k,  20)+6*dataE(k,  7)*dataE(k,  3)^2-3*dataE(k,  3)^4;
  RealT2= dataE(ksc,23)-4*dataE(ksc,3)*dataE(ksc,20)+6*dataE(ksc,7)*dataE(ksc,3)^2-3*dataE(ksc,3)^4;
  RealT1= RealT1/(dataE(k,  7)-dataE(k,  3)^2)^2;
  RealT2= RealT2/(dataE(ksc,7)-dataE(ksc,3)^2)^2;
  dataPrf(k,22)=0.5*(RealT1+RealT2);

  % dudy
  dataPrf(k,23)=0.5*(dudyc(k)-dudyc(ksc))/vorMag;

  % dwdy
  dataPrf(k,24)=0.5*(dwdyc(k)-dwdyc(ksc))/vorMag;

  % wxRms
  RealT1=sqrt(dataE(k,  33)-dwdyp(k  )^2);
  RealT2=sqrt(dataE(k+1,33)-dwdyp(k+1)^2);
  RealT3=sqrt(dataE(ksc,33)-dwdyp(ksc)^2);
  RealT4=sqrt(dataE(ksp,33)-dwdyp(ksp)^2);
  dataPrf(k,25)= 0.25*(RealT1+RealT2+RealT3+RealT4)/vorMag;

  % wyRms
  dataPrf(k,26)=0.5*(sqrt(dataE(k,34))+sqrt(dataE(ksc,34)))/vorMag;

  % wzRms
  RealT1=sqrt(dataE(k,  35)-dudyp(k  )^2);
  RealT2=sqrt(dataE(k+1,35)-dudyp(k+1)^2);
  RealT3=sqrt(dataE(ksc,35)-dudyp(ksc)^2);
  RealT4=sqrt(dataE(ksp,35)-dudyp(ksp)^2);
  dataPrf(k,27)= 0.25*(RealT1+RealT2+RealT3+RealT4)/vorMag;
end
myformat=[repmat('%24.15E',1,27),'\n'];
fid=fopen([dir_statOut,'Profile.txt'],'wt');
strout=' ybar yplus u v w p uRms vRms wRms u''v'' v''w'' u''w'' pRms u''p'' v''p'' w''p'' S(u) S(v) S(w) F(u) F(v) F(w) dudy dwdy wxRms wyRms wzRms';
fprintf(fid,'%s\n',strout);
for k=1:nych
  fprintf(fid,myformat,dataPrf(k,:));
end
fclose(fid);

%% uu budget
uuBud=zeros(nych,8);

% TranTerm
TranTermp=zeros(nyp,1);
TranTermc=zeros(nyc,1);
for k=1:nyc
  TranTermc(k)=dataE(k,15)-2.0*dataE(k,1)*dataE(k,9)-0.5*(dataE(k,2)+dataE(k+1,2))*dataE(k,5)+(dataE(k,2)+dataE(k+1,2))*dataE(k,1)^2;
end
for k=2:nyc
  TranTermp(k)=(TranTermc(k)-TranTermc(k-1))/dyc(k);
end
k=1;
TranTermp(k)= 2.0*TranTermc(1)/dyc(k);
k=nyp;
TranTermp(k)=-2.0*TranTermc(nyc)/dyc(k);
TranTermc=0.5*(TranTermp(1:end-1)+TranTermp(2:end));

% Viscous Transport
VisTranc=zeros(nyc,1);
VisTranp=zeros(nyc+2,1);
for k=1:nyc
  kt=k+1;
  VisTranp(kt)=dataE(k,5)-dataE(k,1)^2;
end
VisTranp(1)    = -VisTranp(2);
VisTranp(nyc+2)= -VisTranp(nyc+1);
for k=1:nyc
  VisTranc(k)= ap2c(k)*VisTranp(k+2) +ac2c(k)*VisTranp(k+1) + am2c(k)*VisTranp(k);
end

% dudyy
dudyyc=zeros(nyc,1);
dudyyp=zeros(nyc+2,1);
for k=1:nyc
  kt=k+1;
  dudyyp(kt)=dataE(k,1);
end
dudyyp(1)    = -dudyyp(2);
dudyyp(nyc+2)= -dudyyp(nyc+1);
for k=1:nyc
  dudyyc(k)= ap2c(k)*dudyyp(k+2) +ac2c(k)*dudyyp(k+1) + am2c(k)*dudyyp(k);
end

for k=1:nych
  ksc=nyc+1-k;ksp=nyp+1-k;
  uuBud(k,1)=yc(k)/height;
  uuBud(k,2)=yc(k)*utaufinal/xnu;
  
  % Production          (P)
  RealT1= dataE(k,  9)-dataE(k,  1)*(dataE(k,  2)+dataE(k+1,2))*0.5;
  RealT2= dataE(ksc,9)-dataE(ksc,1)*(dataE(ksc,2)+dataE(ksp,2))*0.5;
  RealT1= -2.0*RealT1*dudyc(k);  RealT2=-2.0*RealT2*dudyc(ksc);
  uuBud(k,3)= 0.5*(RealT1+RealT2)/budMag;

  % Turbulent Transport (T)
  uuBud(k,4)= -0.5*(TranTermc(k)+TranTermc(ksc))/budMag;

  % Viscous Transport   (D)
  uuBud(k,5)= 0.5*(VisTranc(k)+VisTranc(ksc))*xnu/budMag;

  % Pressure Strain     (pi)
  RealT1=2.0*dataE(k,  24);
  RealT2=2.0*dataE(ksc,24);
  uuBud(k,6)= 0.5*(RealT1+RealT2)/budMag;

  % Pressure Transport  (phi)
  uuBud(k,7)= 0;

  % Viscous Dissipation (epsilon)
  RealT1= 0.5*VisTranc(k  )-dataE(k,  28)+dataE(k,  1)*dudyyc(k  );
  RealT2= 0.5*VisTranc(ksc)-dataE(ksc,28)+dataE(ksc,1)*dudyyc(ksc);
  uuBud(k,8)= 0.5*(RealT1+RealT2)*2.0*xnu/budMag;
end
myformat=[repmat('%24.15E',1,8),'\n'];
fid=fopen([dir_statOut,'uuBudget.txt'],'wt');
fprintf(fid,' ybar yplus P T D pi phi epsilon\n');
for k=1:nych
  fprintf(fid,myformat,uuBud(k,:));
end
fclose(fid);

%% vv budget
vvBud=zeros(nych,8);

% TranTerm
for k=1:nyp
  TranTermp(k)=dataE(k,16)-3.0*dataE(k,2)*dataE(k,6)+2.0*dataE(k,2)^3;
end
for k=1:nyc
  TranTermc(k)= (TranTermp(k+1)-TranTermp(k))/dyp(k);
end

% Viscous Transport
VisTranT=zeros(nyp,1);
for k=1:nyp
  VisTranT(k)=dataE(k,6)-dataE(k,2)^2 ;
end
for k=2:nyc
  VisTranp(k)= am2p(k)*VisTranT(k-1) +ac2p(k)*VisTranT(k)+ ap2p(k)*VisTranT(k+1);
end
k=1;   VisTranp(k)= am2p(k)*VisTranT(1)     +ac2p(k)*VisTranT(2)  + ap2p(k)*VisTranT(3);
k=nyp; VisTranp(k)= am2p(k)*VisTranT(nyc-1) +ac2p(k)*VisTranT(nyc)+ ap2p(k)*VisTranT(nyp);
VisTranc= 0.5*(VisTranp(1:end-1)+VisTranp(2:end));

% PVTerm
PVTermp=zeros(nyp,1);
PVTermc=zeros(nyc,1);
for k=1:nyc
  PVTermc(k)=dataE(k,13)- 0.5*(dataE(k,2)+dataE(k+1,2))*dataE(k,4);
end
for k=2:nyc
  PVTermp(k)=(PVTermc(k)-PVTermc(k-1))/dyc(k);
end
k=1;
PVTermp(k)= 2.0*PVTermc(1,1)/dyc(k);
k=nyp;
PVTermp(k)=-2.0*PVTermc(nyc,1)/dyc(k);
PVTermc=0.5*(PVTermp(1:end-1)+PVTermp(2:end));

for k=1:nych
  ksc=nyc+1-k;ksp=nyp+1-k;
  vvBud(k,1)=yc(k)/height;
  vvBud(k,2)=yc(k)*utaufinal/xnu;

  % Production          (P)
  vvBud(k,3)=0; 

  % Turbulent Transport (T)
  vvBud(k,4)=  -0.5*(TranTermc(k)+TranTermc(ksc))/budMag;

  % Viscous Transport   (D)
  vvBud(k,5)= 0.5*(VisTranc(k)+VisTranc(ksc))*xnu/budMag;

  % Pressure Strain     (pi)
  RealT1=2.0*dataE(k,  25);
  RealT2=2.0*dataE(ksc,25);
  vvBud(k,6)= 0.5*(RealT1+RealT2)/budMag;

  % Pressure Strain     (phi)
  vvBud(k,7)= -2.0*(0.5*(PVTermc(k)+PVTermc(ksc)))/budMag;

  % Viscous Dissipation (epsilon)
  RealT1= 0.5*VisTranp(k  )-dataE(k,  29);
  RealT2= 0.5*VisTranp(k+1)-dataE(k+1,29);  
  RealT3= 0.5*VisTranp(ksc)-dataE(ksc,29);  
  RealT4= 0.5*VisTranp(ksp)-dataE(ksp,29);
  vvBud(k,8)= 0.25*(RealT1+RealT2+RealT3+RealT4)*2.0*xnu/budMag;
end
myformat=[repmat('%24.15E',1,8),'\n'];
fid=fopen([dir_statOut,'vvBudget.txt'],'wt');
fprintf(fid,' ybar yplus P T D pi phi epsilon\n');
for k=1:nych
  fprintf(fid,myformat,vvBud(k,:));
end
fclose(fid);

%% ww budget
wwBud=zeros(nych,8);

% TranTerm
TranTermp=zeros(nyp,1);
TranTermc=zeros(nyc,1);
for k=1:nyc
  TranTermc(k)=dataE(k,17)-2.0*dataE(k,3)*dataE(k,10)-0.5*(dataE(k,2)+dataE(k+1,2))*dataE(k,7)+(dataE(k,2)+dataE(k+1,2))*dataE(k,3)^2;
end
for k=2:nyc
  TranTermp(k)=(TranTermc(k)-TranTermc(k-1))/dyc(k);
end
k=1;
TranTermp(k)= 2.0*TranTermc(1)/dyc(k);
k=nyp;
TranTermp(k)=-2.0*TranTermc(nyc)/dyc(k);
TranTermc=0.5*(TranTermp(1:end-1)+TranTermp(2:end));

% Viscous Transport
VisTranc=zeros(nyc,1);
VisTranp=zeros(nyc+2,1);
for k=1:nyc
  kt=k+1;
  VisTranp(kt)=dataE(k,7)-dataE(k,3)^2;
end
VisTranp(1)    = -VisTranp(2);
VisTranp(nyc+2)= -VisTranp(nyc+1);
for k=1:nyc
  VisTranc(k)= ap2c(k)*VisTranp(k+2) +ac2c(k)*VisTranp(k+1) + am2c(k)*VisTranp(k);
end

% dwdyy
dwdyyc=zeros(nyc,1);
dwdyyp=zeros(nyc+2,1);
for k=1:nyc
  kt=k+1;
  dwdyyp(kt)=dataE(k,3);
end
dwdyyp(1)    = -dwdyyp(2);
dwdyyp(nyc+2)= -dwdyyp(nyc+1);
for k=1:nyc
  dwdyyc(k)= ap2c(k)*dwdyyp(k+2) +ac2c(k)*dwdyyp(k+1) + am2c(k)*dwdyyp(k);
end

for k=1:nych
  ksc=nyc+1-k;ksp=nyp+1-k;
  wwBud(k,1)=yc(k)/height;
  wwBud(k,2)=yc(k)*utaufinal/xnu;

  % Production          (P)
  RealT1= dataE(k,  10)-dataE(k,  3)*(dataE(k,  2)+dataE(k+1,2))*0.5;
  RealT2= dataE(ksc,10)-dataE(ksc,3)*(dataE(ksc,2)+dataE(ksp,2))*0.5;
  RealT1= -2.0*RealT1*dwdyc(k);  RealT2=-2.0*RealT2*dwdyc(ksc);
  wwBud(k,3)= 0.5*(RealT1+RealT2)/budMag;

  % Turbulent Transport (T)
  wwBud(k,4)= -0.5*(TranTermc(k)+TranTermc(ksc))/budMag;

  % Viscous Transport   (D)
  wwBud(k,5)= 0.5*(VisTranc(k)+VisTranc(ksc))*xnu/budMag;

  % Pressure Strain     (pi)
  RealT1=2.0*dataE(k,  26);
  RealT2=2.0*dataE(ksc,26);
  wwBud(k,6)= 0.5*(RealT1+RealT2)/budMag;

  % Pressure Transport  (phi)
  wwBud(k,7)= 0;

  % Viscous Dissipation (epsilon)
  RealT1= 0.5*VisTranc(k  )-dataE(k,  30)+dataE(k,  3)*dwdyyc(k  );
  RealT2= 0.5*VisTranc(ksc)-dataE(ksc,30)+dataE(ksc,3)*dwdyyc(ksc);
  wwBud(k,8)= 0.5*(RealT1+RealT2)*2.0*xnu/budMag;
end
myformat=[repmat('%24.15E',1,8),'\n'];
fid=fopen([dir_statOut,'wwBudget.txt'],'wt');
fprintf(fid,' ybar yplus P T D pi phi epsilon\n');
for k=1:nych
  fprintf(fid,myformat,wwBud(k,:));
end
fclose(fid);

%% uv budget
uvBud=zeros(nych,8);

% TranTerm
TranTermp=zeros(nyp,1);
TranTermc=zeros(nyc,1);
for k=1:nyc
  TranTermc(k)=dataE(k,18)-2.0*(0.5*(dataE(k,2)+dataE(k+1,2)))*dataE(k,9)-dataE(k,1)*0.5*(dataE(k,6)+dataE(k+1,6))+2.0*dataE(k,1)*(0.5*(dataE(k,2)+dataE(k+1,2)))^2;
end
for k=2:nyc
  TranTermp(k)=(TranTermc(k)-TranTermc(k-1))/dyc(k);
end
k=1;
TranTermp(k)= 2.0*TranTermc(1)/dyc(k);
k=nyp;
TranTermp(k)=-2.0*TranTermc(nyc)/dyc(k);
TranTermc=0.5*(TranTermp(1:end-1)+TranTermp(2:end));

% Viscous Transport
VisTranc=zeros(nyc,1);
VisTranp=zeros(nyc+2,1);
for k=1:nyc
  kt=k+1;
  VisTranp(kt)=dataE(k,9)-dataE(k,1)*(dataE(k,2)+dataE(k+1,2))*0.5;
end
VisTranp(1)    = -VisTranp(2);
VisTranp(nyc+2)= -VisTranp(nyc+1);
for k=1:nyc
  VisTranc(k)= ap2c(k)*VisTranp(k+2) +ac2c(k)*VisTranp(k+1) + am2c(k)*VisTranp(k);
end

% PUTerm
PUTermp=zeros(nyp,1);
PUTermc=zeros(nyc,1);
for k=1:nyc
  PUTermc(k)=dataE(k,12)- dataE(k,1)*dataE(k,4);
end
for k=2:nyc
  PUTermp(k)=(PUTermc(k)-PUTermc(k-1))/dyc(k);
end
k=1;
PUTermp(k)= 2.0*PUTermc(1,1)/dyc(k);
k=nyp;
PUTermp(k)=-2.0*PUTermc(nyc,1)/dyc(k);
PUTermc=0.5*(PUTermp(1:end-1)+PUTermp(2:end));

for k=1:nych
  ksc=nyc+1-k;ksp=nyp+1-k;
  uvBud(k,1)=yc(k)/height;
  uvBud(k,2)=yc(k)*utaufinal/xnu;

  % Production          (P)
  RealT1= 0.5*(dataE(k,6)  -dataE(k,  2)^2 +dataE(k+1,6)-dataE(k+1,2)^2);
  RealT2= 0.5*(dataE(ksc,6)-dataE(ksc,2)^2 +dataE(ksp,6)-dataE(ksp,2)^2);
  RealT1= -RealT1*dudyc(k);  RealT2=-RealT2*dudyc(ksc);
  uvBud(k,3)= 0.5*(RealT1-RealT2)/budMag;

  % Turbulent Transport (T)
  uvBud(k,4)= -0.5*(TranTermc(k)-TranTermc(ksc))/budMag;

  % Viscous Transport   (D)
  uvBud(k,5)= 0.5*(VisTranc(k)-VisTranc(ksc))*xnu/budMag;

  % Pressure Strain     (pi)
  RealT1=0.5*(dataE(k,  27)+dataE(k+1,27))-dataE(k,  4)*dudyc(k);
  RealT2=0.5*(dataE(ksc,27)+dataE(ksp,27))-dataE(ksc,4)*dudyc(ksc);
  uvBud(k,6)= 0.5*(RealT1-RealT2)/budMag;

  % Pressure Strain     (phi)
  uvBud(k,7)= -0.5*(PUTermc(k)-PUTermc(ksc))/budMag;

  % Viscous Dissipation (epsilon)
  RealT1= 0.5*(dataE(k,31)-dataE(ksc,31));
  RealT2= 0.25*(dataE(k,32)+dataE(k+1,32)-dataE(ksc,32)-dataE(ksp,32));
  uvBud(k,8)= (RealT1+RealT2)*2*xnu/budMag;
end
myformat=[repmat('%24.15E',1,8),'\n'];
fid=fopen([dir_statOut,'uvBudget.txt'],'wt');
fprintf(fid,' ybar yplus P T D pi phi epsilon\n');
for k=1:nych
  fprintf(fid,myformat,uvBud(k,:));
end
fclose(fid);

%% k budget
kBud=zeros(nych,8);
kBud(:,3:8)=0.5*(uuBud(:,3:8)+vvBud(:,3:8)+wwBud(:,3:8));
for k=1:nych
  ksc=nyc+1-k;ksp=nyp+1-k;
  kBud(k,1)=yc(k)/height;
  kBud(k,2)=yc(k)*utaufinal/xnu;
  
  % Pressure Strain     (pi)
  kBud(k,6)=0;
end
myformat=[repmat('%24.15E',1,8),'\n'];
fid=fopen([dir_statOut,'kBudget.txt'],'wt');
fprintf(fid,' ybar yplus P T D pi phi epsilon\n');
for k=1:nych
  fprintf(fid,myformat,kBud(k,:));
end
fclose(fid);

%% Energy Spectra in x-dir
njSpec=(jSpecEnd-jSpecSet)/jSpecInc+1;
yplus_spec_center=yc(jSpecSet:jSpecInc:jSpecEnd)*utaufinal/xnu;
yplus_spec_center=yplus_spec_center';
yplus_spec_node=yp(jSpecSet:jSpecInc:jSpecEnd)*utaufinal/xnu;
yplus_spec_node=yplus_spec_node';

nxh=nxc/2; nxhp=nxh+1;
SpectraX=zeros(nxhp,njSpec*4);
dir_output=dir(fullfile(dir_statIn,'SpecX*'));
file_names={dir_output.name};
file_num=length(file_names);

if(file_num>0)
file_ave=0;
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))>iTimeSet)
   file_ave=file_ave+1;   
  end   
end
icolumn=zeros(file_ave,3);
columnCount=0;
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))<=iTimeSet);continue;end;
  icolumnt = str2double(datapath(17:20));
  CountFlag=1;
  if(columnCount>0)
    for m=1:columnCount
      if(icolumnt==icolumn(m,1))
        CountFlag=0;
        break
      end
    end
  end
  if(CountFlag==1)
    columnCount=columnCount+1;
    datapath = [dir_statIn,cell2mat(file_names(k)) ];
    SpecData = load(datapath);
    if(size(SpecData,1)~=nxhp)
      error('SpecDataX size wrong');
    end
    icolumn(columnCount,1)=icolumnt;
    icolumn(columnCount,3)=size(SpecData,2);
  end
end
if(mod(file_ave,columnCount)~=0) 
  error('SpecDataX columnCount wrong');
end
ntimeX=file_ave/columnCount;
icolumn=icolumn(1:columnCount,1:3);
icolumn=sortrows(icolumn);
sumColumn=0;
for k=1:columnCount
  sumColumn   = sumColumn+icolumn(k,3);
  icolumn(k,3)= sumColumn;
end
icolumn(1,2)=1;
if(columnCount>1)
  for k=2:columnCount
    icolumn(k,2)=icolumn(k-1,3)+1;
  end
end
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))<=iTimeSet);continue;end;
  icolumnt = str2double(datapath(17:20))+1;
  if(icolumn(icolumnt,1)~=icolumnt-1)
    error('SpecDataX icolumn wrong');
  end
  datapath = [dir_statIn,cell2mat(file_names(k)) ];
  SpecData = load(datapath);
  icolumnSet=icolumn(icolumnt,2);
  icolumnEnd=icolumn(icolumnt,3);  
  SpectraX(:,icolumnSet:icolumnEnd)=SpectraX(:,icolumnSet:icolumnEnd)+SpecData;
  disp( ['read:   ',datapath,'  sucessfully'] );
end
SpectraX=SpectraX/(ntimeX*utaufinal^2);
SpectraX(:,4:4:end)=SpectraX(:,4:4:end)/utaufinal^2;
myformat=[repmat('%24.15E',1,sumColumn/4+1),'\n'];

% ux in x-dir
fid=fopen([dir_statOut,'specx_u.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nxhp
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpectraX(k,1:4:end-3)]);
end
fclose(fid);

% uy in x-dir
fid=fopen([dir_statOut,'specx_v.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_node]);
for k=1:nxhp
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpectraX(k,2:4:end-2)]);
end
fclose(fid);

% uz in x-dir
fid=fopen([dir_statOut,'specx_w.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nxhp
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpectraX(k,3:4:end-1)]);
end
fclose(fid);

% pr in x-dir
fid=fopen([dir_statOut,'specx_p.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nxhp
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpectraX(k,4:4:end)]);
end
fclose(fid);
end

%% Energy Spectra in z-dir
nzh=nzc/2; nzhp=nzh+1;
SpectraZ=zeros(nzhp,njSpec*4);
dir_output=dir(fullfile(dir_statIn,'SpecZ*'));
file_names={dir_output.name};
file_num=length(file_names);

if(file_num>0)
file_ave=0;
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))>iTimeSet)
   file_ave=file_ave+1;   
  end   
end
icolumn=zeros(file_ave,3);
columnCount=0;
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))<=iTimeSet);continue;end;
  icolumnt = str2double(datapath(17:20));
  CountFlag=1;
  if(columnCount>0)
    for m=1:columnCount
      if(icolumnt==icolumn(m,1))
        CountFlag=0;
        break
      end
    end
  end
  if(CountFlag==1)
    columnCount=columnCount+1;
    datapath = [dir_statIn,cell2mat(file_names(k)) ];
    SpecData = load(datapath);
    if(size(SpecData,1)~=nzhp)
      error('SpecDataZ size wrong');
    end
    icolumn(columnCount,1)=icolumnt;
    icolumn(columnCount,3)=size(SpecData,2);
  end
end
if(mod(file_ave,columnCount)~=0) 
  error('SpecDataZ columnCount wrong');
end
ntimeZ=file_ave/columnCount;
icolumn=icolumn(1:columnCount,1:3);
icolumn=sortrows(icolumn);
sumColumn=0;
for k=1:columnCount
  sumColumn   = sumColumn+icolumn(k,3);
  icolumn(k,3)= sumColumn;
end
icolumn(1,2)=1;
if(columnCount>1)
  for k=2:columnCount
    icolumn(k,2)=icolumn(k-1,3)+1;
  end
end
for k=1:file_num
  datapath = cell2mat(file_names(k));
  if(str2double(datapath(6:15))<=iTimeSet);continue;end;
  icolumnt = str2double(datapath(17:20))+1;
  if(icolumn(icolumnt,1)~=icolumnt-1)
    error('SpecDataZ icolumn wrong');
  end
  datapath = [dir_statIn,cell2mat(file_names(k)) ];
  SpecData = load(datapath);
  icolumnSet=icolumn(icolumnt,2);
  icolumnEnd=icolumn(icolumnt,3);  
  SpectraZ(:,icolumnSet:icolumnEnd)=SpectraZ(:,icolumnSet:icolumnEnd)+SpecData;
  disp( ['read:   ',datapath,'  sucessfully'] );
end
SpectraZ=SpectraZ/(ntimeZ*utaufinal^2);
SpectraZ(:,4:4:end)=SpectraZ(:,4:4:end)/utaufinal^2;
myformat=[repmat('%24.15E',1,sumColumn/4+1),'\n'];

% ux in z-dir
fid=fopen([dir_statOut,'specz_u.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nzhp
  fprintf(fid,myformat,[2*pi/zlz*(k-1),SpectraZ(k,1:4:end-3)]);
end
fclose(fid);

% uy in z-dir
fid=fopen([dir_statOut,'specz_v.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_node]);
for k=1:nzhp
  fprintf(fid,myformat,[2*pi/zlz*(k-1),SpectraZ(k,2:4:end-2)]);
end
fclose(fid);

% uz in z-dir
fid=fopen([dir_statOut,'specz_w.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nzhp
  fprintf(fid,myformat,[2*pi/zlz*(k-1),SpectraZ(k,3:4:end-1)]);
end
fclose(fid);

% pr in z-dir
fid=fopen([dir_statOut,'specz_p.txt'],'wt');
fprintf(fid,myformat,[-99,yplus_spec_center]);
for k=1:nzhp
  fprintf(fid,myformat,[2*pi/zlz*(k-1),SpectraZ(k,4:4:end)]);
end
fclose(fid);
end
