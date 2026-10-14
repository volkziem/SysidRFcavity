% pulsed_mode.m, V. Ziemann, 260410
clear all; pause(0.001)
format short
degree=180/pi;
tic
fun=@(t,p)p(3)+p(1).*(1-exp(-p(2)*t));
fun2=@(t,p)p(1).*(1-exp(-p(2)*t));
fun3=@(t,p)p(3)+p(1).*exp(-p(2)*t);

showiq=0;
sysid=0;
pulsefit=1;
do_sim=0;
Nforget=5000000;  
alpha=1-1/Nforget;
x=0;             % previous voltage
qhat=[0;0];      % fit parameters
P=eye(2);
p2=1; qhat2=0; alpha2=alpha; %1-1/100;
PP=eye(2); qqhat=[0;0];
Gamma=[-0.07503; -0.02604]; % determined with sysid=-2 from 150143 
Gamhat=[Gamma(1), -Gamma(2); Gamma(2), Gamma(1)];
Ahat=Gamhat/2;
P3=eye(2); qhat3=[0;0];
P4=eye(2); qhat4=[0;0];

Cavity='VTA data';

%fname = "../data/SDRdata/L2-140Pulse_260224_121629.txt"; % 8/2, 1.6 MV, Qe=2.0e10 
%fname = "../data/SDRdata/L2-140Pulse_260224_122500.txt"; % 2/2, 1.6 MV, Qe=7.8e9
%fname = "../data/SDRdata/L2-140Pulse_260224_145305.txt"; % 2/2, 6 MV, Qe=2.07e10
fname = "../data/SDRdata/L2-140Pulse_260224_150143.txt"; % 8/2, 6 MV, Qe=2.07e10
%fname = "../data/SDRdata/L2-140Pulse_260224_150855.txt"; % 8/2, 6 MV, Qe=2.07e10

yylim=[-12,12];
if contains(fname,'121629')
  phaseoff21=(138)/degree;  % 121629
  phaseoff31=(273)/degree;
  phaseoff41=(212)/degree;
  Rbeta=2.14;
elseif contains(fname,'122500')
  phaseoff21=(140)/degree;  % 121250
  phaseoff31=(300)/degree;
  phaseoff41=(172)/degree;
  Rbeta=2.14;
elseif contains(fname,'145305')
  phaseoff21=(300)/degree;  % 145305
  phaseoff31=(85)/degree;
  phaseoff41=(14)/degree;
  Rbeta=2.25;
elseif contains(fname,'150143') || contains(fname,'150855') 
  phaseoff21=(300)/degree;  % 150143 and 150855
  phaseoff31=(80)/degree;
  phaseoff41=(14)/degree;
  Rbeta=2.25;
end


x=0;             % previous voltage
qhat=[0;0];      % fit parameters
P=eye(2); 

a=1;
b=1;

filter=0;

skip=1000; %115000+4*200*1000;
kmax=skip+1800*1000; % maximum 15'800'000
%kmax=skip+100*1000-1000; 
stride=1;  % use every stride data point
if pulsefit==1 && stride>1, disp('***Set stride to unity in pulsefit mode!'), return; end 

fid = fopen(fname,"r");
if fid == -1, error("Cannot open file: %s", fname); end
line1=fgetl(fid); line1=line1(45:end); disp(line1)
line2=fgetl(fid);
parts = strsplit(strtrim(line2), '\t');
dt = str2double(parts{2});
LO = str2double(parts{4});
disp(['dt = ',num2str(dt),' s   LO = ',num2str(LO),' MHz'])
omegahat=2*pi*LO*1e6;
if stride>1
  dt=stride*dt;
  disp(['dt changed to ',num2str(dt),' s because  stride = ',num2str(stride)])
end
line3=fgetl(fid);

for k=1:skip-1, line=fgetl(fid); end
c=fscanf(fid,'%f',8); 
ampl1_prev=hypot(c(2),c(1));

data=zeros(kmax-skip,25);
iqdata=zeros(kmax-skip,15);
pdata=zeros(100,8); pdata2=zeros(100,6); pulse=0;
if do_sim==1
  simdata=zeros(kmax-skip,4); 
  Vsim=0; Vsim2=0;
  beta=1.09;
  om12dt=0.35*dt;
  QL=omegahat*dt/(2*om12dt);
  QE=QL*(1+1/beta); omegaE=omegahat/QE;
  Q0=beta*QE;
  disp(['Sim: Q_0 = ',num2str(Q0,3),'   Q_E = ',num2str(QE,3), ...
    '   Q_L = ',num2str(QL,3),'   beta = ',num2str(beta,4)])
end

k=0; trigger=-1; schluss=false; 

for kk=1:stride:kmax-skip
  for i1=1:stride
    cc=fscanf(fid,'%f',8);
    if isempty(cc), schluss=true; break; end
    c=(filter*c+cc)/(filter+1);
    c=1000*c;
  end
  if schluss, disp(kk); break; end
  
  ampl1=hypot(c(2),c(1));
  phase1=atan2(c(2),c(1));
  %q=ampl1-ampl1_prev;   [q,q>2e-4]

  if ampl1-ampl1_prev>1, trigger=kk; end
  if ampl1_prev-ampl1>1, trigger=0; end

  k=max(1,k+1);
  
  data(k,1:8)=c';
  ampl2=hypot(c(4),c(3));    Ip=a*ampl2;  % include calibration
  phase2=atan2(c(4),c(3));   
  ampl3=hypot(c(6),c(5));    Im=b*ampl3;  % ditto
  phase3=atan2(c(6),c(5));
  ampl4=hypot(c(8),c(7));    Vc=ampl4;    % cavity voltage
  phase4=atan2(c(8),c(7));
  data(k,9:12)=[ampl1;ampl2;ampl3;ampl4];

  if 1   % IQ data
    phaseref=phase1;
    phase1=mod(phase1-phaseref,2*pi); 
    phase2=mod(phase2-phaseref,2*pi)-phaseoff21;
    phase3=mod(phase3-phaseref,2*pi)-phaseoff31;
    phase4=mod(phase4-phaseref,2*pi)-phaseoff41;
    iqdata(k,1)=ampl1*cos(phase1);
    iqdata(k,2)=ampl1*sin(phase1);
    iqdata(k,3)=ampl2*cos(phase2);
    iqdata(k,4)=ampl2*sin(phase2);  IIp=a*[iqdata(k,3);iqdata(k,4)];
    iqdata(k,5)=ampl3*cos(phase3);
    iqdata(k,6)=ampl3*sin(phase3);  IIm=b*[iqdata(k,5);iqdata(k,6)];
    iqdata(k,7)=ampl4*cos(phase4);
    iqdata(k,8)=ampl4*sin(phase4);  VVc=[iqdata(k,7);iqdata(k,8)];
    iqdata(k,9)=phase1*degree;
    iqdata(k,10)=phase2*degree;
    iqdata(k,11)=phase3*degree;
    iqdata(k,12)=phase4*degree;
  end

  if do_sim==1
    Isim=3.05*ampl1; % make the forward currents agree
    Ifwd=Isim*(1+beta)/(2*beta);
    Vsim=(1-om12dt)*Vsim+Rbeta*Ifwd*omegaE*dt;
    Iref=abs(Vsim-0.5*(1+beta)*Rbeta*Isim/beta)/Rbeta;
    simdata(k,1)=Isim;
    simdata(k,2)=Ifwd;
    simdata(k,3)=Iref;
    simdata(k,4)=Vsim;
    simdata(k,5)=Vsim2;
    %data(k,9:12)=simdata(k,1:4);
  end
  
  if sysid==1 && ampl4>1e-3 % SEL-sysid for omegaE and omega0
    xnew=Vc;
    %y2=xnew-x; 
    if k>1, y=data(k,12)-data(k-1,12); else y=data(k,12); end 
    % G=[-0.5*x+Rbeta*Ip,-0.5*x];
    if k>1, IIIp=data(k-1,10); else IIIp=0; end
    G=[-0.5*x+Rbeta*IIIp,-0.5*x];
    tmp2=eye(2)-P*G'*inv(alpha*eye(1)+G*P*G')*G;    
    Pnew=tmp2*P/alpha;                              
    qhat=tmp2*(qhat+P*G'*y/alpha);
    x=xnew;
    P=Pnew;
    data(k,14:15)=qhat';
    data(k,16)=trace(P'*P)/4;
  else
    x=Vc;
    data(k,14:16)=NaN;
  end

  if sysid==2   % circulator reflections
   % if ampl1<2e-4 && ampl4>2e-2
    if ampl1<2e-4 && ampl4>1e-3   % 3MV, 100% AM
      % Gamma0=2*Ip/Im;
      % if abs(Gamma)<1e-10, Gamma=Gamma0; end
      % Gamma=(maver*Gamma+Gamma0)/(maver+1);
      yy=2*Ip; GG=Im;
      tmp2=1-p2*GG^2/(alpha2+p2*GG^2);
      p2new=tmp2*p2/alpha2;
      qhat2=tmp2*(qhat2+p2*GG*yy/alpha2);
      p2=p2new;
      data(k,19)=1;
    end
    data(k,17)=qhat2;  % circulator reflection Gamma
    data(k,18)=p2;
  end

  if sysid==-2   % circulator reflection with IQ
    if ampl1<2e-4 && ampl4>1e-3   % 3MV, 100% AM
      yy=2*IIp;
      GG=[IIm(1),-IIm(2);IIm(2),IIm(1)];
      tmp2=eye(2)-PP*GG'*inv(alpha*eye(2)+GG*PP*GG')*GG;
      PPnew=tmp2*PP/alpha;
      qqhat=tmp2*(qqhat+PP*GG'*yy/alpha);
      PP=PPnew;
    end
    iqdata(k,13:14)=qqhat';
    iqdata(k,15)=0.25*trace(PP'*PP);
  end

  if sysid==3 % directional coupler crosstalk   
    if ampl1>2e-4 && k-trigger > 200  % fit for a+c and b+d  on flattop
      yy=Vc/Rbeta-Ip-Im;
      G=[Ip,Im];
      tmp3=eye(2)-P3*G'*inv(alpha+G*P3*G')*G;
      P3new=tmp3*P3/alpha;                              
      qhat3=tmp3*(qhat3+P3*G'*yy/alpha);
      P3=P3new;
    end
    data(k,20:21)=qhat3';
    data(k,22)=0.25*trace(P3'*P3);
    if ampl1<2e-4 && Vc>0.001  % fit for a and b on falling edge
      y4=qhat3(1)*Ahat*Ip-Ip+Ahat*Im*(1+qhat3(2));
      G=[Ip+Ahat*Ip,Im+Ahat*Im]
      tmp4=eye(2)-P4*G'*inv(alpha+G*P4*G')*G;   % bombs G
      P4new=tmp4*P4/alpha;
      qhat4=tmp4*(qhat4+P4*G'*y4/alpha);
      P4=P4new; 
    end
    data(k,23:24)=qhat4';
    data(k,25)=0.25*trace(P4'*P4);
  end

  if sysid==-3   % directional coupler crosstalk with IQ data
    if ampl1>2e-4 && k-trigger > 200   % fit for a+c and b+d  on flattop
      yy3=VVc/(Rbeta)-IIp-IIm;
      G3=[IIp,IIm];
      tmp3=eye(2)-P3*G3'*inv(alpha*eye(2)+G3*P3*G3')*G3;   
      qhat3=tmp3*(qhat3+P3*G3'*yy3/alpha);
      P3=tmp3*P3/alpha;                                   
    end
    data(k,20:21)=qhat3';
    data(k,22)=0.25*trace(P3'*P3);
    if ampl1<2e-4 && Vc>0.001  % fit for a and b on falling edge
      yy4=qhat3(1)*Ahat*IIp-IIp+Ahat*IIm*(1+qhat3(2));
      G4=[IIp+Ahat*IIp,IIm+Ahat*IIm];
      tmp4=eye(2)-P4*G4'*inv(alpha*eye(2)+G4*P4*G4')*G4;
      P4new=tmp4*P4/alpha;
      qhat4=tmp4*(qhat4+P4*G4'*yy4/alpha);
      P4=P4new;   
    end
    data(k,23:24)=qhat4';
    data(k,25)=0.25*trace(P4'*P4); 
  end

 if pulsefit==1 && trigger>0    % do the omega12 fit
    data(k,13)=k-trigger;
    if k-trigger==380000/4   % acquired data for a new pulse  
      pulse=pulse+1;
      mm1=trigger:k; ttt=(mm1-trigger)*dt;
      chisq_fun=@(p)sum((data(mm1,12)-fun(ttt,p)').^2);
      [pfit,fval]=fminsearch(chisq_fun,[8,0.37,8]);
      f12=pfit(2)/(2*pi);
      am1=pfit(1)+pfit(3);  % Cavity asymptotic
      disp(['Pulse ',num2str(pulse),':  Data from k=',num2str(trigger),' to ',num2str(k), ... 
        '   f_{12} = ',num2str(f12),' Hz, max_{asymp} = ',num2str(pfit(1)+pfit(3))])
      disp(['   Vcav   ',num2str(pfit),' ',num2str(am1)])
      pdata(pulse,1)=f12;
      pdata(pulse,2)=pfit(1);
      pdata(pulse,3)=pfit(3);
      QL=omegahat/(2*pfit(2));

      chisq_fun3=@(p)sum((data(mm1,11)-fun3(ttt,p)').^2);
      [pfit3,fval3]=fminsearch(chisq_fun3,[4,0.37,0.4]);
      am3=pfit3(3);           % Reflected asymptotic
      disp(['   Vref   ',num2str(pfit3),' ',num2str(am3)])
      am2=mean(data(mm1,10)); % Forward asymptotic
      GAM=am3/am2;            % reflection coefficient
      BET=(1+GAM)/(1-GAM);    % coupling
      disp(['   Vfwd   ',num2str(GAM),' ',num2str(BET),' ',num2str(am2)])     
      RRb=am1/(am2+am3);
      disp(['   Rbeta = ',num2str(RRb)])
      pdata(pulse,4)=RRb;
      pdata(pulse,5)=am2+am3; 
      pdata(pulse,6)=GAM;
      pdata(pulse,7)=BET;
      pdata(pulse,8)=QL;
      pdata2(pulse,1:3)=pfit;
      pdata2(pulse,4:6)=pfit3;
      
      %
      if 0   % fit with two parameters, pedestal subtracted
        mmoff=trigger-10:trigger-2;
        offset=mean(data(mmoff,12));
        chisq_fun2=@(p)sum((data(mm1,12)-offset-fun2(ttt,p)').^2);
        [pfit,fval]=fminsearch(chisq_fun2,[0.05,220]);
        f12=pfit(2)/(2*pi);
        disp(['Pulse ',num2str(pulse),':  Data from k=',num2str(trigger),' to ',num2str(k), ... 
          '   f_{12} = ',num2str(f12),' Hz, max_{asymp} = ',num2str(pfit(1)+pfit(3))])  
        pdata(pulse,1)=f12;
        pdata(pulse,2)=pfit(1);
        pdata(pulse,3)=offset; % pfit(3);
        pdata(pulse,4)=min(data(mm1,12))/max(data(mm1,12));
      end    
     
      trigger=0;
      %if abs(f12-0.05)>0.02, pulse=pulse-1; disp('*** Pulse dropped ***'), end
      %if abs(f12-0.003)>0.0015, pulse=pulse-1; disp('*** Pulse dropped ***'), end
    end
  end
  ampl1_prev=ampl1;
end
fclose(fid);
kmax=k;
mm=1:kmax;
mm1=1:kmax/2; mm2=floor(kmax/2)+1:kmax; tt1=mm1*dt; tt2=mm2*dt;
tt=mm*dt;
toc
if pulse>0
  tmp=strsplit(fname,'L2-140Pulse_260224');
  outnam='pdata' + tmp(2);
  outnam2='pdata2' + tmp(2);
  dlmwrite(outnam,[(1:pulse)',pdata(1:pulse,:)],'\t');
  dlmwrite(outnam2,[(1:pulse)',pdata2(1:pulse,:)],'\t');
  om12=mean(pdata2(1:pulse,2));
  QL=mean(pdata(1:pulse,8));
  Rbeta_mean=mean(pdata(1:pulse,4));
  beta_mean=mean(pdata(1:pulse,7));
  disp([om12,QL/1e10,Rbeta_mean,beta_mean])
end
%dlmwrite('tmp.txt',[tt',data(:,10),data(:,12)],'\t')

if 0  % Raw signals I and Q
  figure(1); clf  
  subplot(4,1,1)
  plot(tt,data(mm,1),'k',tt,data(mm,2),'r')
  ylabel('Drive'); xlim([min(tt),max(tt)])
  title(line1);
  subplot(4,1,2)
  plot(tt,data(mm,3),'k',tt,data(mm,4),'r')
  ylabel('Fwd'); xlim([min(tt),max(tt)])
  title(fname,'Interpreter','none')
  subplot(4,1,3)
  plot(tt,data(mm,5),'k',tt,data(mm,6),'r')
  ylabel('Ref'); xlim([min(tt),max(tt)])
  subplot(4,1,4)
  plot(tt,data(mm,7),'k',tt,data(mm,8),'r')
  ylabel('Probe');  xlabel('time [s]'); xlim([min(tt),max(tt)])
  sgtitle('Raw I+Q signals')
end

if 1  % Raw amplitudes of Gen, FwD, Ref, Probe
  figure(1); clf
  subplot(3,1,1); 
  plot(tt,data(mm,9),'k')
  if do_sim==1, hold on, plot(tt,simdata(mm,1),'b'); hold off, end
  xlim([min(tt),max(tt)]);
  xlabel('Time [s]'); ylabel('Gen')
  title(line1);

  subplot(3,1,2); 
  plot(tt,data(mm,10),'k',tt,data(mm,11),'r') 
  if do_sim==1, hold on, plot(tt,simdata(mm,2),'b',tt,simdata(mm,3),'m'); hold off, end
  xlim([min(tt),max(tt)])
  %m=mean(data(:,10)); ylim([m-0.007,m+0.007])
  %ylim([0.03,0.036])
  xlabel('Time [s]'); ylabel('Raw amplitude')
  legend('Fwd','Ref')
  title(fname,'Interpreter','none')

  subplot(3,1,3); 
  plot(tt,data(mm,12),'k')
  if do_sim==1, hold on, plot(tt,simdata(mm,4),'b',tt,simdata(mm,5),'r'); hold off; end
  xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('Probe')
end

if 0  % pulse histogram
  figure(5); clf
  meanf=mean(pdata(1:pulse,1)); stdf=std(pdata(1:pulse,1));
  f12mean=meanf;
  QL=omegahat/(4*pi*f12mean);
  edges=meanf-3*stdf:6*stdf/30:meanf+3*stdf;
  subplot(6,1,1); histogram(pdata(1:pulse,1),edges); xlabel('f_{12} [Hz]')
  legend(['mean = ',num2str(meanf),'+/-',num2str(stdf),' Hz'])
  title(line1)
  subplot(6,1,2); plot(1:pulse,pdata(1:pulse,1),'k'); 
  xlim([0,pulse]); 
  title(fname,'Interpreter','none')
  xlabel('pulse'); ylabel('f_{12} [Hz]')
  meanf=mean(pdata(1:pulse,2)); stdf=std(pdata(1:pulse,2));
  edges=meanf-3*stdf:6*stdf/30:meanf+3*stdf;
  subplot(6,1,3); histogram(pdata(1:pulse,2),edges); xlabel('Amplitude')
  legend(['mean = ',num2str(meanf),'+/-',num2str(stdf)])
  subplot(6,1,4); plot(1:pulse,pdata(1:pulse,2),'k'); xlim([0,pulse])
  xlabel('pulse'); ylabel('Ampl')
  meanf=mean(pdata(1:pulse,3)); stdf=std(pdata(1:pulse,3));
  edges=meanf-3*stdf:6*stdf/30:meanf+3*stdf;
  subplot(6,1,5); histogram(pdata(1:pulse,3),edges); xlabel('Pedestal')
  legend(['mean = ',num2str(meanf),'+/-',num2str(stdf)])
  subplot(6,1,6); plot(1:pulse,pdata(1:pulse,3),'k'); xlim([0,pulse])
  xlabel('pulse'); ylabel('Pedestal')
  title(['AM ratio LO/HI = ',num2str(100*mean(pdata(:,4)),3),' %'])
end

if 1
  figure(2); clf
  subplot(3,1,1);
  plot(tt,data(mm,14),'k')
  xlim([min(tt),max(tt)]); ylim([0,1e-3])
  xlabel('Time [s]'); ylabel('\omega_E \Deltat')

  subplot(3,1,2);
  plot(tt,data(mm,15),'k')
  xlim([min(tt),max(tt)]); ylim([0,1e-3])
  xlabel('Time [s]'); ylabel('\omega_0 \Deltat')

  subplot(3,1,3);
  loglog(tt,data(mm,16),'k')
  xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('|P|')
  pause(0.3)
end

if sysid==1
  phases=mean(data(mm2,14:15));
  figure(3); clf
  subplot(2,1,1);
  QE=omegahat*dt./data(mm,14);
  QQE=mean(QE(end-100:end));
  semilogy(tt,QE,'k')
  xlim([min(tt),max(tt)]); ylim([1e8,1e12])
  xlabel('Time [s]'); ylabel('Q_E')
  title(line1);
  text(0.2*max(tt),5e9,['\omega_E\Deltat = ',num2str(phases(1),4), ...
    ' and  \omega_0\Deltat = ',num2str(phases(2),3)])

  subplot(2,1,2);
  Q0=omegahat*dt./data(mm,15);
  z=sum(Q0<0);
  z2=sum(Q0(mm2)<0);
  QQ0=mean(Q0(end-100:end));
  if z2>0 
    semilogy(tt,abs(Q0),'r')
    text(0.5*max(tt),1e11,'Q_0 is negative!')
  else
    semilogy(tt1,abs(Q0(mm1)),'r',tt2,abs(Q0(mm2)),'b') 
  end
  xlim([min(tt),max(tt)]); ylim([1e8,1e12])
  xlabel('Time [s]'); ylabel('Q_0')
  title(fname,'Interpreter','none')

  figure(4); clf
  if z2>0 
    semilogy(tt,QE,'k',tt,abs(Q0),'r')
    text(0.5*max(tt),5e8,'Q_0 is negative!')
  else
    %semilogy(tt,QE,'k',tt,Q0,'b')
    semilogy(tt,QE,'k',tt2,abs(Q0(mm2)),'b',tt1,abs(Q0(mm1)),'r') 
  end
  xlim([min(tt),max(tt)]); ylim([1e8,1e12])
  xlabel('Time [s]'); ylabel('Q_E, Q_0')
  legend(['Q_E = ',num2str(QQE,'%9.3e')],['Q_0 = ',num2str(QQ0,'%9.3e')]);
  QQL=1/(1/QQE+1/QQ0);
  text(0.5*max(tt),1e9,['QQ_L = ',num2str(QQL,'%8.2e')])
  if exist('QL') 
    text(0.5*max(tt),5e8,['Q_L = ',num2str(QL,'%8.2e'),' or QQ_L = ',num2str(QQL,'%8.2e')]); 
  end
  title(line1);
end

if sysid==2   % circulator reflections
  figure(2); clf
  subplot(3,1,1);
  plot(tt,data(mm,17),'k'); xlim([min(tt),max(tt)])
  ylim([0,1.05*max(data(mm,17))])
  xlabel('Time [s]'); ylabel('\Gamma')
  Gamma=mean(data(k/2:k,17));
  title(['\Gamma = ',num2str(Gamma)])
  subplot(3,1,2);
  plot(tt,data(mm,19),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]');
  subplot(3,1,3);
  loglog(tt,data(mm,18),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('p_2')
end

if sysid==-2   % circulator with I and Q
  figure(3); clf
  subplot(2,1,1)
  plot(tt,iqdata(mm,13),'k',tt,iqdata(mm,14),'r');
  xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('\Gamma')
  Gammar=mean(iqdata(floor(k/2):k,13));
  Gammai=mean(iqdata(floor(k/2):k,14));
  title(['\Gamma = ',num2str(Gammar),' + ',num2str(Gammai),' i'])
  subplot(2,1,2);
  loglog(tt,iqdata(mm,15),'k'); xlim([min(tt),max(tt)])
  xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('|PP|')
end

if abs(sysid)==3 % directional coupler
  figure(2)
  subplot(3,1,1)
  plot(tt,data(mm,20),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('a+c'); ylim([-0.009,0.009])
  s1=mean(data(k/2:k,20)); title(['s_1 = a+c = ',num2str(s1)])
  subplot(3,1,2)
  plot(tt,data(mm,21),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]');  ylabel('b+d'); ylim([-0.009,0.009])
  s2=mean(data(k/2:k,21)); title(['s_2 = b+d = ',num2str(s1)])
  subplot(3,1,3)
  loglog(tt,sqrt(data(mm,22)),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]');ylabel('|P_3|')
  figure(3)
  subplot(3,1,1)
  plot(tt,data(mm,23),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('a'); ylim([-0.12,0.012])
  a=mean(data(k/2:k,23)); c=s1-a;
  title(['a = ',num2str(a),'     c = ',num2str(c)])
  subplot(3,1,2)
  plot(tt,data(mm,24),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('b'); ylim([-0.012,0.012])
  b=mean(data(k/2:k,24)); d=s2-b;
  title(['b = ',num2str(b),'     d = ',num2str(d)])
  subplot(3,1,3)
  loglog(tt,sqrt(data(mm,25)),'k'); xlim([min(tt),max(tt)])
  xlabel('Time [s]'); ylabel('|P_4|')
end

if showiq==1  % Demodulated signals
  figure(2); clf
  subplot(4,1,1)
  plot(tt,iqdata(mm,1),'k',tt,iqdata(mm,2),'r')
  ylabel('Drive I,Q'); xlim([min(tt),max(tt)])
  title(line1);
  subplot(4,1,2)
  plot(tt,iqdata(mm,3),'k',tt,iqdata(mm,4),'r')
  ylabel('Fwd I,Q'); xlim([min(tt),max(tt)]); % ylim([0.0255,0.0268])
  title(fname,'Interpreter','none')
  subplot(4,1,3)
  plot(tt,iqdata(mm,5),'k',tt,iqdata(mm,6),'r')
  ylabel('Ref I,Q'); xlim([min(tt),max(tt)])
  subplot(4,1,4)
  plot(tt,iqdata(mm,7),'k',tt,iqdata(mm,8),'r')
  ylabel('Probe I,Q'); xlabel('time [s]'); xlim([min(tt),max(tt)])
  %ylim([0.0995,0.103])
  sgtitle('I and Q signals')
end
if showiq==1   % Phases versus time
  figure(4); clf
  subplot(4,1,1); plot(tt,iqdata(mm,9),'k'); ylabel('Phase(Drive)')
  xlim([min(tt),max(tt)])
  title(line1);
  subplot(4,1,2); plot(tt,iqdata(mm,10),'k'); ylabel('Phase(Fwd)')
  title(fname,'Interpreter','none')
  xlim([min(tt),max(tt)]); ylim(yylim)
  subplot(4,1,3); plot(tt,iqdata(mm,11),'k'); ylabel('Phase(Ref)')
  xlim([min(tt),max(tt)]); ylim(yylim)
  subplot(4,1,4); plot(tt,iqdata(mm,12),'k'); ylabel('Phase(Probe)')
  xlim([min(tt),max(tt)]); ylim(yylim)
  xlabel('time [s]'); sgtitle('Phases');
end

if 0 && pulse>0
  p=[(1:pulse)',pdata(1:pulse,:)];
  figure(4); clf

  subplot(5,1,1);
  plot(p(:,1),p(:,6),'k')
  ylabel('Itot')
  %title(['Mean = ',num2str(mean(p(:,6)))])
  title(fname,'Interpreter','none')

  subplot(5,1,2);
  plot(p(:,1),p(:,5),'k')
  ylabel('Rbeta')
  title(['Mean = ',num2str(mean(p(:,5)))])

  subplot(5,1,3);
  plot(p(:,1),p(:,7),'k')
  ylabel('\Gamma')
  title(['Mean = ',num2str(mean(p(:,7)))])

  subplot(5,1,4);
  plot(p(:,1),p(:,8),'k')
  ylabel('\beta')
  title(['<\beta> = ',num2str(mean(p(:,8)))])

  subplot(5,1,5);
  plot(p(:,1),p(:,9),'k')
  xlabel('Pulse'); ylabel('Q_L')
  title(['<Q_L> = ',num2str(mean(p(:,9)),3)])

  beta=mean(p(:,8))
  QL=mean(p(:,9))
  QE=QL*(1+1/beta)
  Q0=beta*QE
  title(['<Q_L> = ',num2str(mean(p(:,9)),3),'    Q_0 = ',num2str(Q0,3),'    Q_E = ',num2str(QE,3)])
end

toc