% calculate_mu_PI.m, V. Ziemann, 230811
  CC=eye(2)+Areal+Breal*Bplus-Breal*Ki;
  DD=Areal+Breal*Bplus;  
  Z1=0.5*CC+sqrtm(0.25*CC^2-DD);
  Z2=0.5*CC-sqrtm(0.25*CC^2-DD);
  SS=inv(Z1-Z2);
  R1=SS*(Z1-eye(2));
  R2=SS*(Z2-eye(2)); 
  A11=inv(eye(2)-Z1*Z1');
  A12=inv(eye(2)-Z1*Z2');
  A21=inv(eye(2)-Z2*Z1');
  A22=inv(eye(2)-Z2*Z2');

%   Evv=R1*A11*R1'-R1*A12*R2'-R2*A21*R1'+R2*A22*R2'
%   Evy=SS*(A11*(eye(2)-Z1')-A12*(eye(2)-Z2')-A21*(eye(2)-Z1')+A22*(eye(2)-Z2'))*SS'
%   Eyv=SS*((eye(2)-Z1)*A11-(eye(2)-Z1)*A12-(eye(2)-Z2)*A21+(eye(2)-Z2)*A22)*SS'
%   Eyy=SS*(A11-A12-A21+A22)*SS'

  
  Evv=R1*A11*R1'-R1*A12*R2'-R2*A21*R1'+R2*A22*R2';
  Evy=R1*A11*Z1'*SS'-R1*A12*Z2'*SS'-R2*A21*Z1'*SS'+R2*A22*Z2'*SS';
  Eyv=SS*Z1*A11*R1'-SS*Z1*A12*R2'-SS*Z2*A21*R1'+SS*Z2*A22*R2';
  Eyy=SS*(Z1*A11*Z1'-Z1*A12*Z2'-Z2*A21*Z1'+Z2*A22*Z2')*SS';

  h=zeros(2,2,2);
  h(1,:,:)=[-1+R*Kp,0;0,1]; 
  h(2,:,:)=[0,-1;-1+R*Kp,0];
  hb=zeros(2,2,2); 
  hb(:,:,1)=[-1+R*Kp,0;0,1];
  hb(:,:,2)=[0,-1+R*Kp;-1,0];

  jj=zeros(2,2,2);
  jj(1,:,:)=R*Ki*[1,0;0,0];
  jj(2,:,:)=R*Ki*[0,0;1,0];
  jjb=zeros(2,2,2);
  jjb(:,:,1)=R*Ki*[1,0;0,0];
  jjb(:,:,2)=R*Ki*[0,1;0,0];

  tmp=tensorprod(hb,Evv,3,1);
  hvvh=tensorprod(tmp,h,[3,2],[1,2]);
  tmp=tensorprod(hb,Evy,3,1);
  hvyj=tensorprod(tmp,jj,[3,2],[1,2]);
  tmp=tensorprod(jjb,Eyv,3,1);
  jyvh=tensorprod(tmp,h,[3,2],[1,2]);
  tmp=tensorprod(jjb,Eyy,3,1);
  jyyj=tensorprod(tmp,jj,[3,2],[1,2]);

  Q=hvvh-hvyj-jyvh+jyyj;
%   mu1=Q(1,1)
%   mu2=Q(2,2)
  ev=eig(Q);
  mu1=max(ev);
  mu2=min(ev);

  