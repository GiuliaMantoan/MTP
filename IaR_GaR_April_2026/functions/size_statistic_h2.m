function[QVrejvecs,CVMrejvecs,QVrejvecs_bonf,CVMrejvecs_bonf]=size_statistic_h2(z,Kv,CVMv,H,stream1,rvec,el, bootMC)
% across the three different critical values - 1%, 5%, 10%. Not sure. 
QVrejvecs=zeros(1,3);
CVMrejvecs=QVrejvecs; 
QVrejvecs_bonf=QVrejvecs;
CVMrejvecs_bonf=QVrejvecs;
QVrejvecs_bonfh=zeros(1,H); CVMrejvecs_bonfh=QVrejvecs_bonfh;

% This requires pretty much the same inputs as Rossi-Sekhposyan (2019). 

%% Check what P should be here - index it in a way that the size of P is automatically taken into account. 
% Or is P always the same? The out-of-sample part should be fixed to equal
% the horizon?

indmaxKv=find(Kv==max(Kv), 1 );
indmaxCvM=find(CVMv==max(CVMv), 1 );
gstep=(1/H)/((H+1)/2);

% Some agnostic weighting scheme 
for h=1:H
    wi(h)=gstep*h; %#ok<AGROW>
end
w=flip(wi',1);

[tableboot_sup, tableboot_bonf, ~]  = CVfinalbootstrsuptest2(el, bootMC, z, rvec, w', wi,stream1);

%% Sup tests

QVrejvecs(1,1) = Kv(indmaxKv) > tableboot_sup(1,1);
CVMrejvecs(1,1) = CVMv(indmaxCvM) > tableboot_sup(1,2);

QVrejvecs(1,2) = max(Kv.*w') > tableboot_sup(2,1);
CVMrejvecs(1,2) = max(CVMv.*w') > tableboot_sup(2,2);

QVrejvecs(1,3) = max(Kv.*wi) > tableboot_sup(3,1);
CVMrejvecs(1,3) = max(CVMv.*wi) > tableboot_sup(3,2);

% %% Horizon by horizon
% QVrejvecs_h=zeros(1,H); CVMrejvecs_h=QVrejvecs_h;
% QVrejvecs_bonfh=QVrejvecs_h; CVMrejvecs_bonfh=QVrejvecs_h;
% for h=1:H
%     QVrejvecs_h(h) = Kv(h) > tableboot_h(1,h);
%     CVMrejvecs_h(h) = CVMv(h) > tableboot_h(2,h);
% end

%% Bonferroni Correction
for h=1:H
    QVrejvecs_bonfh(h)=Kv(h) > tableboot_bonf(1,h);
    CVMrejvecs_bonfh(h) = CVMv(h) > tableboot_bonf(2,h);
end

QVrejvecs_bonf(1,1)=Kv(indmaxKv) > tableboot_bonf(1,indmaxKv);
CVMrejvecs_bonf(1,1)=CVMv(indmaxCvM) > tableboot_bonf(2,indmaxCvM);

QVrejvecs_bonf(1,2) = max(Kv.*w') > max(tableboot_bonf(1,:).*w');
CVMrejvecs_bonf(1,2) = max(CVMv.*w') > max(tableboot_bonf(2,:).*w');

QVrejvecs_bonf(1,3) = max(Kv.*wi) > max(tableboot_bonf(1,:).*wi);
CVMrejvecs_bonf(1,3) = max(CVMv.*wi) > max(tableboot_bonf(2,:).*wi);

end
