function[QVrejvecs,CVMrejvecs,QVrejvecs_bonf,CVMrejvecs_bonf]=size_statistic_h2_shifted_h1rankfix(z,Kv,CVMv,H,stream1,rvec,el, bootMC,Kv_sim_h1,CVMv_sim_h1)
% 2026-08-13: real-data counterpart of size_statistic_h_shifted2_h1rankfix.m,
% built to match the ORIGINAL size_statistic_h2.m's exact interface and
% output shape (z, Kv, CVMv already computed from the observed/real PIT
% data -- no DGP simulation here, unlike the size_statistic_h* family used
% for the Monte Carlo size study) plus two extra inputs for the horizon-1
% fix (Kv_sim_h1, CVMv_sim_h1 -- precompute once via
% closedform_h1_cv_anchored.m using the SAME rvec passed here, it doesn't
% depend on the model/data).
%
% Swaps the original's CVfinalbootstrsuptest2 (no cross-horizon
% correlation modeling at all -- each horizon's bootstrap column built
% independently) for CVfinalbootstrsuptest_shifted_h1rankfix (index-
% aligned full-sharing cross-horizon correlation, no theta needed since
% there's no known DGP for real data, plus the anchored horizon-1 rank
% fix, which IS generic/safe for any correctly-specified model). See
% CVfinalbootstrsuptest_shifted_h1rankfix.m's header for the fuller
% rationale on why "shifted" (not "shifted2") is the appropriate choice
% here.
%
% Everything else (weighting scheme, Sup/Bonferroni test construction) is
% unchanged from the original size_statistic_h2.m.

QVrejvecs=zeros(1,3);
CVMrejvecs=QVrejvecs;
QVrejvecs_bonf=QVrejvecs;
CVMrejvecs_bonf=QVrejvecs;
QVrejvecs_bonfh=zeros(1,H); CVMrejvecs_bonfh=QVrejvecs_bonfh;

indmaxKv=find(Kv==max(Kv), 1 );
indmaxCvM=find(CVMv==max(CVMv), 1 );
gstep=(1/H)/((H+1)/2);
for h=1:H
    wi(h)=gstep*h; %#ok<AGROW>
end
w=flip(wi',1);

[tableboot_sup, tableboot_bonf, ~]  = CVfinalbootstrsuptest_shifted_h1rankfix(el, bootMC, z, rvec, w', wi, stream1, Kv_sim_h1, CVMv_sim_h1);

%% Sup tests
QVrejvecs(1,1) = Kv(indmaxKv) > tableboot_sup(1,1);
CVMrejvecs(1,1) = CVMv(indmaxCvM) > tableboot_sup(1,2);

QVrejvecs(1,2) = max(Kv.*w') > tableboot_sup(2,1);
CVMrejvecs(1,2) = max(CVMv.*w') > tableboot_sup(2,2);

QVrejvecs(1,3) = max(Kv.*wi) > tableboot_sup(3,1);
CVMrejvecs(1,3) = max(CVMv.*wi) > tableboot_sup(3,2);

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
