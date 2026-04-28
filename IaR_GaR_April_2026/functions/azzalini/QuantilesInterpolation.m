function [lc,sc,sh,df,lc0,sc0,sh0,df0] = QuantilesInterpolation(qqTarg,QQ,lc0,sc0,sh0,df0);

% A*beta  <= B
B   = [];%-mR*eqw;
A   = [];%-mR;

% Aeq*beta  = Beq
Beq   = [];
Aeq   = [];

% LB <= beta <= UB 
LB = [-20  1   -20  2 ];
UB = [ 20  50   20  30 ];

options =  optimset('MaxFunEvals',1e+4,'MaxIter',1e+4,'TolX',1e-6,'TolCon',1e-6,'Display','off','Algorithm','interior-point' ,'LargeScale','on');


[temp,jq50] = min(abs(QQ-.50));
[temp,jq25] = min(abs(QQ-.25));
[temp,jq75] = min(abs(QQ-.75));
[temp,jq05] = min(abs(QQ-.05));
[temp,jq95] = min(abs(QQ-.95));



if nargin<3
    iqn = norminv(.75)-norminv(.25);    
    lc0 = qqTarg(jq50);
    sc0 = (qqTarg(jq75)-qqTarg(jq25))/iqn;
    sh0 = 0;
    df0 = UB(end);
end;


X0 = [lc0,sc0,sh0,df0];


Select = [jq05 jq25 jq50 jq75 jq95];

X = fmincon(@DistST,X0,A,B,Aeq,Beq,LB,UB,[],options,QQ(Select),qqTarg(Select));


lc = X(1);
sc = X(2);
sh = X(3);
df = round(X(4));

function dist = DistST(X,QQ,qqTarg);
location = X(1);
scale = X(2);
shape = X(3);
df = round(X(4));

qq = qskt(QQ,location,scale,shape,df);


dist = norm(qq-qqTarg);


%disp([qq;qqTarg])