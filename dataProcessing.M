clear
filenameSP='Historical S&P 500 Price.xlsx';
[price,pr_txt, ~]=xlsread(filenameSP);
time_tmp=[];
time_tmp=string(pr_txt(:,:));
time_pr=[];
for n=2:1812
time_pr(n-1)=datenum(time_tmp(n,1));
end
filenameIR='Historical Interest Rates.xls';
[rate, ir_txt, ~]=xlsread(filenameIR);
time_tmp=[];
time_tmp=string(ir_txt(:,:));
time_ir=[];
for n=12:871
time_ir(n-11)=datenum(time_tmp(n,1));
end

time_pr_min=min(time_pr); % 683370
time_pr_max=max(time_pr); % 738451
time_ir_min=min(time_ir); % 712224
time_ir_max=max(time_ir); % 738369

ID=find(time_pr>=time_ir_min & time_pr<=time_ir_max);
time_pr=time_pr(ID);price=price(ID);

23

time_pr=flipud(time_pr(:));price=flipud(price(:));
price=interp1(time_pr,price,time_ir);
time_pr=time_ir;

figure(1)
plot(time_pr,price,'r-','linewidth',1)
ylabel('S&P price','fontsize',10,'fontweight','bold')
xlabel('time (mm/dd/yy)','fontsize',10,'fontweight','bold')
grid on

xtic=[min(time_pr):floor((max(time_pr)-min(time_pr))/5):max(time_pr)];
set(gca,'xtick',xtic);

xtic_txt=[];
for n=1:length(xtic)
tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];
end

set(gca,'xticklabel',xtic_txt)

xlim([time_ir_min time_ir_max]);
figure(2)
plot(time_ir,rate,'b-','linewidth',1)
ylabel('interest rate','fontsize',10,'fontweight','bold')
xlabel('time (mm/dd/yy)','fontsize',10,'fontweight','bold')
grid on

24

xtic=[min(time_ir):floor((max(time_ir)-min(time_ir))/5):max(time_ir)];
set(gca,'xtick',xtic);

xtic_txt=[];
for n=1:length(xtic)
tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];
end

set(gca,'xticklabel',xtic_txt)
xlim([time_ir_min time_ir_max]);

figure(3)
plot(time_pr,log(price),'r-','linewidth',1)
ylabel('log(S&P price)','fontsize',10,'fontweight','bold')
xlabel('time (mm/dd/yy)','fontsize',10,'fontweight','bold')
grid on
xtic=[min(time_pr):floor((max(time_pr)-min(time_pr))/5):max(time_pr)];
set(gca,'xtick',xtic);

xtic_txt=[];
for n=1:length(xtic)
tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];
end

set(gca,'xticklabel',xtic_txt)

25

xlim([time_ir_min time_ir_max]);

figure(4)
plot(time_pr,detrend(log(price)),'r-','linewidth',1)
ylabel('log(S&P price)','fontsize',10,'fontweight','bold')
xlabel('time (mm/dd/yy)','fontsize',10,'fontweight','bold')
grid on

xtic=[min(time_pr):floor((max(time_pr)-min(time_pr))/5):max(time_pr)];
set(gca,'xtick',xtic);

xtic_txt=[];
for n=1:length(xtic)
tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];

end

set(gca,'xticklabel',xtic_txt)
xlim([time_ir_min time_ir_max]);
title('Logged (S&P) data with the trend being removed','fontsize',12,'fontweight','bold')

figure(5)
[AX,h1,h2]= plotyy(time_pr,detrend(log(price)),time_ir,rate)
xtic=[min(time_pr):floor((max(time_pr)-min(time_pr))/5):max(time_pr)];
set(AX(1),'xtick', xtic);
xtic_txt=[];
for n=1:length(xtic)

26

tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];
end
set(AX(1),'xticklabel',xtic_txt);

xtic=[min(time_ir):floor((max(time_ir)-min(time_ir))/5):max(time_ir)];
set(AX(2),'xtick', xtic);
xtic_txt=[];
for n=1:length(xtic)
tmp=datetime(xtic(n),'ConvertFrom','datenum');
formatOut = 'mm/dd/yy';
xtic_txt=[xtic_txt;datestr(tmp,formatOut)];
end
set(AX(2),'xticklabel',xtic_txt);
grid on
set(h1,'color','r','linestyle','-','linewidth',1)
set(h2,'color','b','linestyle','-','linewidth',1)
set(AX(1),'ycolor','r','fontsize',10);
set(AX(2),'ycolor','b','fontsize',10);
set(get(AX(1),'ylabel'),'string','log(S&P) anomalies');
set(get(AX(2),'ylabel'),'string','interest rate anomalies');
set(get(AX(1),'xlabel'),'string','time (mm/dd/yy)');
set(AX(1),'xlim',[time_ir_min time_ir_max]);
set(AX(2),'xlim',[time_ir_min time_ir_max]);

price=detrend(log(price));
[value, ID]=arrange(rate);
num=length(ID);

27

[pcon, prate ,pprice]=condp(rate(ID),price(ID),-1,-1);
write_out(time_ir(ID),rate(ID),price(ID));

close all
clear
load raw_data.mat
num=length(price_raw)
for n=2:num
P(n-1)=price_raw(n)/price_raw(n-1);
end
lamda=[-5:0.001:5];
Vtrans=boxcox(P,lamda,mean(P(:)));
ERROR=[];
for n=1:length(lamda)
Pa=[];Pa=P(:)-mean(P(:));
tmp=Vtrans(:,n);tmp=tmp(:);
ERROR(n)=sum((Pa-tmp).^2);
end
figure
plot(lamda, ERROR,'b-','linewidth',1)
xlabel('\lambda');
ylabel('{\Sigma} {(p_i-\phi_i)^2}')
grid minor
[Y, ID]=min(ERROR(:));
set(gca,'xtick',[-5 -4 -3 -2 -1 0 1 2 3 4 5]);
Lamda_opm=lamda(ID);

28

xlim([-5 5]);
hold on
ylim([0.03 0.13])
plot(Lamda_opm,ERROR(ID),'r*','linewidth',3);
plot([Lamda_opm Lamda_opm],[0.03 ERROR(ID)],'r-');
plot(Lamda_opm,0.03,'ro','linewidth',3);
text(-1.5,0.041,'minimum value');
text(-0.5,0.033,'(-0.577)','color','r','fontsize',4);
x=[-5:0.0001:5];
for n=1:length(x(:))
y(n)=erf(x(n));
end
y0=0.037594;
for n=1:length(y(:))-1
if (y(n+1)>=y0 & y(n)<=y0)
N=n;
break;
end
end
figure
Y=[y(N) y(N+1)];X=[x(N) x(N+1)];
x0=interp1(Y,X,y0);
plot(x,y,'b-');
xlabel('x');
ylabel('erf(x)');
hold on
set(gca,'xtick',[-5 -4 -3 -2 -1 0 1 2 3 4 5]);
grid on
figure;
r=normrnd(-0.00856,0.04085660912,29);
figure

29

histogram(r)
lam_op=-0.577;
figure;
y=(1+lam_op*r).^(1/lam_op);
histogram(y)
y=y(:);
P0=100;
for n=1:length(y)
if n==1
Pr(n)=P0*y(n);
else
Pr(n)=Pr(n-1)*y(n);
end
end
PPR=zeros(500,860);
for RM=1:500
for n=1:860
R(n)=randi([1 860],1);
end
ID=find(10<R<=243);
num=length(ID);
rate=zeros(1, 860);
for n=1:num
rate(1,ID(n))=normrnd(1.0078,0.174411,1);
while rate(1,ID(n)) == 1
rate(1,ID(n))=normrnd(1.0078,0.174411,1);
end
end
ID2=find(R<=10);
rate(1,ID2)=1;
Prr=zeros(1,860);
ID0=find(rate==0);

30

num1=length(ID0);
for n=1:num1
Prr(1,ID0(n))=normrnd(0.003340148321,0.04085660912,1);
while Prr(1,ID0(n)) >0.003340148321 + 0.04085660912*3
Prr(1,ID0(n))=normrnd(0.003340148321,0.04085660912,1);
end
end
IDg=find(rate>1);
num1=length(IDg);
for n=1:num1
Prr(1,IDg(n))=normrnd(0.00526,0.04085660912,1);
while Prr(1,IDg(n)) >0.00526 + 0.04085660912*3
Prr(1,IDg(n))=normrnd(0.00526,0.04085660912,1);
end
end
IDs=find(rate<1);
num1=length(IDs);
for n=1:num1
Prr(1,IDs(n))=normrnd(0.01524,0.04085660912,1);
while Prr(1,IDs(n)) > 0.01524 + 0.04085660912*3
Prr(1,IDs(n))=normrnd(0.01524,0.04085660912,1);
end
end
MIN=min(Prr(:));
IDe=find(rate==1);
num1=length(IDe);
for n=1:num1
Prr(1,IDe(n))=1.5*MIN;
end
Pr=zeros(1,860);Pr=Pr(:);
NN=length(Prr(:));
P0=10;
Prr=Prr(:);
for n=1:NN
y=(1+lam_op*Prr(n)).^(1/lam_op);
if n==1
Pr(n)=P0*y;
else

31

Pr(n)=Pr(n-1)*y;
end
end
PPR(RM,:)=Pr(:);
end
for RM=1:500
tmp=PPR(RM,:);tmp=tmp(:);
E(RM)=std(tmp-price_raw(:));
end
[Y I]=min(E(:));
[YY II]=max(E(:));
figure;
[AX h1 h2]=plotyy([1:860],price_raw,[1:860],PPR(II,:));
legend([h1 h2],'raw data','modelled data')
grid on
set(h1,'color','r','linestyle','-','linewidth',2)
set(h2,'color','b','linestyle','-','linewidth',2)
set(AX(1),'ycolor','r','fontsize',10);
set(AX(2),'ycolor','b','fontsize',10);
set(get(AX(1),'ylabel'),'string','raw data','fontsize',10);
set(get(AX(2),'ylabel'),'string','modelled one','fontsize',10);
