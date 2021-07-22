clc
clear
pi=3.14159265;
HzAmp=zeros;
HzAng=zeros;
Sx=zeros;
w=zeros;
H_z=zeros;
Yi=zeros;
fenzi1= [0,0];
fenmu1= [0,0];
N=input('请选择采样点数N:');
p=input('请选择AR模型阶数p(N/3<p<N/2):');
n=1:N;
t98=csvread('t98.csv');
bandpower=zeros;
for t=1:20
    for l=1:10
        wn=0.0001*randn(1,N);
        xn=t98((N*(10*(t-1)+l-1+00)+1):N*(10*(t-1)+l+00),:)'+wn;
        rou=zeros(p,p);
        e(1,:)=xn;
        b(1,:)=xn;
        sita(1,1)=mean(xn.^2);
        for m=1:p
            if(m==1)
                km=-2*(sum(e(m,m+1:N).*(b(m,(m:N-1)))))/sum(e(m,m+1:N).^2+b(m,m:N-1).^2);
            else
                km=-2*(sum(e(m,m:N).*(b(m,m-1:N-1))))/sum(e(m,m:N).^2+b(m,m-1:N-1).^2);
            end
            for i=1:m-1 
                rou(m,i)=rou(m-1,i)+km*rou(m-1,m-i);
            end
            rou(m,m)=km;
            for i=1:N
                if(i==1)
                    e(m+1,i)=e(m,i);
                    b(m+1,i)=km*e(m,i);
                else
                    e(m+1,i)=e(m,i)+km*b(m,i-1);
                    b(m+1,i)=b(m,i-1)+km*e(m,i);
                end
            end
            sita(m+1)=(1-km.^2)*sita(m);
        end
        G2=sita(p);
        ap=zeros(1,p);
        for k=1:p
            ap(k)=rou(p,k);
        end
        b1=[G2^0.5,0,0,0,0];
        a=[1,ap];
        for i=1:N
            w(i)=2*pi*(i-1)/N;
            fenzi1(1)=0;
            fenzi1(2)=0;
            fenmu1(1)=0;
            fenmu1(2)=0;
            for j=1:(p+1)
                fenzi1(1)=fenzi1(1)+b1(j)*cos(-j*w(i));
                fenzi1(2)=fenzi1(2)+b1(j)*sin(-j*w(i));
                fenmu1(1)=fenmu1(1)+a(j)*cos(-j*w(i));
                fenmu1(2)=fenmu1(2)+a(j)*sin(-j*w(i));
            end
            tmp=fenmu1(1)*fenmu1(1)+fenmu1(2)*fenmu1(2);
            H_z(1)=(fenzi1(1)*fenmu1(1)+fenzi1(2)*fenmu1(2))/tmp;
            H_z(2)=(fenzi1(2)*fenmu1(1)-fenzi1(1)*fenmu1(2))/tmp;
            HzAmp(i)=sqrt(H_z(1)*H_z(1)+H_z(2)*H_z(2)); 
            Sx(i)= HzAmp(i)^2;
        end
        db=20*log10(Sx);
        band(1:12)=db(1:12);
        bandpower(l)=mean(band);
    end
    Yi(t)=sum(bandpower);
end
Y_base=max(Yi);
t99=csvread('t99.csv');
tp=0;
for v=1:20
    for r=1:10
        wn=0.0001*randn(1,N);
        xn=t99((N*(10*(v-1)+r-1+200)+1):N*(10*(v-1)+r+200),:)'+wn;
        rou=zeros(p,p);
        e(1,:)=xn;
        b(1,:)=xn;
        sita(1,1)=mean(xn.^2);
        for m=1:p
            if(m==1)
                km=-2*(sum(e(m,m+1:N).*(b(m,(m:N-1)))))/sum(e(m,m+1:N).^2+b(m,m:N-1).^2);
            else
                km=-2*(sum(e(m,m:N).*(b(m,m-1:N-1))))/sum(e(m,m:N).^2+b(m,m-1:N-1).^2);
            end
            for i=1:m-1 
                rou(m,i)=rou(m-1,i)+km*rou(m-1,m-i);
            end
            rou(m,m)=km;
            for i=1:N
                if(i==1)
                    e(m+1,i)=e(m,i);
                    b(m+1,i)=km*e(m,i);
                else
                    e(m+1,i)=e(m,i)+km*b(m,i-1);
                    b(m+1,i)=b(m,i-1)+km*e(m,i);
                end
            end
            sita(m+1)=(1-km.^2)*sita(m);
        end
        G2=sita(p);
        ap=zeros(1,p);
        for k=1:p
            ap(k)=rou(p,k);
        end
        b1=[G2^0.5,0,0,0,0];
        a=[1,ap];
        for i=1:N
            w(i)=2*pi*(i-1)/N;
            fenzi1(1)=0;
            fenzi1(2)=0;
            fenmu1(1)=0;
            fenmu1(2)=0;
            for j=1:(p+1) 
                fenzi1(1)=fenzi1(1)+b1(j)*cos(-j*w(i));
                fenzi1(2)=fenzi1(2)+b1(j)*sin(-j*w(i));
                fenmu1(1)=fenmu1(1)+a(j)*cos(-j*w(i));
                fenmu1(2)=fenmu1(2)+a(j)*sin(-j*w(i));
            end
            tmp=fenmu1(1)*fenmu1(1)+fenmu1(2)*fenmu1(2);
            H_z(1)=(fenzi1(1)*fenmu1(1)+fenzi1(2)*fenmu1(2))/tmp;
            H_z(2)=(fenzi1(2)*fenmu1(1)-fenzi1(1)*fenmu1(2))/tmp;
            HzAmp(i)=sqrt(H_z(1)*H_z(1)+H_z(2)*H_z(2)); 
            Sx(i)= HzAmp(i)^2;
        end
        db=20*log10(Sx);
        %figure(1);plot(w/(2*pi),db);%绘制功率谱
        %title('功率谱密度');
        band(1:12)=db(1:12);
        bandpower(r)=mean(band);
    end
    Y=sum(bandpower)-Y_base;
    if (Y>5)
        tp=tp+1;
    end
    if(Y<=5)
    end
end
if (tp>=3)
    str=sprintf('发生拉弧')
end
if (tp<3)
    str=sprintf('未发生拉弧')
end