function [ct] = tumorSimulation(STAGE, GENE_LVL)

tumorSimulation(3, [10.4, 2.7, 3.5, 1.2, 7.4])



threshold = 10;
if STAGE == 1
    endct = 21;
elseif STAGE == 2
    endct = 34;
elseif STAGE == 3
    endct = 48;
elseif STAGE == 4
    endct = 65;
elseif STAGE == 5
    endct = 85;
end

Thold = [0.5, 0.5, 0.5, 0.5, 0.5];
GeneHot = [];
for i=1:5
    if GENE_LVL(i)>Thold(i)
        ACTIVATED = true;
    else
        ACTIVATED = false;
    end
    GeneHot = [GeneHot, ACTIVATED];
end
geneScore=sum(GeneHot);
m_prob=[0.2, 0.4, 0.6, 0.75, 0.85, 0.95]; % The Prob that needs to be checked 
if geneScore == 0
    M_prob=m_prob(1);
elseif geneScore == 1
    M_prob=m_prob(2);
elseif geneScore == 2
    M_prob=m_prob(3);
elseif geneScore == 3
    M_prob=m_prob(4);
elseif geneScore == 4
    M_prob=m_prob(5);
elseif geneScore == 5
    M_prob=m_prob(6);
end

%original probablity.m 
b=3*10^12/10000;
a=0.00286;
t=linspace(1,4000,100000);
N=b.^(1-exp(-a*t));
Dif=diff(N);
norm_Dif=Dif/sum(Dif);

x1=[0:0.01:1.4];
norm1=normpdf(x1,1.4,1.2);

x2=[1.4:0.01:3];
norm2=normpdf(x2,1.4,0.5);
norm2=norm1(end)*(norm2/norm2(1));

X=[x1 x2];
Norm=[norm1 norm2];
Norm(303:305)=0;

IM = ones(100);
IM(10,15)=0.1;
IM(10,25)=0.1;
IM(10,30)=0.1;
IM(10,70)=0.1;
IM(10,80)=0.1;
IM(10,90)=0.1;
IM(20,30)=0.1;
IM(20,80)=0.1;
IM(30,20)=0.1;
IM(30,70)=0.1;
IM(40,30)=0.1;
IM(40,90)=0.1;
IM(50,20)=0.1;
IM(50,90)=0.1;
IM(70,20)=0.1;
IM(70,90)=0.1;
IM(90,20)=0.1;
IM(90,30)=0.1;

% Creating Square border for LN
for i=2:99 
    for j=2:99
        if IM(i,j)== 0.1
            IM(i-1,j) = 0.2;
            IM(i+1,j) = 0.2;
            IM(i,j+1) = 0.2;
            IM(i, j-1) = 0.2;
            IM(i-1, j+1) = 0.2;
            IM(i-1, j-1) = 0.2;
            IM(i+1, j+1) = 0.2;
            IM(i+1, j-1) = 0.2;
        end
    end
end

Original_IM = IM;

% Original Tumor location
tx = randi([5,95]);
ty = randi([5,95]);

N0=1;
ct=N0;

N0=1;
ct=N0;

Xgrowth=[tx];
Ygrowth=[ty];

IM(tx,ty)=0;

thold=1;
for i=1:endct
    i;
    AN(i,1)=i-1;
    AN(i,2)=ct;
    int = round(ct/(10^3));
    if int>304
        int=303;
    end

    Prob = Norm(int+1);
    for j=1:ct 
        if rand<Prob
            ct=ct+1;
            thold=thold+1;
        else
            ct=ct;
            thold=thold;
        end
    end
    ct;
    if mod(i,5) == 0
        for m=2:99
            for n=2:99
                if IM(m,n) == 0
                    prob=rand;
                    if prob<0.25
                        IM(m+1,n) = 0;
                        IM(m-1,n) = 0;
                    elseif 0.25<prob<0.5
                        IM(m,n+1) = 0;
                        IM(m,n-1) = 0;
                    elseif 0.5<prob<0.75
                        IM(m+1,n+1) = 0;
                        IM(m+1,n-1) = 0;
                    elseif 0.75<prob<1
                        IM(m-1,n+1) = 0;
                        IM(m-1,n-1) = 0;
                    end
                end
            end
        end
        thold=0;

        %print((num2str(i)),'-djpeg')
    end
%         if thold > threshold
%             for m=2:99
%                 for n=2:99
%                     if IM(m,n) == 0
%                         prob=rand;
%                         if prob<0.25
%                             IM(m+1,n) = 0;
%                         elseif 0.25<prob<0.5
%                             IM(m-1,n) = 0;
%                         elseif 0.5<prob<0.75
%                             IM(m,n+1) = 0;
%                         elseif 0.75<prob<1
%                             IM(m,n-1) = 0;
%                         end
%                     end
%                 end
%             end
%             thold=0;
%                    
%             %print((num2str(i)),'-djpeg')
%         end
    %if random number is above certain threshold, invade LYMPH NODE! 


    CT=0;
    for c=1:100
        for v=1:100
            if Original_IM(c,v)==0.1 && IM(c,v)~=0.1 
                rm_prob = rand;
                if rm_prob<M_prob
                    IM(c,v)=100;
                end
            end           
        end
    end

%         if thold ~= 0
%             
%             CT=0;
%             for c=1:100
%                 for v=1:100
%                     if Original_IM(c,v)==0.1 && IM(c,v)~=0.1 
%                         rm_prob = rand;
%                         if rm_prob<M_prob
%                             IM(c,v)=100;
%                         end
%                     end           
%                 end
%             end
%         end


end
ct=0;
for i=1:100
    for j=1:100
        if IM(i,j)== 100
            ct=ct+1;
        end
    end           
end
ct;
end

% function [ct,endct] = tumorSimulation(STAGE, GENE_LVL)
%     threshold = 10;
%     if STAGE == 1
%         endct = 21;
%     elseif STAGE == 2
%         endct = 34;
%     elseif STAGE == 3
%         endct = 48;
%     elseif STAGE == 4
%         endct = 65;
%     elseif STAGE == 5
%         endct = 120;
%     end
%     
%     Thold = [0.5, 0.5, 0.5, 0.5, 0.5];
%     GeneHot = [];
%     for i=1:5
%         if GENE_LVL(i)>Thold(i)
%             ACTIVATED = true;
%         else
%             ACTIVATED = false;
%         end
%         GeneHot = [GeneHot, ACTIVATED];
%     end
%     geneScore=sum(GeneHot);
%     m_prob=[0.5, 0.6, 0.7, 0.75, 0.8, 0.9]; % The Prob that needs to be checked 
%     if geneScore == 0
%         M_prob=m_prob(1);
%     elseif geneScore == 1
%         M_prob=m_prob(2);
%     elseif geneScore == 2
%         M_prob=m_prob(3);
%     elseif geneScore == 3
%         M_prob=m_prob(4);
%     elseif geneScore == 4
%         M_prob=m_prob(5);
%     elseif geneScore == 5
%         M_prob=m_prob(6);
%     end
%     
%     %original probablity.m 
%     b=3*10^12/10000;
%     a=0.00286;
%     t=linspace(1,4000,100000);
%     N=b.^(1-exp(-a*t));
%     Dif=diff(N);
%     norm_Dif=Dif/sum(Dif);
%     
%     x1=[0:0.01:1.4];
%     norm1=normpdf(x1,1.4,1.2);
% 
%     x2=[1.4:0.01:3];
%     norm2=normpdf(x2,1.4,0.5);
%     norm2=norm1(end)*(norm2/norm2(1));
% 
%     X=[x1 x2];
%     Norm=[norm1 norm2];
%     Norm(303:305)=0;
%     
%     IM = ones(100);
%     IM(10,15)=0.1;
%     IM(10,25)=0.1;
%     IM(10,30)=0.1;
%     IM(10,70)=0.1;
%     IM(10,80)=0.1;
%     IM(10,90)=0.1;
%     IM(20,30)=0.1;
%     IM(20,80)=0.1;
%     IM(30,20)=0.1;
%     IM(30,70)=0.1;
%     IM(40,30)=0.1;
%     IM(40,90)=0.1;
%     IM(50,20)=0.1;
%     IM(50,90)=0.1;
%     IM(70,20)=0.1;
%     IM(70,90)=0.1;
%     IM(90,20)=0.1;
%     IM(90,30)=0.1;
% 
%     % Creating Square border for LN
%     for i=2:99 
%         for j=2:99
%             if IM(i,j)== 0.1
%                 IM(i-1,j) = 0.2;
%                 IM(i+1,j) = 0.2;
%                 IM(i,j+1) = 0.2;
%                 IM(i, j-1) = 0.2;
%                 IM(i-1, j+1) = 0.2;
%                 IM(i-1, j-1) = 0.2;
%                 IM(i+1, j+1) = 0.2;
%                 IM(i+1, j-1) = 0.2;
%             end
%         end
%     end
% 
%     Original_IM = IM;
% 
%     % Original Tumor location
%     tx = randi([20,80]);
%     ty = randi([20,80]);
% 
%     N0=1;
%     ct=N0;
% 
%     N0=1;
%     ct=N0;
% 
%     Xgrowth=[tx];
%     Ygrowth=[ty];
% 
%     IM(tx,ty)=0;
% 
%     thold=1;
%     for i=1:endct
%         i;
%         AN(i,1)=i-1;
%         AN(i,2)=ct;
%         int = round(ct/(10^3));
%         if int>304
%             int=303;
%         end
% 
%         Prob = Norm(int+1);
%         for j=1:ct 
%             if rand<Prob
%                 ct=ct+1;
%                 thold=thold+1;
%             else
%                 ct=ct;
%                 thold=thold;
%             end
%         end
%         ct;
%         if thold > threshold
%             for m=2:99
%                 for n=2:99
%                     if IM(m,n) == 0
%                         prob=rand;
%                         if prob<0.25
%                             IM(m+1,n) = 0;
%                             IM(m-1,n) = 0;
%                         elseif 0.25<prob<0.5
%                             IM(m,n+1) = 0;
%                             IM(m,n-1) = 0;
%                         elseif 0.5<prob<0.75
%                             IM(m+1,n+1) = 0;
%                             IM(m+1,n-1) = 0;
%                         elseif 0.75<prob<1
%                             IM(m-1,n+1) = 0;
%                             IM(m-1,n-1) = 0;
%                         end
%                     end
%                 end
%             end
%             thold=0;
%                    
%             %print((num2str(i)),'-djpeg')
%         end
% %         if thold > threshold
% %             for m=2:99
% %                 for n=2:99
% %                     if IM(m,n) == 0
% %                         prob=rand;
% %                         if prob<0.25
% %                             IM(m+1,n) = 0;
% %                         elseif 0.25<prob<0.5
% %                             IM(m-1,n) = 0;
% %                         elseif 0.5<prob<0.75
% %                             IM(m,n+1) = 0;
% %                         elseif 0.75<prob<1
% %                             IM(m,n-1) = 0;
% %                         end
% %                     end
% %                 end
% %             end
% %             thold=0;
% %                    
% %             %print((num2str(i)),'-djpeg')
% %         end
%         %if random number is above certain threshold, invade LYMPH NODE! 
%         
% 
%         CT=0;
%         for c=1:100
%             for v=1:100
%                 if Original_IM(c,v)==0.1 && IM(c,v)~=0.1 
%                     rm_prob = rand;
%                     if rm_prob<M_prob
%                         IM(c,v)=100;
%                     end
%                 end           
%             end
%         end
% 
% %         if thold ~= 0
% %             
% %             CT=0;
% %             for c=1:100
% %                 for v=1:100
% %                     if Original_IM(c,v)==0.1 && IM(c,v)~=0.1 
% %                         rm_prob = rand;
% %                         if rm_prob<M_prob
% %                             IM(c,v)=100;
% %                         end
% %                     end           
% %                 end
% %             end
% %         end
% 
% 
%     end
%     ct=0;
%     for i=1:100
%         for j=1:100
%             if IM(i,j)== 100
%                 ct=ct+1;
%             end
%         end           
%     end
%     ct;
% end