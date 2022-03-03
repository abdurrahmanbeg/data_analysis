clear all
y1 = [
1.1686
1.5645
1.4983
1.411
1.1738
1.7631
1.83
1.9448
2.0379
1.9024
1.8689
1.8455
2.0313
2.1698
2.0842
1.9124
1.8453
2.0754
2.2264
2.1634];

y2 = [
1.7583
1.7794
1.706
1.6897
1.3691
2.0883
2.6321
2.567
2.4516
2.2608
2.2407
2.7293
2.5966
2.6692
2.539
2.3073
2.7743
2.6216
2.759
2.6646
];

y3 = [
 1.763
1.8416
1.8361
1.8312
1.403
4.5701
6.2028
5.8989
5.8177
5.1758
5.5439
8.1817
7.1365
6.4261
7.1545
6.0315
9.4456
7.3666
6.108
7.2779   
];

y(:,:,1) = [    1.1686    1.7631    1.8689    1.9124;
    1.5645    1.8300    1.8455    1.8453;
    1.4983    1.9448    2.0313    2.0754;
    1.4110    2.0379    2.1698    2.2264;
    1.1738    1.9024    2.0842    2.1634]

y(:,:,2) = [    1.7583    2.0883    2.2407    2.3073;
    1.7794    2.6321    2.7293    2.7743;
    1.7060    2.5670    2.5966    2.6216;
    1.6897    2.4516    2.6692    2.7590;
    1.3691    2.2608    2.5390    2.6646]

y(:,:,3) = [

    1.7630    4.5701    5.5439    6.0315;
    1.8416    6.2028    8.1817    9.4456;
    1.8361    5.8989    7.1365    7.3666;
    1.8312    5.8177    6.4261    6.1080;
    1.4030    5.1758    7.1545    7.2779]

%[y(:,:,1),y(:,:,2),y(:,:,3)] = ScriptMatrix(y1, y2, y3)

[b,a,r] = size(y);
u= mean(mean(mean(y)));

for i=1:b
bEffect(i) = mean(mean(y(i,:,:))) - u;
end
for i=1:a
aEffect(i) = mean(mean(y(:,i,:))) - u;
end
for i=1:r
rEffect(i) = mean(mean(y(:,:,i))) - u;
end

SSY = sum(sum(sum(y.^2)));
SS0 = a*b*r*u^2;
SST = SSY - SS0;
SSA = b*r*sum(aEffect.^2);
SSB = a*r*sum(bEffect.^2);
SSR = a*b*sum(rEffect.^2);
MF = (SSA+SSB+SSR)/SST;
for i=1:b
    for j=1:a
        intAB(i,j) = mean(y(i,j,:)) - (u + bEffect(i) + aEffect(j)); 
    end
end    
SSAB = r*sum(sum(intAB.^2));
for i=1:a
    for j=1:r
        intAR(i,j) = mean(y(:,i,j)) - (u + aEffect(i) + rEffect(j)); 
    end
end    
SSAR = b*sum(sum(intAR.^2));
for i=1:b
    for j=1:r
        intBR(i,j) = mean(y(i,:,j)) - (u + bEffect(i) + rEffect(j)); 
    end
end    
SSBR = a*sum(sum(intBR.^2));
FE = (SSAB+SSBR+SSAR)/SST;
%SecondOrder ABR
for i=1:b
    for j=1:a
        for k=1:r
             intABR(i,j,k) = mean(y(i,j,k)) - (mean(mean(y(i,:,:))) + mean(mean(y(:,j,:))) + mean(mean(y(:,:,k))) ) + u ...
                   - (intAB(i,j)+intAR(j,k)+intBR(i,k)) + u; 
        end
    end   
end
SSABR = sum(sum(sum(intABR.^2)));
SE = SSABR/SST;

fprintf("================================================================================\n");
%-------------------------ANOVA----------------------------%
fprintf("--------------------------------------------------------------------------------\n");
fprintf("|\tComponent\t\tSumofSquares\t\t%%Variation\t\tDF\t\tMeanSquare\t|\n");
fprintf("--------------------------------------------------------------------------------\n");
fprintf("|\t%10s\t\t\t%7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'y', SSY, '', b*a*r, '');
fprintf("|\t%10s\t\t\t%7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'ybar', SS0, '', 1, '');
fprintf("|\t%10s\t\t\t%7.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'y-ybar', SST, 100, a*b*r-1, '');
fprintf("|\t%10s\t\t%7.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t|\n",...
            '[Main Effects]', SSA+SSB+SSR, MF*100, (b-1)+(a-1)+(r-1), (SSA+SSB+SSR)/((b-1)+(a-1)+(r-1)));      
fprintf("|\t%10s\t\t%7.4f\t\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'No. of Stations', SSA, '', (a-1), '');  
fprintf("|\t%10s\t\t%7.4f\t\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'MMAC Protocols', SSB, '', (b-1), '');  
fprintf("|\t%10s\t\t%7.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'No. of Channels', SSR, '', (r-1), '');  
fprintf("|%10s %7.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            '  [FirstOrderInteract]', SSAB+SSBR+SSAR, FE*100, (a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1), (SSAB+SSBR+SSAR)/((a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1)));          
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Protocols', SSAB, '', (a-1)*(b-1), '');
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Channels', SSAR, '', (a-1)*(r-1), '');        
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Protocols*Channels', SSBR, '', (b-1)*(r-1), ''); 
fprintf("|%10s %7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            '  [SecondOrderInteract]', SSABR, SE*100, (a-1)*(b-1)*(r-1), SSABR/((a-1)*(b-1)*(r-1)));
fprintf("|\t%10s %7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stats*Proto*Channels', SSABR, '', (a-1)*(b-1)*(r-1),'');
fprintf("================================================================================\n");
%-----------------ERRORS--And--QQPlot----------------------%
index=1; 
errors=zeros(1,60);
for i=1:b
    for j=1:a
        for k=1:r
%                 errors(index) = mean(Ys(i,j,k,l)) - (mainEffects(1,i) + mainEffects(2,j) + mainEffects(3,k) + mainEffects(4,l) + 4*u) ...
%                     - interactionsDM(j,l) + 3*u;
                errors(index) =  mean(y(i,j,k)) - (mean(mean(y(i,:,:))) + mean(mean(y(:,j,:))) + mean(mean(y(:,:,k))) ) + 2*u ...
                      - intAR(j,k);
                index=index+1;
        end
    end
end
[x,yt]=QQplot_normal(errors); 
plot(x,sort(errors),'*',x,yt); 
%axis([-4 4 -1.5 2])
% xticks([-4:2:4]); yticks([-0.5:0.25:0.5]);
xlabel('Normal quantile'); ylabel('Residual quantile');