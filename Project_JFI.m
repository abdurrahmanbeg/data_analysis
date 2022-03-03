clear all
y1 = [0.8069567
0.7977091
0.750852
0.7110085
0.7368779
0.7904404
0.5488091
0.6478886
0.6146575
0.647567
0.7898148
0.5298515
0.6393644
0.6073565
0.6386571
0.7898233
0.5383611
0.6371013
0.6011949
0.633528];

y2 = [
0.9934713
0.7892591
0.7667958
0.759892
0.7657992
0.8143268
0.6522625
0.6567529
0.6215822
0.654103
0.8274237
0.5710262
0.6474815
0.6101001
0.6433079
0.8240555
0.542972
0.6434413
0.6094636
0.639345
];

y3 = [0.9936895
0.7787996
0.7795374
0.7818171
0.7701067
0.9884716
0.6937148
0.6867246
0.6857782
0.6876925
0.9873182
0.6819697
0.6759531
0.6630887
0.6786368
0.9875801
0.675866
0.6714494
0.6538152
0.6700614
];

[y(:,:,1),y(:,:,2),y(:,:,3)] = ScriptMatrix(y1,y2,y3);

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
fprintf("|\t%10s\t\t\t%12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'y', SSY, '', b*a*r, '');
fprintf("|\t%10s\t\t\t%12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'ybar', SS0, '', 1, '');
fprintf("|\t%10s\t\t\t%12.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'y-ybar', SST, 100, a*b*r-1, '');
fprintf("|\t%10s\t\t%12.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t|\n",...
            '[Main Effects]', SSA+SSB+SSR, MF*100, (b-1)+(a-1)+(r-1), (SSA+SSB+SSR)/((b-1)+(a-1)+(r-1)));      
fprintf("|\t%10s\t\t%12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'No. of Stations', SSA, '', (a-1), '');  
fprintf("|\t%10s\t\t%12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'MMAC Protocols', SSB, '', (b-1), '');  
fprintf("|\t%10s\t\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'No. of Channels', SSR, '', (r-1), '');  
fprintf("|%10s %12.4f\t\t\t%.1f\t\t%d\t\t\t%.5f\t|\n",...
            '  [FirstOrderInteract]', SSAB+SSBR+SSAR, FE*100, (a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1), (SSAB+SSBR+SSAR)/((a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1)));          
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Protocols', SSAB, '', (a-1)*(b-1), '');
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Channels', SSAR, '', (a-1)*(r-1), '');        
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Protocols*Channels', SSBR, '', (b-1)*(r-1), ''); 
fprintf("|%10s %12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.5f\t|\n",...
            '  [SecondOrderInteract]', SSABR, SE*100, (a-1)*(b-1)*(r-1), SSABR/((a-1)*(b-1)*(r-1)));
fprintf("|\t%10s %12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
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
                      - intBR(i,k);
                index=index+1;
        end
    end
end
[x,yt]=QQplot_normal(errors); 
plot(x,sort(errors),'*',x,yt); 
%axis([-4 4 -1.5 2])
% xticks([-4:2:4]); yticks([-0.5:0.25:0.5]);
xlabel('Normal quantile'); ylabel('Residual quantile');