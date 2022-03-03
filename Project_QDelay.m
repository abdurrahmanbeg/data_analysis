clear all
y1 = [
31.9649
17.2823
18.4422
20.1416
26.0319
113.2934
103.3374
96.6363
91.8245
99.0138
194.5549
191.5202
172.9494
161.3628
168.3504
275.8265
280.6731
248.0905
230.6428
237.6968];

y2 = [
23.0141
14.1096
15.102
15.3325
21.033
95.892
69.107
71.0692
74.8223
81.8951
162.4923
126.5819
133.3925
129.5107
136.6349
228.8616
183.6705
194.6338
184.5011
191.302
];

y3 = [22.9672
13.3236
13.3912
13.4517
20.306
44.6234
24.1418
25.8433
26.337
30.7154
66.3773
36.2226
42.8449
48.5763
42.7147
88.2319
47.5783
63.5375
78.4834
64.4173 
];

[y(:,:,1),y(:,:,2),y(:,:,3)] = ScriptMatrix(y1,y2, y3)

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
fprintf("|%10s %12.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t|\n",...
            '  [FirstOrderInteract]', SSAB+SSBR+SSAR, FE*100, (a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1), (SSAB+SSBR+SSAR)/((a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1)));          
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Protocols', SSAB, '', (a-1)*(b-1), '');
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Stations*Channels', SSAR, '', (a-1)*(r-1), '');        
fprintf("|\t%10s\t%12.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Protocols*Channels', SSBR, '', (b-1)*(r-1), ''); 
fprintf("|%10s %12.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t|\n",...
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