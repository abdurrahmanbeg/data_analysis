clear all

y(:,:,1) = [ 73.6666667	84.333333;
69	79;
73.6666667	84;
70	79;
73.6666667	84.333333]



y(:,:,2) = [  121	124.33333;
114	99.333333;
112.66667	123.66667;
113	115.33333;
111.66667	124.66667]

y(:,:,3) = [121.66667	122.666667;
112	114.333333;
121	123.333333;
119.33333	121;
121.33333	126]

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
            'Density (D)', SSA, (SSA/SST)*100, (a-1), '');  
fprintf("|\t%10s\t\t%7.4f\t\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'Type Dist.(T)', SSB, (SSB/SST)*100, (b-1), '');  
fprintf("|\t%10s\t\t%7.4f\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Routes (R)', SSR, (SSR/SST)*100, (r-1), '');  
fprintf("|%10s %7.4f\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            '  [FirstOrderInteract]', SSAB+SSBR+SSAR, FE*100, (a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1), (SSAB+SSBR+SSAR)/((a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1)));          
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Density*TypeDist', SSAB, '', (a-1)*(b-1), '');
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'Density*Routes', SSAR, '', (a-1)*(r-1), '');        
fprintf("|\t%10s\t%7.4f\t\t\t\t\t%.1f\t\t%d\t\t\t%.1f\t\t|\n",...
            'TypeDist*Routes', SSBR, '', (b-1)*(r-1), ''); 
fprintf("|%10s %7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            '  [SecondOrderInteract]', SSABR, SE*100, (a-1)*(b-1)*(r-1), SSABR/((a-1)*(b-1)*(r-1)));
fprintf("|\t%10s %7.4f\t\t\t%.1f\t\t\t%d\t\t\t%.1f\t\t|\n",...
            'Density*TypeDist*Routes', SSABR, '', (a-1)*(b-1)*(r-1),'');
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