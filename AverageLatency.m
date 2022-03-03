clear all

y(:,:,1) = [ 1.303035683	1.714659321	3.072276295	4.655262929;
1.303081924	1.736031626	3.09066212	4.688298742;
1.3029963	1.735993572	3.047199397	4.594438561;
1.303019909	1.736010623	3.042905231	3.555321345;
1.303019484	1.736001338	3.099192105	4.498108279]



y(:,:,2) = [  1.259252466	1.619351623	2.881401682	5.131500625;
1.259250413	1.57473902	3.371134109	7.060530904;
1.25916389	1.69216611	3.114062612	4.703607516;
1.303019909	1.736010623	3.042905231	3.555321345;
1.303019484	1.736001338	3.099192105	4.498108279]

y(:,:,3) = [1.303035683	1.714659321	3.072276295	4.655262929;
1.303081924	1.736031626	3.09066212	4.688298742;
1.3029963	1.735993572	3.047199397	4.594438561;
1.303019909	1.736010623	3.042905231	3.555321345;
1.303019484	1.736001338	3.099192105	4.498108279]

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
fprintf("|\tComponent\tSumofSquares\t%%Variation\tDF\tMeanSquare\t|\n");
fprintf("--------------------------------------------------------------------------------\n");
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'y', SSY, '', b*a*r, '');
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'ybar', SS0, '', 1, '');
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'y-ybar', SST, 100, a*b*r-1, '');
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            '[Main Effects]', SSA+SSB+SSR, MF*100, (b-1)+(a-1)+(r-1), (SSA+SSB+SSR)/((b-1)+(a-1)+(r-1)));      
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'Density (D)', SSA, (SSA/SST)*100, (a-1), '');  
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'Type Dist.(T)', SSB, (SSB/SST)*100, (b-1), '');  
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'R. Protocols (RP)', SSR, (SSR/SST)*100, (r-1), '');  
fprintf("|\t%10s %7.4f\t%.1f\t%d\t%.1f\t|\n",...
            '  [FirstOrderInteract]', SSAB+SSBR+SSAR, FE*100, (a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1), (SSAB+SSBR+SSAR)/((a-1)*(b-1)+(a-1)*(r-1)+(b-1)*(r-1)));          
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'Density*TypeDist', SSAB, '', (a-1)*(b-1), '');
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'Density*Protocols', SSAR, '', (a-1)*(r-1), '');        
fprintf("|\t%10s\t%7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'TypeDist*Protocols', SSBR, '', (b-1)*(r-1), ''); 
fprintf("|\t%10s %7.4f\t%.1f\t%d\t%.1f\t|\n",...
            '  [SecondOrderInteract]', SSABR, SE*100, (a-1)*(b-1)*(r-1), SSABR/((a-1)*(b-1)*(r-1)));
fprintf("|\t%10s %7.4f\t%.1f\t%d\t%.1f\t|\n",...
            'Density*TypeDist*Protocols', SSABR, '', (a-1)*(b-1)*(r-1),'');
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