function [y1t,y2t,y3t] = ScriptMatrix(y, y2, y3)

k=1;
yt = zeros (4,5);
for i=1:5:20   
    j=1;
    fprintf ("%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\n", y(i),y(i+1),y(i+2),y(i+3),y(i+4));
    yt(k,j) = y(i); %[y(i) y(i+1) y(i+2) y(i+3) y(i+4)];
    yt(k,j+1) = y(i+1);
    yt(k,j+2) = y(i+2);
    yt(k,j+3) = y(i+3);
    yt(k,j+4) = y(i+4);
    k = k +1;
end
y1t = yt'; 
fprintf ("=====================================\n");

k=1;
yt = zeros (4,5);
for i=1:5:20   
    j=1;
    fprintf ("%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\n", y2(i),y2(i+1),y2(i+2),y2(i+3),y2(i+4));
    yt(k,j) = y2(i); %[y(i) y(i+1) y(i+2) y(i+3) y(i+4)];
    yt(k,j+1) = y2(i+1);
    yt(k,j+2) = y2(i+2);
    yt(k,j+3) = y2(i+3);
    yt(k,j+4) = y2(i+4);
    k = k +1;
end
fprintf ("=====================================\n");
y2t = yt'; 
k=1;
yt = zeros (4,5);
for i=1:5:20   
    j=1;
    fprintf ("%7.5f\t%7.5f\t%7.5f\t%7.5f\t%7.5f\n", y3(i),y3(i+1),y3(i+2),y3(i+3),y3(i+4));
    yt(k,j) = y3(i); %[y(i) y(i+1) y(i+2) y(i+3) y(i+4)];
    yt(k,j+1) = y3(i+1);
    yt(k,j+2) = y3(i+2);
    yt(k,j+3) = y3(i+3);
    yt(k,j+4) = y3(i+4);
    k = k +1;
end
fprintf ("=====================================\n");
y3t = yt'; 
%-------------------------------------------------------%

end
