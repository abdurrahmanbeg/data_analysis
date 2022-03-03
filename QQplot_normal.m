function [x,y1] = QQplot_normal(y)
n=length(y);
y=sort(y);
qi=([1:n]-0.5)/n;
x=norminv(qi);
b1=(sum(x*y') - n*mean(x)*mean(y))/(sum(x.^2)-n*mean(x)^2);
b0= mean(y) - b1*mean(x);
y1=b0 + b1*x;
end