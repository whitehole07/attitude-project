function [Mat_g, Mat_h, Mat_gsv,  Mat_hsv]=ghnorm()
[m] = readtable( 'igrfcoefficients_2020-2025.xlsx');
col1 = m.Var1; 
col1 = col1(2:end, :);
gh=col1;
col2 = m.Var2;
col3 = m.Var3;
col4 = m.Var4;
col5 = m.Var5;
mat = [col2, col3, col4, col5];
mat = mat (2:end, :);

n = mat(:, 1);
m = mat (:, 2);
val = mat(:, 3);
sv = mat(:, 4);
N=max(mat (:, 2));

g=zeros(N,N+1);
h=zeros(N,N+1);
hsv=zeros(N,N+1);
gsv=zeros(N,N+1);

for x=1:length(gh)
if strcmp(gh(x),'g')
g(n(x),m(x)+1) = val(x);
gsv(n(x),m(x)+1) = sv(x);
else
h(n(x),m(x)+1) = val(x);
hsv(n(x),m(x)+1) = sv(x);
end
end
count=1;
% starts quasi-normalization of Gram - Schmidt
S = zeros(N,N+1);

for ii=1:N
    for jj=0:ii
        if jj>1
            S(ii,jj+1) = S(ii,jj)*((ii-jj+1)/(ii+jj))^0.5;
        elseif jj>0
            S(ii,jj+1) = S(ii,jj)*(2*(ii-jj+1)/(ii+jj))^0.5;
        elseif ii==1
            S(ii,1) = 1;
        else
            S(ii,1) = S(ii-1,1)*(2*ii-1)/(ii);
        end
        
%         gS(count,1) = n; 
%         gS(count,2)=m;
%         gS(count,3)=g(n,m+1)*S(n,m+1); 
%         gS(count,4)=gsv(n,m+1)*S(n,m+1);
% 
%         
% 
%         hS(count,1) = n; 
%         hS(count,2)=m;
%         hS(count,3)=h(n,m+1)*S(n,m+1); 
%         hS(count,4)=hsv(n,m+1)*S(n,m+1);


                count=count+1;
    end
end

 Mat_g= g .* S * 1e-09;
 Mat_h = h .* S * 1e-09;

 Mat_gsv= gsv .* S*  1e-09;
 Mat_hsv = hsv .* S*  1e-09;


% 
% dlmwrite('igrfSg.txt',gS,'\t')
% dlmwrite('igrfSh.txt',hS,'\t')
end