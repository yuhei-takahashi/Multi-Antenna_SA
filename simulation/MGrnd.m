%{
function h = MGrnd(a, b, c, n)
    %MATLABの関数
    N = length(a); % ガンマ成分の数
    h = zeros(1, n); % 出力の初期化
    sum_weights = sum(a .* gamma(b) ./ (c.^b)); % 重みの合計を計算
    
    for j = 1:n
        for i = 1:N
            w = a(i) * gamma(b(i)) / (c(i)^b(i));
            h(j) = h(j) + w * gamrnd(b(i), 1 / c(i));
        end
    end
    h = h / sum_weights; % 重みの合計で正規化
end
%}



%{
function h=MGrnd(a,b,c)
%MATLABの関数
h=0;
n=length(a);
sum=0;
for i=1:n
    w=a(i)*gamma(b(i))/(c(i)^b(i));
    sum=sum+w;
    h=h+w*gamrnd(b(i),1/c(i));
end
%disp(sum);
%disp(a)
%disp(b)
%disp(c)
end
%}


%{
function h=MGrnd(a,b,c)
    w=a.*gamma(b)./(c.^b);
    %w=[0,1];
    gam=gamrnd(b,1./c);
    h=sum(w.*gam);
    %disp(gam)
    %disp(w)
    %h
end
%}

function h=MGrnd(a,b,c)
%MATLABの関数
h=0;
w=0;
n=length(a);
rand_number=rand(1);
for i=1:n
    w=w+a(i)*gamma(b(i))/(c(i)^b(i));
    if rand_number<w
        h=gamrnd(b(i),1/c(i));
        break
    end
end
end
