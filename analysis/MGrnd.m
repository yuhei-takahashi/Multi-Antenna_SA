function h=MGrnd(a,b,c)
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
