function h=Ricrnd(a,sigmah)
   
    u=randn();
    u=(u*sqrt(sigmah));
    u=u+a;
    u=u.^2;
    v=randn();
    v=(v*sqrt(sigmah)).^2;
    
    h=(u+v); %|h|^2
   
end

