global x; x= zeros(4,1);
global A; A= zeros(4,1);
global B; B= zeros(4,1);
global T; T= zeros(4,1);
global E; E= zeros(4,1);
global x_est; x_est= [];
global d_c;d_c = 10;


function bounds(E(i),d_c,T(i),i,R(i,i))
    if d_c < T(i)
        increment(i,m);
    else
        compute= (y(i)-E(i)+sqrt(d_c-T(i)))/R(i,i);
        A(i)=max(0,ceil(compute));
        B(i)=min(Q-1,floor(compute));
        x(i)=A(i)-1;
        natural_spanning(x(i),B(i));
    end
end


function natural_spanning(x(i),B(i))
    x(i)=x(i)+1;
    if x(i) <= B(i)
        decrement(R,x,T(i),E(i),y(i));      
    else
        increment(i,m);
    end
end



function increment(i,m)
    if i==m
        return
    else
        i=i+1;
        natural_spanning(x(i),B(i));
    end
end



function decrement(R,x,T(i),E(i),y(i))
    if i>1
        for j = i:m
            E(i-1) = E(i-1) + R(i-1,j)*x(j);
        end
        T(i-1) = T(i) + (y(i)-E(i)-R(i,i)*x(i))^2;
        i = i-1;
        bounds(y(i),E(i),d_c,T(i),i,R(i,i));
        
    else
        x_est = (x_est x);
        d' = T(1)+(y(1)-E(1)-R(1,1)*x(1))^2;
        if d'<d_c
            d_c = d';
            for L = 1:m
                B(L) = min(Q-1,floor(y(L)-E(L)+sqrt(d_c-T(L)))/R(L,L)));
            end
        natural_spanning(x(i),B(i));
        
        end
    
    end

end