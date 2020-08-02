%global R; R=[-4.38112687193725,0,0.995272654476084,1.30856086206953;0,-4.38112687193725,-1.30856086206953,0.995272654476084;0,0,5.07631532632225,-4.44089209850063e-16;0,0,0,-5.07631532632225];
%global d_c;d_c = 100;
%y =[7.73505848559650;1.27232920478601;-7.97395403822037;13.0184198562958];
%x_est = sphere(d_c,R,y)
function x_est = sphere_dec(d_c,R_mat,y);

    global x; x= zeros(4,1);
    global A; A= zeros(4,1);
    global B; B= zeros(4,1);
    global T; T= zeros(4,1);
    global E; E= zeros(4,1);
    global x_est; x_est= [];
    global m; m=4;
    global i; i=m;
    global Q, Q=4;

    
    bounds(x,A,B,T,E,d_c, m,i,Q,R_mat,y);
    function bounds(x,A,B,T,E,d_c,m,i,Q,R_mat,y)
        if d_c < T(i)
            increment(x,A,B,T,E,d_c, m,i,Q,R_mat,y);
        else
            compute1= (y(i)-E(i)-sqrt(d_c-T(i)))/R_mat(i,i);
            compute2= (y(i)-E(i)+sqrt(d_c-T(i)))/R_mat(i,i);
            A(i)=max(0,ceil(compute1));
            B(i)=min(Q-1,floor(compute2));
            x(i)=A(i)-1;
            natural_spanning(x,A,B,T,E,d_c, m,i,Q,R_mat,y);
        end
    end


    function natural_spanning(x,A,B,T,E,d_c, m,i,Q,R_mat,y)
        x(i)=x(i)+1;
        if x(i) <= B(i)
            decrement(x,A,B,T,E,d_c, m,i,Q,R_mat,y);      
        else
            increment(x,A,B,T,E,d_c, m,i,Q,R_mat,y);
        end
    end



    function increment(x,A,B,T,E,d_c, m,i,Q,R_mat,y)
        if i==m
            if isempty(x_est)
                d_c = 2*d_c;
                x_est = sphere_dec(d_c,R_mat,y);
            else 
                return
            end
        else
            i=i+1;
            natural_spanning(x,A,B,T,E,d_c, m,i,Q,R_mat,y);
        end
    end



    function  decrement(x,A,B,T,E,d_c, m,i,Q,R_mat,y)
        if i>1
            for j = i:m
                E(i-1) = E(i-1) + R_mat(i-1,j)*x(j);
            end
            T(i-1) = T(i) + (y(i)-E(i)-R_mat(i,i)*x(i))^2;
            i = i-1;
            bounds(x,A,B,T,E,d_c, m,i,Q,R_mat,y);

        else
            x_est = [x_est x];           
            d = T(1)+(y(1)-E(1)-R_mat(1,1)*x(1))^2;
            if d<d_c
                d_c = d;
                for L = 1:m
                    B(L) = min(Q-1,floor((y(L)-E(L)+sqrt(d_c-T(L)))/R_mat(L,L)));
                end
            natural_spanning(x,A,B,T,E,d_c, m,i,Q,R_mat,y);

            end

        end

    end
end