nodes = load ('breast2.nod');
elements = load('breast2.elm');
boundary_elms = load("breast2.bel");
boundary_cond = load('breast2.bcs');
data = load('breast2.data');

Ud = zeros(870*2,1);

for q = 1:length(data)
    node = data(q,1);
    Ud(node*2-1,1) = data(q,2);
    Ud(node*2,1) = data(q,3);
end

incTri = elements(:,2:4);
x = nodes(:,2);
y = nodes(:,3);

E = [36000; 36000; 36000];
z = 0.45;

J = zeros(2*870,3);

E_matrix = zeros(20,3);

E3E1 = zeros(20,1);
E2E1 = zeros(20,1);

L2 = zeros(20,1);

for iter = 1:20

    Ag = zeros(2*870,2*870);
    Bg = zeros(870*2,1);

    for L = 1:length(elements)

        
        [dpx, dpy, Area] = basisfunc3(L,x,y,incTri);
        
        for i = 1:3
            iix = 2*incTri(L,i) - 1;
            iiy = iix + 1;

            for j = 1:3
                jjx = 2*incTri(L,j) - 1;
                jjy = jjx + 1;

                K = E(elements(L,5));
                Ke = [((K*(1-z)/((1+z)*(1-2*z)))*(dpx(i)*dpx(j))*Area)+(K/(2*(1+z))*(dpy(i)*dpy(j))*Area)... 
                    (((K*z)/((1+z)*(1-2*z)))*(dpy(j)*dpx(i))*Area)+(K/(2*(1+z))*(dpx(j)*dpy(i))*Area);
                    ((K*z)/((1+z)*(1-2*z))*(dpx(j)*dpy(i))*Area)+(K/(2*(1+z))*(dpy(j)*dpx(i))*Area)... 
                    ((K*(1-z)/((1+z)*(1-2*z)))*(dpy(j)*dpy(i))*Area)+(K/(2*(1+z))*(dpx(j)*dpx(i))*Area)];
                
                Ag(iix,jjx) = Ag(iix,jjx) + Ke(1,1);
                Ag(iix,jjy) = Ag(iix,jjy) + Ke(1,2);

                Ag(iiy,jjx) = Ag(iiy,jjx) + Ke(2,1);
                Ag(iiy,jjy) = Ag(iiy,jjy) + Ke(2,2);
            end
        end

    end
    for N = 1:length(boundary_cond) % BCs
        x_bc = boundary_cond(N,5);
        y_bc = boundary_cond(N,7);
        node = boundary_cond(N,2);
        Ag(node*2 - 1,:) = 0;
        Ag(node*2,:) = 0;
        Ag(node*2 - 1,node*2 - 1) = 1;
        Ag(node*2,node*2) = 1;
        Bg(node*2 - 1) = x_bc;
        Bg(node*2) = y_bc;
    end

    Um = Ag\Bg;

    for materials = 1:3

        Kg = zeros(2*870,2*870);

        for L = 1:length(elements)

            element_type = elements(L,5);
            [dpx, dpy, Area] = basisfunc3(L,x,y,incTri);

            if element_type == materials

                for i = 1:3
                    iix = 2*incTri(L,i) - 1;
                    iiy = iix + 1;

                    for j = 1:3
                        jjx = 2*incTri(L,j) - 1;
                        jjy = jjx + 1;

                        Kd = [(((1-z)/((1+z)*(1-2*z)))*(dpx(i)*dpx(j))*Area)+(1/(2*(1+z))*(dpy(i)*dpy(j))*Area)... 
                               (((z)/((1+z)*(1-2*z)))*(dpy(j)*dpx(i))*Area)+(1/(2*(1+z))*(dpx(j)*dpy(i))*Area);
                               ((z)/((1+z)*(1-2*z))*(dpx(j)*dpy(i))*Area)+(1/(2*(1+z))*(dpy(j)*dpx(i))*Area)... 
                               (((1-z)/((1+z)*(1-2*z)))*(dpy(j)*dpy(i))*Area)+(1/(2*(1+z))*(dpx(j)*dpx(i))*Area)];

                        Kg(iix,jjx) = Kg(iix,jjx) + Kd(1,1);
                        Kg(iix,jjy) = Kg(iix,jjy) + Kd(1,2);

                        Kg(iiy,jjx) = Kg(iiy,jjx) + Kd(2,1);
                        Kg(iiy,jjy) = Kg(iiy,jjy) + Kd(2,2);
                    end
                end
            else
                continue
            end

        end

        B = -Kg*Um;
        dUdE = Ag\B;
        J(:,materials) = dUdE;
        
    end

    for clearing = 1:length(nodes)
        if rem(clearing,5) ~= 0

            Um(clearing*2 - 1,:) = 0;
            Um(clearing*2,:) = 0;

            J(clearing*2 - 1,:) = 0;
            J(clearing*2,:) = 0;
        else
            continue
        end
    end

    GE = J'*(Ud - Um);
    LHS = (J'*J);
    Ek = LHS\GE;
    E = Ek + E;
    
    E_matrix(iter,:) = E';

    E3E1(iter) = E(3)/E(1);
    E2E1(iter) = E(2)/E(1);

    L2(iter) = norm(Um - Ud);

end

plot(E_matrix(:,1))
hold on
plot(E_matrix(:,2))
plot(E_matrix(:,3))
ylabel('Stiffness (Pascals)')
xlabel('Iterations')
grid on
legend('Material 1', 'Material 2', 'Material 3')
hold off
title('Material Properties of 3 Tissue types in a Breast')


plot(E3E1)
hold on
plot(E2E1)
grid on
xlabel('Ratio')
ylabel('Iterations')
legend('Material 3 / Material 1', 'Material 2 / Material 1')
hold off
title('Ratios of Material 3 to 1 and Material 2 to 1 over Iterations')


plot(L2)
grid on
xlabel('L2 norm difference')
ylabel('Iterations')
hold off
title('L2 Norm Differences between Data and Model')


