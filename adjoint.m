nodes = load ('box.nod');
elements = load('box.elm');
data = load('box.data');
boundary_elm = load('box.bel');
property = load('box.prop');
boundary_cond = load('box.bcs');

Ud = zeros(54*2,1);

for q = 1:length(data)
    node = data(q,1);
    Ud(node*2-1,1) = data(q,2);
    Ud(node*2,1) = data(q,3);
end

incTri = elements(:,2:4);
x = nodes(:,2);
y = nodes(:,3);

J = zeros(2*54,54);

E = zeros(length(nodes),1);
E(1:end) = 18000;
z = 0.45;

Ag = zeros(2*length(nodes),2*length(nodes));
Bg = zeros(length(nodes)*2,1);

for L = 1:length(elements)

    [dpx, dpy, Area] = basisfunc3(L,x,y,incTri);

    for i = 1:3
        iix = 2*incTri(L,i) - 1;
        iiy = iix + 1;

        for j = 1:3
            jjx = 2*incTri(L,j) - 1;
            jjy = jjx + 1;

            for kx = 1:3
                NOI = elements(L,kx + 1);
                K = E(NOI);
                Ke = [((K*(1-z)/((1+z)*(1-2*z)))*(dpx(i)*dpx(j))*Area/3)+(K/(2*(1+z))*(dpy(i)*dpy(j))*Area/3)...
                    (((K*z)/((1+z)*(1-2*z)))*(dpy(j)*dpx(i))*Area/3)+(K/(2*(1+z))*(dpx(j)*dpy(i))*Area/3);
                    ((K*z)/((1+z)*(1-2*z))*(dpx(j)*dpy(i))*Area/3)+(K/(2*(1+z))*(dpy(j)*dpx(i))*Area/3)...
                    ((K*(1-z)/((1+z)*(1-2*z)))*(dpy(j)*dpy(i))*Area/3)+(K/(2*(1+z))*(dpx(j)*dpx(i))*Area/3)];

                Ag(iix,jjx) = Ag(iix,jjx) + Ke(1,1);
                Ag(iix,jjy) = Ag(iix,jjy) + Ke(1,2);

                Ag(iiy,jjx) = Ag(iiy,jjx) + Ke(2,1);
                Ag(iiy,jjy) = Ag(iiy,jjy) + Ke(2,2);
            end
        end
    end
end

for N = 1:length(boundary_cond) % BCs
    if boundary_cond(N,3) == 1
        x_bc = boundary_cond(N,5);
        y_bc = boundary_cond(N,8);
        node = boundary_cond(N,2);
        Ag(node*2 - 1,:) = 0;
        Ag(node*2,:) = 0;
        Ag(node*2 - 1,node*2 - 1) = 1;
        Ag(node*2,node*2) = 1;
        Bg(node*2 - 1) = x_bc;
        Bg(node*2) = y_bc;
    else
        continue
    end
end

Um = Ag\Bg;

for materials = 1:length(nodes)

    Kg = zeros(2*length(nodes),2*length(nodes));

    for L = 1:length(elements)

        [dpx, dpy, Area] = basisfunc3(L,x,y,incTri);

        for kx = 1:3

            NOI = elements(L,kx + 1);

            if NOI == materials
                for i = 1:3
                    iix = 2*incTri(L,i) - 1;
                    iiy = iix + 1;

                    for j = 1:3
                        jjx = 2*incTri(L,j) - 1;
                        jjy = jjx + 1;

                        Kd = [(((1-z)/((1+z)*(1-2*z)))*(dpx(i)*dpx(j))*Area/3)+(1/(2*(1+z))*(dpy(i)*dpy(j))*Area/3)...
                            (((z)/((1+z)*(1-2*z)))*(dpy(j)*dpx(i))*Area/3)+(1/(2*(1+z))*(dpx(j)*dpy(i))*Area/3);
                            ((z)/((1+z)*(1-2*z))*(dpx(j)*dpy(i))*Area/3)+(1/(2*(1+z))*(dpy(j)*dpx(i))*Area/3)...
                            (((1-z)/((1+z)*(1-2*z)))*(dpy(j)*dpy(i))*Area/3)+(1/(2*(1+z))*(dpx(j)*dpx(i))*Area/3)];

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
    end

    B = Kg*Um;
    dUdE = Ag\B;
    J(:,materials) = dUdE;

    adj = Ag'\(-2*(Ud-Um));
    gradient(materials) = adj'*(-Kg*Um);

end

dOdE = -2*J'*(Ud-Um); %newton 

gradient = gradient'; %adjoint

plot(gradient)
hold on 
plot(-dOdE, 'o')
grid on 
legend('newton method','adjoint method')
xlabel('node number')
ylabel('dOmega/dStiffness')


