%% part b
Be = [0; 0; 0];
nodes = load('vessel.nod');
elements = load('vessel.elm');
bound = load('vessel.bel');

x = nodes(:,2);
y=nodes(:,3);
inc = elements(:,2:4);

Bg = zeros(28,1);
Ag = zeros(28,28);


for L = 1:38
    [dpx, dpy, Area] = basisfunc3(L,x,y,inc);

    Ae = [(dpx(1)^2)*Area + (dpy(1)^2)*Area (dpx(1)*dpx(2))*Area + (dpy(1)*dpy(2))*Area (dpx(1)*dpx(3))*Area + (dpy(1)*dpy(3))*Area;
        (dpx(1)*dpx(2))*Area + (dpy(1)*dpy(2))*Area (dpx(2)^2)*Area + (dpy(2)^2)*Area (dpx(2)*dpx(3))*Area + (dpy(2)*dpy(3))*Area;
        (dpx(1)*dpx(3))*Area + (dpy(1)*dpy(3))*Area (dpx(2)*dpx(3))*Area + (dpy(2)*dpy(3))*Area (dpx(3)^2)*Area + (dpy(3)^2)*Area];
    
    for i = 1:3
        ii = inc(L,i);
        Bg(ii) = Bg(ii) + Be(i);

        for j = 1:3
            jj=inc(L,j);
            Ag(ii,jj) = Ag(ii,jj) + Ae(i,j);
        end
    end
end

for P = 1:28
    if P == 1||P ==2||P ==3||P ==4||P ==5||P ==10||P ==15||P ==20||P ==24||P ==27||P ==28||P ==6||P ==11||P ==16
        Ag(P,:) = 0;
        Ag(P,P) = 1;
        if P == 5||P ==10||P ==15||P ==20||P ==24||P ==27||P ==28
            Bg(P,1) = 1;
        elseif P == 4
            Bg(P,1) = 0.75;
        elseif P == 3
            Bg(P,1) = 0.5;
        elseif P == 2
            Bg(P,1) = 0.25;
        end
    else
        continue
    end
end

psi = Ag\Bg;

dPSIx = zeros(38,1);
dPSIy = zeros(38,1);

for N = 1:38
    [dpx, dpy, Area] = basisfunc3(N,x,y,inc);
    dPSIx(N) = psi(inc(N,1))*dpx(1) + psi(inc(N,2))*dpx(2) + psi(inc(N,3))*dpx(3);
    dPSIy(N) = psi(inc(N,1))*dpy(1) + psi(inc(N,2))*dpy(2) + psi(inc(N,3))*dpy(3);
end

%% part c
Be = [0; 0; 0];
nodes = load('vessel.nod');
elements = load('vessel.elm');
bound = load('vessel.bel');

x = nodes(:,2);
y=nodes(:,3);
inc = elements(:,2:4);

Bg = zeros(28,1);
Ag = zeros(28,28);


for L = 1:38
    [dpx, dpy, Area] = basisfunc3(L,x,y,inc);

    Ae = [(dpx(1)^2)*Area + (dpy(1)^2)*Area (dpx(1)*dpx(2))*Area + (dpy(1)*dpy(2))*Area (dpx(1)*dpx(3))*Area + (dpy(1)*dpy(3))*Area;
        (dpx(1)*dpx(2))*Area + (dpy(1)*dpy(2))*Area (dpx(2)^2)*Area + (dpy(2)^2)*Area (dpx(2)*dpx(3))*Area + (dpy(2)*dpy(3))*Area;
        (dpx(1)*dpx(3))*Area + (dpy(1)*dpy(3))*Area (dpx(2)*dpx(3))*Area + (dpy(2)*dpy(3))*Area (dpx(3)^2)*Area + (dpy(3)^2)*Area];
    
    for i = 1:3
        ii = inc(L,i);
        Bg(ii) = Bg(ii) + Be(i);

        for j = 1:3
            jj=inc(L,j);
            Ag(ii,jj) = Ag(ii,jj) + Ae(i,j);
        end
    end
end

for P = 1:28
    if P == 1||P ==2||P ==3||P ==4||P ==5||P ==10||P ==15||P ==20||P ==24||P ==27||P ==28||P ==6||P ==11||P ==16 || P==25
        Ag(P,:) = 0;
        Ag(P,P) = 1;
        if P == 5||P ==10||P ==15||P ==20||P ==24||P ==27||P ==28 || P==25
            Bg(P,1) = 1;
        elseif P == 4
            Bg(P,1) = 0.75;
        elseif P == 3
            Bg(P,1) = 0.5;
        elseif P == 2
            Bg(P,1) = 0.25;
        end
    else
        continue
    end
end

psi25 = Ag\Bg;

dPSIx = zeros(38,1);
dPSIy = zeros(38,1);

for N = 1:38
    [dpx, dpy, Area] = basisfunc3(N,x,y,inc);
    dPSIx(N) = psi25(inc(N,1))*dpx(1) + psi25(inc(N,2))*dpx(2) + psi25(inc(N,3))*dpx(3);
    dPSIy(N) = psi25(inc(N,1))*dpy(1) + psi25(inc(N,2))*dpy(2) + psi25(inc(N,3))*dpy(3);
end

%%

[p,e,t] = pet;






