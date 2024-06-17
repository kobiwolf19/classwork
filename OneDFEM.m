% numerical
pErrorVector = zeros(6,1);
hVector = zeros(6,1);
NODES = [6 11 21 41 81 161];

for k = 1:numel(NODES)

    a = 0.5;
    r = -10;
    N = 6; %%NODES(k); % change for number of nodes in between 1cm and 2cm
    h = 1/(N-1);

    for i = 0:N-2
        xEL(i+1) = 1 + h/2 + h*i;
    end

    incList = [1:N - 1]';
    incList(:,2) = [2:N]';

    Bg = zeros(N,1);
    Ag = zeros(N);

    for L = 1:N-1
        Ae = [(a*xEL(L).^2)/h,-(a*(xEL(L)).^2)/h;-(a*xEL(L).^2)/h,(a*(xEL(L)).^2)/h];
        Be = [-10*h/2;-10*h/2];
        for i = 1:2
            ii = incList(L,i);
            Bg(ii) = Bg(ii) + Be(i);
            for j = 1:2
                jj = incList(L,j);
                Ag(ii,jj) = Ag(ii,jj) + Ae(i,j);
            end
        end
    end

    Ag(1,:) = 0;
    Ag(1,1) = 1;
    Ag(N,:) = 0;
    Ag(N,N) = 1;

    Bg(1,1) = 10;
    Bg(N,1) = 5;

    p = Ag\Bg;

    % analytical

    p1 = 10;
    p2 = 5;

    AEq = -p1 + 2*p2 + ((2*r)/a)*log(2);
    BEq = 2*(p1 - p2 - (r/a)*log(2));

    x=1:h:2;

    for L = 1:N-1
        p0(L) = AEq + BEq/x(L) - (r/a)*log(x(L));
    end

    p0 = p0';
    p0(1) = 10;
    p0(N) = 5;

    for q = 1:N
        pError(q,1) = abs(p(q) - p0(q));
    end

    pErrorVector(k,1) = max(pError);
    hVector(k,1) = h;
end


plot(p,xEL)
xlabel('p values')
ylabel('x values')
title('p vs x values at N = 161')
grid on


plot(pErrorVector,hVector)
xlabel('Max P Error Values')
ylabel('Interval Values')
title('Maximum Error Values of p at Different Interval Sizes')
grid on

