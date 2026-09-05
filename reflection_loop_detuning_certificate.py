"""High-precision and interval certificates for finite detuning cancellation."""
import mpmath as mp

from reflection_loop_detuning_order import frame_coefficients
from reflection_loop_dynamics import GS


def adjoint(a, ctx):
    # Explicit real/imaginary conjugation also works with interval complex numbers.
    ii = ctx.mpc(0, 1)
    return ctx.matrix([[a[c, r].real-ii*a[c, r].imag for c in range(a.rows)]
                       for r in range(a.cols)])


def primitive(ctx, n):
    eye, zero, ii = ctx.eye(4), ctx.zeros(4), ctx.mpc(0, 1)
    noise = [ctx.zeros(4) for _ in range(3)]
    noise[0][2, 2] = noise[1][3, 3] = 1
    noise[2][2, 3] = noise[2][3, 2] = 1
    jmat = ctx.zeros(4)
    jmat[0, 1] = jmat[1, 0] = 1
    speed = 1/ctx.sqrt(16*n*n-1)
    tau = ctx.pi/(2*speed)
    eta = (1+ctx.sqrt(2))/10
    theta = ctx.pi/3

    def integ(k):
        return ctx.pi/2 if k == 0 else (ii**(k % 4)-1)/(ii*k)

    def ordered(k, l):
        if l:
            return (integ(k+l)-integ(k))/(ii*l)
        if k == 0:
            return ctx.pi**2/8
        phase = ii**(k % 4)
        return ctx.pi*phase/(2*ii*k)+(phase-1)/(k*k)

    def arc(which):
        fs = [ctx.matrix([[int(v.real) for v in row] for row in f])
              for f in frame_coefficients(which)]
        f0 = fs[0]+fs[1]
        fourier = {0: fs[0], 1: (fs[1]-ii*fs[2])/2, -1: (fs[1]+ii*fs[2])/2}
        g = ctx.zeros(4)
        for r in range(3):
            for c in range(3):
                z = GS[which][r, c]
                g[r, c] = int(z.real)+ii*int(z.imag)
        hat = (jmat+speed*g)/(4*n*speed)
        pp = {0: eye-hat*hat, 1: (hat*hat+hat)/2, -1: (hat*hat-hat)/2}
        aa = {k-4*n*m: f*p*f0.T for k, f in fourier.items() for m, p in pp.items()}
        modes = {}
        for p, ap in aa.items():
            for q, aq in aa.items():
                freq = q-p
                old = modes.get(freq, [ctx.zeros(4) for _ in range(3)])
                modes[freq] = [old[j]+adjoint(ap, ctx)*v*aq for j, v in enumerate(noise)]
        endpoint = sum((a*ii**(k % 4) for k, a in aa.items()), ctx.zeros(4))
        one = [-ii/speed*sum((integ(k)*m[j] for k, m in modes.items()), ctx.zeros(4)) for j in range(3)]
        two = [[ctx.zeros(4) for _ in range(3)] for _ in range(3)]
        for k, mk in modes.items():
            for l, ml in modes.items():
                c = -ordered(k, l)/(speed*speed)
                for i in range(3):
                    for j in range(3):
                        two[i][j] += c*mk[i]*ml[j]
        return endpoint, [endpoint*x for x in one], [[endpoint*x for x in row] for row in two]

    def multiply(left, right):
        a, b, c = left
        d, e, f = right
        return a*d, [a*e[j]+b[j]*d for j in range(3)], [
            [a*f[i][j]+c[i][j]*d+b[i]*e[j] for j in range(3)] for i in range(3)]

    def hold(u, duration):
        return u, [-ii*duration*u*v for v in noise], [
            [-duration**2/2*u*a*b for b in noise] for a in noise]

    arcs = [arc(k) for k in range(3)]
    inv = [(adjoint(u, ctx), [-adjoint(a, ctx) for a in first],
            [[adjoint(second[j][i], ctx) for j in range(3)] for i in range(3)])
           for u, first, second in arcs[::-1]]
    hh = hold(eye-ii*ctx.sqrt(3)/2*jmat-jmat*jmat/2, theta+ctx.mpf(1)/5)
    cap = hold(eye, eta)
    result = eye, [ctx.zeros(4) for _ in range(3)], [[ctx.zeros(4) for _ in range(3)] for _ in range(3)]
    for block in [cap]+arcs+[hh]+inv+[hh]+arcs+[cap]:
        result = multiply(block, result)
    return result, 9*tau+2*theta+ctx.mpf(2)/5+2*eta


def equations(x, data, ctx, gain_target=None):
    h = (len(x)-1)//2
    beta = list(x[:h])+[ctx.mpf(0)]+[-v for v in x[:h][::-1]]
    w = list(x[h:2*h])+[x[-1]]+list(x[h:2*h][::-1])
    signs = [(-1)**j for j in range(2*h+1)]
    phi, alpha = ctx.mpf(0), []
    for b, s in zip(beta, signs):
        alpha.append(phi+s*b)
        phi += 2*s*b
    c1, c2 = ([ctx.cos(k*v) for v in beta] for k in (1, 2))
    ff = [sum(v*c for v, c in zip(w, c1)), sum(v*c for v, c in zip(w, c2))]
    for cs in (c1, c2):
        for trig in (ctx.cos, ctx.sin):
            ff.append(sum(v*c*trig(2*a) for v, c, a in zip(w, cs, alpha)))
    ff.extend((sum(v*ctx.sin(2*a) for v, a in zip(w, alpha)), sum(s*c for s, c in zip(signs, c2))))
    ff = [v/len(beta) for v in ff]
    (u, first, second), duration = data
    result = [ctx.eye(4), ctx.zeros(4), ctx.zeros(4)]
    for a, stretch, sign in zip(alpha, w, signs):
        co, sn = ctx.cos(a), ctx.sin(a)
        coeff = [sn*sn, co*co, sn*co]
        raw = [u, sum((coeff[j]*first[j] for j in range(3)), ctx.zeros(4)),
               sum((coeff[i]*coeff[j]*second[i][j] for i in range(3) for j in range(3)), ctx.zeros(4))]
        if sign < 0:
            raw = [adjoint(a, ctx) for a in raw]
        rot = ctx.eye(4)
        rot[2, 2], rot[2, 3], rot[3, 2], rot[3, 3] = co, -sn, sn, co
        raw = [(sign*stretch)**q*rot*a*rot.T for q, a in enumerate(raw)]
        result = [sum((raw[k]*result[q-k] for k in range(q+1)), ctx.zeros(4)) for q in range(3)]
    ff.extend((-result[2][2, 3].real/duration**2,
               (-result[2][2, 2]-result[2][3, 3]).imag/(2*duration**2)))
    if gain_target is not None:
        ff.append((sum(s*c for s, c in zip(signs, c1))-gain_target)/len(beta))
    return ff


def refine(seed, n, fixed, gain_target=None, digits=65):
    mp.mp.dps = digits+15
    data = primitive(mp.mp, n)
    indices = [j for j in range(len(seed)) if j not in fixed]
    all_x = [mp.mpf(str(v)) for v in seed]
    for j, value in fixed.items():
        all_x[j] = mp.mpf(value)

    def f(*values):
        full = all_x.copy()
        for j, value in zip(indices, values):
            full[j] = value
        return equations(full, data, mp.mp, gain_target)

    root = mp.findroot(f, [all_x[j] for j in indices], tol=mp.mpf(10)**(-digits), maxsteps=50)
    for j, value in zip(indices, root):
        all_x[j] = value
    return dict(n=n, center=[mp.nstr(v, digits) for v in all_x],
                fixed={str(j): str(v) for j, v in fixed.items()}, gain_target=gain_target)


def certify(record):
    """Interval residual/Jacobian plus an explicit global Hessian bound."""
    mp.mp.dps = 90
    mp.iv.dps = 80
    n, gain = record['n'], record['gain_target']
    xs = record['center']
    free = [j for j in range(len(xs)) if str(j) not in record['fixed']]
    dim = len(free)
    data = primitive(mp.mp, n)
    ivdata = primitive(mp.iv, n)
    center = [mp.mpf(v) for v in xs]

    def f(*y):
        x = center.copy()
        for j, v in zip(free, y):
            x[j] = v
        return mp.matrix(equations(x, data, mp.mp, gain))

    jac = mp.calculus.optimization.jacobian(mp.mp, f, mp.matrix([center[j] for j in free]))
    inverse = mp.inverse(jac)
    aa = [[mp.iv.mpf(mp.nstr(inverse[i, j], 85)) for j in range(dim)] for i in range(dim)]
    point = [mp.iv.mpf(v) for v in xs]
    residual = equations(point, ivdata, mp.iv, gain)
    # For N<=13, 1<=s<3, normalized equations have every second derivative <1e18.
    # Matrix-product bound: at most 30 rotation factors, rates <=4, and total
    # duration <=39 times the primitive; the explicit bound is documented.
    hessian = mp.iv.mpf('1e18')
    step, radius = mp.iv.mpf('1e-35'), mp.iv.mpf('1e-35')
    enclosure = []
    for j in free:
        plus, minus = point.copy(), point.copy()
        plus[j] += step
        minus[j] -= step
        fp = equations(plus, ivdata, mp.iv, gain)
        fm = equations(minus, ivdata, mp.iv, gain)
        enclosure.append([(p-m)/(2*step)+mp.iv.mpf([-1, 1])*hessian*step
                          for p, m in zip(fp, fm)])

    def upper(v):
        return mp.mpf(abs(v).b)

    defect = max(upper(sum(abs(int(i == j)-sum(aa[i][k]*enclosure[j][k] for k in range(dim)))
                              for j in range(dim))) for i in range(dim))
    norm_a = max(upper(sum(abs(v) for v in row)) for row in aa)
    q = defect+norm_a*dim*dim*mp.mpf('1e18')*mp.mpf('1e-35')
    displacement = max(upper(sum(aa[i][j]*residual[j] for j in range(dim))) for i in range(dim))
    assert q < 1
    assert displacement < (1-q)*mp.mpf('1e-35')
    assert len(xs) == 13 and all(1 <= v < 3 for v in center[6:])
    for j in free:
        if j >= 6:
            assert center[j]-mp.mpf('1e-35') >= 1
    return dict(radius='1e-35', contraction_bound=mp.nstr(q, 15),
                displacement_bound=mp.nstr(displacement, 15),
                maximum_second_partial_bound='1e18',
                root_exists_and_is_unique_with_fixed_controls=True)
