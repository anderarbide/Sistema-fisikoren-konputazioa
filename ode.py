def leapfrog(lfdiffeq,r0,v0,t,h):
    """Ekuazio diferentzial arruntak ebazten ditu era bektorialean:
        dr[i]/dt=f[i] if lfdiffect(0,,,) edo dv[i]/dt=g[i] if lfdiffect(1,,,) bueltazen ditu"""
    hh=h/2.0
    r1=r0+hh*lfdiffeq(0,r0,v0,t)
    v1=v0+h*lfdiffeq(1,r1,v0,t+hh)
    r1=r1+hh*lfdiffeq(0,r0,v1,t+h)
    
    return r1,v1

def Euler(diffeq,y0,t,h):
    """y-ren balioa t-n hemanda y bueltatzen du t+h denboran"""
    dydt=diffeq(y0,t)

    return y0+h*dydt

def rk4(diffeq,y0,t,h):
    """Ekuazio diferentzialak ebazteko metodoa runge kuta 4 metodoa erabiliz:
            y0 t aldiunean emanda, y1 t+h aldiunean bueltatzen digu."""
    
    k1=h*diffeq(y0,t)
    k2=h*diffeq(y0+k1/2.0,t+h/2)
    k3=h*diffeq(y0+k2/2.0,t+h/2.0)
    k4=h*diffeq(y0+k3,t+h)
    
    return y0+k1/6.0+k2/3.0+k3/3.0+k4/6.0

def rk4b(diffeq,y0,t,h):
    """Ekuazio diferentzialak ebazteko metodoa runge kuta 4 era bektorialean erabiliz:
            y0[] t aldiunean emanda, y1[] t+h aldiunean bueltatzen digu."""
    
    {k1}=h*{diffeq({y0},t)}
    {k2}=h*{f({y0}+{k1}/2.0,t+h/2.0)}
    {k3}=h*{f({y0}+{k2}/2.0,t+h/2.0)}
    {k4}=h*{f({y0}+{k3},t+h)}

    return {y0}+{k1}/6.0+{k2}/3.0+{k3}/3.0+{k4}/6.0

def rk4b2(diffeq,y0,t,h):
    """Ekuazio diferentzialak ebazteko metodoa runge kuta 4 era bektorialean erabiliz:
            y0[] t aldiunean emanda, y1[] t+h aldiunean bueltatzen digu."""

    for i in y0:
        
        k1(i)=h*diffeq(y0(i),t)
        k2(i)=h*diffeq(y0(i)+k1(i)/2.0,t+h/2)
        k3(i)=h*diffeq(y0(i)+k2(i)/2.0,t+h/2.0)
        k4(i)=h*diffeq(y0(i)+k3(i),t+h)
        y0(i)=y0(i)+k1(i)/6.0+k2(i)/3.0+k3(i)/3.0+k4(i)/6.0

    return y0
