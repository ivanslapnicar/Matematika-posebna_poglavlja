### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 8a4b7c05-081c-4913-8278-838f865bd2eb
md"""
# Greenova funkcija

---

Promotrimo regularni Sturm-Liouvilleov problem:

\begin{align*}
Au& \equiv -(p(x)\,u')'+q(x)\,u = f, \quad a\leq x\leq b,\\
B_1u(a)&=\alpha_1u(a)+\alpha_2u'(a)=0,\\
B_2u(b)&=\beta_1 u(b)+\beta_2u'(b)=0,
\end{align*}

gdje je 

$$
\quad p,p',q,f,u,u',u''\in C[a,b],\quad p(x)>0.
$$

Zadatak je naći rješenje $u(x)$.

Operator $A$ je linearan, a kompletni problem (zajedno s rubnim uvjetima) možemo prikazati u obliku $Lu=f$, pri čemu je $L$ linearan operator __koji u sebi sadrži i rubne uvjete__.

Ako se radi o matricama, onda sustav jednadžbi $Lu=f$, pri čemu je matrica $L$ regularna ($\lambda=0$ nije svojstvena vrijednost matrice $L$), možemo riješiti invertiranjem:

$$
x=L^{-1}f. \tag{1}
$$

__Pitanje.__ Možemo li slično napraviti i za regularni SLP, odnosno možemo li naći inverzni operator 
$L^{-1}$ tako da vrijedi (1)? 

__Odgovor je potvrdan!__ Inverzni operator je integralni operator

$$
u(x)\equiv (L^{-1} f)(x)=\int_a^b g(x,\xi) f(\xi) \, d\xi, \tag{2}
$$

pri čemu je __Greenova funkcija__ $g(x,\xi)$ njegova jezgra.
"""

# ╔═╡ f8ebd0bc-794b-4018-b023-f72f890aba08
md"""
__Teorem.__ Neka je operator $L$ regularan ($\lambda=0$ nije svojstvena vrijednost). Tada $L^{-1}$ postoji i definiran je s (2), gdje je 

$$
g(x,\xi)=\begin{cases} -\displaystyle\frac{u_1(x)u_2(\xi)}{p(\xi)W(\xi)}, \quad x\leq\xi, \\
-\displaystyle\frac{u_1(\xi)u_2(x)}{p(\xi)W(\xi)},\quad x>\xi.
\end{cases}
$$

Ovdje su $u_1$ i $u_2$ linearno nezavisna rješenja homogene diferencijalne jednadžbe $Au=0$ uz rubne uvjete $B_1u_1(a)=0$ i $B_2u_2(b)=0$, a 

$$
W(x)=\begin{vmatrix} u_1 & u_2 \\ u_1' & u_2' \end{vmatrix}=u_1u_2'-u_1'u_2
$$

je Wronskijan.
"""

# ╔═╡ 61a891d8-6f5c-40a1-a6e1-335525f60c40
md"""
## Svojstva Greenove funkcije

__Heavisideova step funkcija__ $H(x)$ je

$$
H(x)=\begin{cases}1, \quad x > 0,\\ \frac{1}{2},\quad x=0,\\ 0,\quad x<0. \end{cases}
$$

__Greenovu funkcija__ je

$$
g(x,\xi)=-\displaystyle\frac{1}{p(\xi)W(\xi)}[u_1(x)u_2(\xi) H(\xi-x)
+u_1(\xi)u_2(x)H(x-\xi)]. \tag{3}
$$

Greenova funkcija je __neprekidna__: $u_1$, $u_2$ i njihove derivacije su neprekidne
funkcije po pretpostavci. Funkcija $p$ je neprekidna i različita od nule, a Wronskijan je neprekidan i različit od nule jer su funkcije $u_1$ i $u_2$ linearno nezavisne. Posebno, Greenova funkcija je očito neprekidna u točki $x=\xi$ (limesi s lijeva i zdesna su jednaki!).

Derivacija $g_x(x,\xi)$ ima skok u točki $x=\xi$: vrijedi

\begin{align*}
g_x(\xi^+,\xi)&=-\displaystyle\frac{1}{p(\xi)W(\xi)}[u_1'(\xi)u_2(\xi) H(0^-)
+u_1(\xi)u_2(\xi)H'(0^-)(-1)\\
&+u_1(\xi)u_2'(\xi)H(0^+)+u_1(\xi)u_2(\xi)H'(0^+)]\\
&=-\displaystyle\frac{u_1(\xi)u_2'(\xi)}{p(\xi)W(\xi)}.
\end{align*}

Ovdje smo koristili činjenice da je 

$$
H(0^+)=1,\quad H(0^-)=0,\quad H'(0^+)=0,\quad H'(0^-)=0.
$$

Slično se pokaže da je 

$$
g_x(\xi^-,\xi)=-\displaystyle\frac{u_1'(\xi)u_2(\xi)}{p(\xi)W(\xi)}
$$

pa je 

\begin{align*}
g_x(\xi^+,\xi)-g_x(\xi^-,\xi)&=-\frac{1}{p(\xi)W(\xi)}[u_1(\xi)u_2'(\xi)
-u_1'(\xi)u_2(\xi)]\\ &=-\frac{1}{p(\xi)W(\xi)}W(\xi)=-\frac{1}{p(\xi)}\neq 0.
\end{align*}
"""

# ╔═╡ d172800a-341c-4224-b933-51cda3121014
md"""
## Primjer 

Riješimo problem rubnih vrijednosti

$$
-u''(x)=f(x),\quad 0<x<1,\quad u(0)=0,\quad u(1)=0,
$$

pomoću Greenove funkcije.

Vrijedi $p(x)=1$. Homogeni problem $Au=0$ glasi $u''(x)=0$ pa je $u=ax+b$.
Tražimo dva linearno nezavisna rješenja tako da je $u_1(0)=0$ i $u_2(1)=0$.

Možemo uzeti $u_1=x$. Uvjet $u_2(1)=0$ povlači $a\cdot 1+b=0$ pa je $a=-b$, odnosno
$u_2=a(x-1)$. Konstanta $a$ je proizvoljna pa uzmimo $u_2=x-1$. Vrijedi

$$
W(x)=\begin{vmatrix} x & x-1 \\ 1 & 1 \end{vmatrix}=x-(x-1)=1\neq 0,
$$

pa su funkcije $u_1$ i $u_2$ linearno nezavisne.

Prema formuli (3) je (nacrtajte funkciju $g(x,\xi)$ kao funkciju od $x$ za neki $\xi$)

$$
g(x,\xi)=-x(\xi-1)H(\xi-x)-\xi(x-1)H(x-\xi)
$$

pa je, prema (2), 

\begin{align*}
u(x)&=\int_0^1 x(1-\xi)H(\xi-x)\,f(\xi)\, d\xi
+\xi(1-x)H(x-\xi)\,f(\xi)\,d\xi\\
&=\int_0^1 [x(1-\xi)(1-H(x-\xi))+\xi(1-x)H(x-\xi)]\,f(\xi)\, d\xi\\
&=x\int_0^1 (1-\xi)\,f(\xi)\, d\xi+\int_0^1 [\xi(1-x)-x(1-\xi)]H(x-\xi)\,f(\xi)\,d\xi\\
&=x\int_0^1 (1-\xi)\,f(\xi)\, d\xi+\int_0^x (\xi-x)\,f(\xi)\,d\xi.
\end{align*}
"""

# ╔═╡ 649f28a3-f4f2-48d1-a1e4-b79b378c522c
using SymPy, Plots, QuadGK

# ╔═╡ 159a90f8-2bf4-40ad-9c31-de1692cea6c4
x,y=symbols("x,y",real=true)

# ╔═╡ c7483547-316d-4c82-affc-2aaf5a203cbf
# ξ=0.5
g(x)=x*(1-ξ)*Heaviside(ξ-x)+ξ*(1-x)*Heaviside(x-ξ)

# ╔═╡ 05c5f5ca-9931-4419-ae67-aa1b8236a599
ξ=0.5
plot(g(x),0,1)

# ╔═╡ 0db3d836-fd58-4597-a183-7043059d20dd
ξ=0.8
plot(g(x),0,1)

# ╔═╡ 3f1178ff-83a4-42d1-8c6a-fb39dd4de4ff
md"""
Fizikalna interpretacije je sljedeća:

> Green-ova funkcija je odgovor linearnog operatora $Au$ na impuls u točki $x=\xi$ (neovisno o funkciji $f$) koji zadovolja i rubne uvjete, pa integral (2) daje rješenje problema rubnih vrijednosti za zadanu funkciju $f$

Provjerimo zadovoljava li 

$$u(x)=x\int_0^1 (1-\xi)\,f(\xi)\, d\xi+\int_0^x (\xi-x)\,f(\xi)\,d\xi$$ 

diferencijalnu jednadžbu i rubne uvjete: 

$$
-u''(x)=f(x),\quad 0<x<1,\quad u(0)=0,\quad u(1)=0.
$$

Leibnitzovo pravilo daje

\begin{align*}
u'(x)&=\int_0^1 (1-\xi)\,f(\xi)\,d\xi+\int_0^x (-1)\,f(\xi)\,d\xi +(x-x)\,f(x)\cdot 1-
(0-x)\,f(0)\cdot 0 \\
&=\int_0^1 (1-\xi)\,f(\xi)\,d\xi-\int_0^x f(\xi)\,d\xi\\
u''(x)&=0-\int_0^x 0\cdot d\xi -f(x)\cdot 1 + 0=-f(x).
\end{align*}

U rubovima vrijedi

\begin{align*}
u(0)&=0\cdot \int_0^1 (1-\xi)\,f(\xi)\, d\xi+\int_0^0 (\xi-0)\,f(\xi)\,d\xi=0\\
u(1)&=1\cdot \int_0^1 (1-\xi)\,f(\xi)\, d\xi+\int_0^1 (\xi-1)\,f(\xi)\,d\xi=0.
\end{align*}
"""

# ╔═╡ be021115-1978-49fa-81a4-434c383883b0
md"""
Na primjer, za $f(x)=\sin(x)$ lako provjerimo da rješenje problema glasi

$$
u(x)=\sin x -x\sin 1.
$$
"""

# ╔═╡ f122c3c3-1e29-47e7-8fc6-e7e013920962
u(x)=sin(x)-x*sin(1)

# ╔═╡ 88107568-140d-4485-a680-ef1e74241f78
plot(u,0,1)

# ╔═╡ d923f52a-b7b1-428e-9e81-40c45db0deba
diff(u(x),x,x)

# ╔═╡ 2b88038f-09e7-4722-9360-170c06fa1eaf
u(0), u(1)

# ╔═╡ 80613a67-6843-44ed-be63-51c7cd62693e
# Rješenje pomoću Greenove funkcije
f(x)=sin(x)
ug(x)=x*integrate((1-y)*f(y),(y,0,1))+integrate((y-x)*f(y),(y,0,x))

# ╔═╡ 95f14d64-1c0a-4f67-a6fa-82616b57fc01
ug(0.5)

# ╔═╡ e30fd6a4-7153-4689-b343-abd8898c5410
# Ovaj simbolički račun traje duže
plot(ug,0,1)

# ╔═╡ 119ba9d5-ced7-4ebd-8cf3-c127ef2f0aec
md"""
## Delta funkcija

__[Diracova Delta funkcija](https://en.wikipedia.org/wiki/Dirac_delta_function)__ je funkcija na skupu realnih brojeva koja je nula svugdje osim u ishodištu, gdje je beskonačna (idealizirani impuls),

$$\displaystyle \delta (x)=\begin{cases}+\infty,\quad  x=0\\0,\qquad  x\neq 0,\end{cases}$$

te zadovoljava ograničenje

$$\displaystyle \int\limits _{-\infty }^{\infty }\delta (x)\,dx=1.$$

Diracova delta nije funkcija u klasičnom smislu jer ni jedna funckija definirana na skupu realnih brojeva ne može imati ova svojstva. Međutim, iz ove definicije 
slijedi __svojstvo sampliranja__:

$$f(x)=\displaystyle \int\limits _{-\infty }^{\infty }f(t)\delta (t-x)\,dt$$

za svaku realnu funkciju $f$.

__Heavisideova step funkcija__ je definirana s

$$H(x)=\begin{cases}0,\quad  x<0,\\ \displaystyle\frac{1}{2},\quad  x=0,\\
1,\quad  x>0.\end{cases}$$

Vrijedi 

$$\frac{d}{dx} H(x)=\delta(x).$$

Zaista, za $x\neq 0$ vrijedi

$$H'(x)=\lim_{\Delta x\to 0} \frac{H(x+\Delta x)-H(x)}{\Delta x}=
\lim_{\Delta x\to 0} \frac{0}{\Delta x}=0,$$ 

dok za $x=0$ vrijedi

$$
H'(0_-)=\lim_{\Delta x\to 0_-} \frac{H(\Delta x)-H(0)}{\Delta x}=
\lim_{\Delta x\to 0_-} \frac{0-\frac{1}{2}}{\Delta x}=+\infty
$$

i

$$
H'(0_+)=\lim_{\Delta x\to 0_+} \frac{H(\Delta x)-H(0)}{\Delta x}=
\lim_{\Delta x\to 0_+} \frac{1-\frac{1}{2}}{\Delta x}=+\infty.
$$
"""

# ╔═╡ 2e4a1f96-005f-441d-beab-5b6a927c1e2e
md"""
__Primjer.__ Riješimo problem rubnih vrijednosti

$$
u''+u'=f(x),\quad 0<x<1,\quad u(0)=u_0,\quad u(1)=u_1.
$$

Rubni uvjeti nisu homogeni pa ne možemo koristiti teorem. Želimo rješenje u obliku

$$
u(x)=\int_0^1 g(x,\xi)\, f(\xi)\, d\xi + B \tag{1}
$$

gdje $B$ daje utjecaj rubnih uvjeta. Moramo odrediti problem rubnih vrijednosti koji funkcija $g$ zadovoljava.

Diracova delta funkcija daje

$$
u(x)=\int_0^1 \delta(\xi-x)\, u(\xi)\, d\xi. \tag{2}
$$

Uvrštavanjem DJ u (1) imamo

$$
u(x)=\int_0^1 g(x,\xi)\,\frac{d}{d\xi}\bigg(\frac{d}{d\xi}+1\bigg)\,u(\xi)\, d\xi + B
$$

Želimo postići da diferencijalni operatori djeluju na funkciju $g$. Primijenimo dva puta parcijalnu integraciju: prva parcijalna integracija daje

$$
u(x)=g(x,\xi)\,\bigg(\frac{d}{d\xi}+1\bigg)\,u(\xi) \, \Bigg|_0^1 -
\int_0^1 \frac{d}{d\xi}g(x,\xi)\cdot\bigg(\frac{d}{d\xi}+1\bigg)\, u(\xi)\, d\xi +B \\
$$

Parcijalna integracija preostalog integrala daje

$$
\int_0^1 \frac{\partial g}{\partial \xi}\frac{du}{d\xi}\, d\xi+
\int_0^1 u\frac{\partial g}{\partial \xi}\, d\xi
=\frac{\partial g}{\partial \xi} u\, \bigg|_0^1 -\int_0^1 u\frac{\partial^2 g}{\partial \xi^2}\, d\xi + \int_0^1 u\frac{\partial g}{\partial \xi}\, d\xi.
$$

Uvrštavanjem i korištenjem (2) imamo

\begin{align*}
u(x)&=\int_0^1 \bigg(\frac{\partial^2 g}{\partial \xi^2}-
\frac{\partial g}{\partial \xi}\bigg)\, u\, d\xi +
\bigg[g\bigg(\frac{du}{d\xi}+u\bigg)-\frac{\partial g}{\partial \xi} u\bigg]_0^1+B
\\ &=\int_0^1 \delta(\xi-x)\, u(\xi)\, d\xi.
\end{align*}

Jednakost će biti zadovoljena ako je 

\begin{align*}
\frac{\partial^2 g}{\partial \xi^2}-
\frac{\partial g}{\partial \xi}&=\delta(\xi-x) \tag{3} \\
B&=\bigg[\frac{\partial g}{\partial \xi} u-
g\bigg(\frac{du}{d\xi}+u\bigg)\bigg]_0^1 \tag{4}
\end{align*}
"""

# ╔═╡ aa4614d8-04ee-4de9-979b-eab95a84f1e3
md"""
Vrijednosti $u(0)$ i $u(1)$ su zadane, dok $u'(0)$ i $u'(1)$ nisu, pa stavimo

$$
g(x,0)=0,\quad g(x,1)=0, \tag{5}
$$

dok ćemo vrijednosti $g_\xi(x,0)$ i 
$g_\xi(x,1)$ izračunati.

Integriranje (3) po $\xi$ daje (derivacija Heavisideove funkcije je delta funkcija)

$$
\frac{\partial g}{\partial \xi}-g=K_1+H(\xi-x).
$$

Ovo je [linearna diferencijalna jednadžba prvog reda](http://lavica.fesb.unist.hr/matematika2/predavanja/node91.html) čije rješenje glasi

$$
g(x,\xi)=e^\xi \left(K_2 + \int e^{-\xi}(K_1+H(\xi-x))\, d\xi\right).
$$


Osnovni teorem integralnog računa, $\displaystyle\int f(x)\, dx=\int_0^x f(t)\, dt$, daje

\begin{align*}
g(x,\xi)&=e^\xi \left(K_2 + \int_0^\xi e^{-\zeta}(K_1+H(\zeta-x))\, d\zeta\right)\\
&=e^\xi \left(K_2+K_1\big(-e^{-\zeta}\bigg|_0^\xi \big)+
\int_0^\xi e^{-\zeta} H(\zeta-x)\, d\zeta \right)\\
&=K_2 e^\xi + K_1(e^\xi-1)+\frac{1}{2} e^\xi
\left( |e^{-x}-e^{-\xi}|+e^{-x}-e^{-\xi}\right). \tag{6}
\end{align*}

Zaista, za $x<\xi$ je

\begin{align*}
\int_0^\xi e^{-\zeta} H(\zeta-x)\, d\zeta&= \int_0^x e^{-\zeta} H(\zeta-x)\, d\zeta
+\int_x^\xi e^{-\zeta} H(\zeta-x)\, d\zeta\\ &=\int_0^x e^{-\zeta} 0\, d\zeta
+\int_x^\xi e^{-\zeta} 1\, d\zeta=-e^{-\zeta}\, \bigg|_x^\xi\\
&=-e^{-\xi}+e^{-x}=e^{-x}-e^{-\xi},
\end{align*}

a za $x>\xi$ je integral jednak nuli. S druge strane je 

$$
\frac{1}{2}
\left( |e^{-x}-e^{-\xi}|+e^{-x}-e^{-\xi}\right)=\begin{cases} 
e^{-x}-e^{-\xi},\quad & x<\xi,\\
0,\quad & x>\xi. \end{cases}
$$

Iz (5) slijedi

\begin{align*}
0&=g(x,0)=K_2\cdot 1 + K_1(1-1)+\frac{1}{2} \cdot 1
\left( |e^{-x}-1|+e^{-x}-1\right)=K_2\\
0&=g(x,1)=K_1(e-1)+e(e^{-x}-e^{-1})=K_1(e-1)+e^{1-x}-1
\end{align*}

pa je 

$$
K_1=\frac{1-e^{1-x}}{e-1}.
$$
"""

# ╔═╡ 6030f0b9-677e-4fde-b9d1-d01541628c04
md"""
Uvrštavanje u (6) daje

$$
g(x,\xi)=\frac{1-e^{1-x}}{e-1}(e^\xi-1)+\frac{1}{2} 
e^\xi\left( |e^{-x}-e^{-\xi}|+e^{-x}-e^{-\xi}\right). \tag{7}
$$

Da bismo odredili $B$ iz (4) trebamo još izračunati $g_\xi(x,0)$ i 
$g_\xi(x,1)$. Vrijedi

\begin{align*}
g_\xi(x,\xi)=&\frac{1-e^{1-x}}{1-e}e^\xi+
\frac{1}{2} e^\xi\left( |e^{-x}-e^{-\xi}|+e^{-x}-e^{-\xi}\right)\\
&+ \frac{1}{2} e^\xi\left( \frac{e^{-x}-e^{-\xi}}{|e^{-x}-e^{-\xi}|}
e^{-\xi}+e^{-\xi}\right)
\end{align*}

pa je 

\begin{align*}
g_\xi(x,0)&=\frac{1-e^{1-x}}{e-1}+
\frac{1}{2}\left( |e^{-x}-1|+e^{-x}-1\right)+\frac{1}{2}\left( \frac{e^{-x}-1}{|e^{-x}-1|}+1\right)\\
&=\frac{1-e^{1-x}}{e-1}\\
g_\xi(x,1)&=e\frac{1-e^{1-x}}{e-1}+
e(e^{-x}-e^{-1}) +\frac{1}{2}e(e^{-1}+e^{-1})\\
&=e\frac{1-e^{1-x}}{e-1}+e^{1-x}.
\end{align*}

Uvrštavanjem u (4) konačno slijedi

$$
u(x)=\int_0^1 g(x,\xi)\, f(\xi)\, d\xi + u_1 
\bigg(e\frac{1-e^{1-x}}{e-1}+e^{1-x}\bigg)-u_0\frac{1-e^{1-x}}{e-1},
$$

gdje je $g(x,\xi)$ dano sa (7).
"""

# ╔═╡ b63428a6-cb0a-4b8f-8688-eb0a72fc6774
md"""
Nacrtajmo rješenje za $f(x)=\sin(2x)$. Rješenje ćemo izračunati pomoću numeričke integracije koristeći paket [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).
"""

# ╔═╡ 01fef681-189d-4963-89b7-f8d06b6ce588
x=Float64
y=Float64
f(x)=sin(2*x)
u₀=1
u₁=2

# ╔═╡ 321cfc8a-e5ca-476d-8393-f8d23eb04f64
using Base.MathConstants

# ╔═╡ f7409dfc-e679-411f-9fd9-79a3ff9fe835
g(x,ξ)=(1-exp(1-x))*(exp(ξ)-1)/(e-1) + 
    exp(ξ)*(abs(exp(-x)-exp(-ξ))+exp(-x)-exp(-ξ))/2

# ╔═╡ 78cd1016-f21f-4eb7-b81a-a831e2487452
# Provjerimo rubne uvjete za g
g(0.7,0),g(0.7,1)

# ╔═╡ 5be1ee41-b2c5-4bd6-880c-7a2cf79ee970
# Plot
ξ=0.8
plot(x->g(x,ξ),0,1)

# ╔═╡ dc91b239-5aa8-4e98-8103-e3b446d984d6
# Plot
x=0.6
plot(ξ->g(x,ξ),0,1)

# ╔═╡ 42e4fba9-a66c-4ec8-b798-1c6d880ba1d1
X=range(0,stop=1,length=101)
Ξ=range(0,stop=1,length=101)
surface(X,Ξ,g)

# ╔═╡ 226397c1-ef69-47c4-9d0a-31a1ee154f91
# Rješenje
u(x)=quadgk(ξ->g(x,ξ)*f(ξ),0,1)[1]+u₁*(e*(1-exp(1-x)) /
    (e-1)+exp(1-x))-u₀*(1-exp(1-x))/(e-1)

# ╔═╡ c90e0e61-ef01-4e85-a343-7414b1e59579
# Provjerimo zadovoljava li rješenje rubne uvjete
u(0),u(1)

# ╔═╡ 150e7125-d68a-4eb2-981e-18e3dc428d1f
# Nacrtajmo rješenje
plot(u,0,1)

# ╔═╡ 8dc78f2e-3b0c-46ec-a7b0-aee3e34b6345
md"""
__Zadatak.__ Nacrtajte rješenja za razne funkcije $f(x)$.
"""

# ╔═╡ 2477f9d3-562e-4560-9020-9798e4421374
md"""
## Bilinearni razvoj

Regularni Sturm-Liouvilleov problem glasi

$$
Lu=f, \quad a<x<b
$$

pri čemu operator $L$ uključuje diferencijalnu jednadžbu i rubne uvjete. Operator $L$ ima beskonačno svojstvenih vrijednosti $\lambda_n$ i pripadnih ortonormiranih svojstvenih funkcija $\phi_n(x)$, tako da vrijedi

$$
L\phi_n(x)=\lambda_n \phi_n(x).
$$

Funkcije $\phi_n$ imaju normu jedan i tvore bazu prostora. Greenova funkcija jednaka je 

$$
g(x,\xi)=\sum \frac{\phi_n(x)\phi_n(\xi)}{\lambda_n}.
$$

Dokažimo ovu tvrdnju. Funkcija $f$ se može prikazati pomoću elemenata baze:

$$
f(x)=\sum f_n \phi_n(x), \quad f_n=(f,\phi_n) \equiv \int\limits_a^b f(x)\phi_n(x)\, dx.
$$

Neka je

$$
u(x)=\sum c_n \phi_n(x)
$$
rješenje. Zbog linearnosti operatora $L$ vrijedi 

$$
Lu=L\big(\sum c_n\phi_n\big)=\sum c_n L(\phi_n)=\sum c_n \lambda_n \phi_n= \sum f_n \phi_n
$$

pa izjednačavanje koeficijenata uz $\phi_n$ daje

$$
c_n=\frac{f_n}{\lambda_n}.
$$

Dakle,

\begin{align*}
u(x)&=\sum c_n\phi_n(x)=\sum \frac{1}{\lambda_n} \bigg( \int_a^b f(\xi)\phi_n(\xi)\, d\xi \bigg) \phi_n(x)\\
&= \int_a^b \left(\sum \frac{\phi_n(x) \phi_n(\xi)}{\lambda_n} \right) f(\xi)\, d\xi
\end{align*}

pa je Greenova funkcija $g(x,\xi)$ dana izrazom u zagradama. 
"""

# ╔═╡ e0b01dd6-37e9-40b3-aea8-28b3209ac352
md"""
### Primjer Greenove funkcije u 2D
"""

# ╔═╡ aa92ae02-756e-4995-bcf4-0c8d5511aa4d


# ╔═╡ Cell order:
# ╟─8a4b7c05-081c-4913-8278-838f865bd2eb
# ╟─f8ebd0bc-794b-4018-b023-f72f890aba08
# ╟─61a891d8-6f5c-40a1-a6e1-335525f60c40
# ╟─d172800a-341c-4224-b933-51cda3121014
# ╠═649f28a3-f4f2-48d1-a1e4-b79b378c522c
# ╠═159a90f8-2bf4-40ad-9c31-de1692cea6c4
# ╠═c7483547-316d-4c82-affc-2aaf5a203cbf
# ╠═05c5f5ca-9931-4419-ae67-aa1b8236a599
# ╠═0db3d836-fd58-4597-a183-7043059d20dd
# ╟─3f1178ff-83a4-42d1-8c6a-fb39dd4de4ff
# ╟─be021115-1978-49fa-81a4-434c383883b0
# ╠═f122c3c3-1e29-47e7-8fc6-e7e013920962
# ╠═88107568-140d-4485-a680-ef1e74241f78
# ╠═d923f52a-b7b1-428e-9e81-40c45db0deba
# ╠═2b88038f-09e7-4722-9360-170c06fa1eaf
# ╠═80613a67-6843-44ed-be63-51c7cd62693e
# ╠═95f14d64-1c0a-4f67-a6fa-82616b57fc01
# ╠═e30fd6a4-7153-4689-b343-abd8898c5410
# ╟─119ba9d5-ced7-4ebd-8cf3-c127ef2f0aec
# ╟─2e4a1f96-005f-441d-beab-5b6a927c1e2e
# ╟─aa4614d8-04ee-4de9-979b-eab95a84f1e3
# ╟─6030f0b9-677e-4fde-b9d1-d01541628c04
# ╟─b63428a6-cb0a-4b8f-8688-eb0a72fc6774
# ╠═01fef681-189d-4963-89b7-f8d06b6ce588
# ╠═321cfc8a-e5ca-476d-8393-f8d23eb04f64
# ╠═f7409dfc-e679-411f-9fd9-79a3ff9fe835
# ╠═78cd1016-f21f-4eb7-b81a-a831e2487452
# ╠═5be1ee41-b2c5-4bd6-880c-7a2cf79ee970
# ╠═dc91b239-5aa8-4e98-8103-e3b446d984d6
# ╠═42e4fba9-a66c-4ec8-b798-1c6d880ba1d1
# ╠═226397c1-ef69-47c4-9d0a-31a1ee154f91
# ╠═c90e0e61-ef01-4e85-a343-7414b1e59579
# ╠═150e7125-d68a-4eb2-981e-18e3dc428d1f
# ╟─8dc78f2e-3b0c-46ec-a7b0-aee3e34b6345
# ╟─2477f9d3-562e-4560-9020-9798e4421374
# ╟─e0b01dd6-37e9-40b3-aea8-28b3209ac352
# ╠═aa92ae02-756e-4995-bcf4-0c8d5511aa4d
