### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 544b7449-12a0-48a1-ba7b-0dd429ee580b
md"""
# Integralne transformacije

---

__Integralna transformacija__ funkcije $f(t)$ na intervalu $[a,b]$ je

$$
F(s)=\int_a^b K(s,t)\,f(t) \, dt.
$$

Funkcija $K(s,t)$ je __jezgra__ transformacije.

## Laplaceova transformacija

Za $a=0$, $b=\infty$ i $K(s,t)=e^{-st}$ imamo __Laplaceovu transformaciju__:

$$
(\mathcal{L}u)(s)\equiv U(s)=\int_0^\infty u(t)\,e^{-st}\, dt.
$$

Funkcije koje su 

* po djelovima neprekidne na svakom konačnom intervalu i 
* koje su __eksponencijalnog rasta__, odnosno za koje postoje konstante  $M>0$ i $a>0$ takve da je

$$
\big|f(t)\big|\leq Me^{at}
$$

sigurno imaju Laplaceovu transformaciju. Ovo su __dovoljni uvjeti__, ali __ne i nužni__. 

Laplaceova transformacija je __linearni operator__.

Lapleaceova transformacija ima __inverz__:

$$
\mathcal{L}^{-1}U(s)=u(t)=\frac{1}{2\pi i}\int_{a-i\infty}^{a+i\infty} U(s)\,e^{st}\, ds.
$$

Parovi transformacija i njihovih inverza se nalaze u [tablicama](http://integral-table.com/downloads/LaplaceTable.pdf).

Posebno su važne formule za deriviranje:

\begin{align*}
(\mathcal{L}u')(s)&=sU(s)-u(0), \\
(\mathcal{L}u'')(s)&=s^2U(s)-su(0)-u'(0).
\end{align*}
"""

# ╔═╡ 5bbb766f-dbf7-4555-834f-952c4c8aa180
md"""
### Problem početnih vrijednosti

Riješimo problem

$$
u''+u=0, \quad t>0, \quad u(0)=0, \quad u'(0)=1.
$$

Laplaceova transformacija cijele jednadžbe daje

$$
s^2U(s)-su(0)-u'(0)+U(s)=0.
$$

Uvrštavanje početnih uvjeta daje

$$
s^2U(s)-1+U(s)=0
$$

pa je 

$$
U(s)=\frac{1}{1+s^2}.
$$

Primjena inverzne transformacije daje rješenje

$$
u(t)=\mathcal{L}^{-1}\bigg(\frac{1}{1+s^2}\bigg)= \sin t.
$$
"""

# ╔═╡ b8d575dd-985d-4d5b-90fa-c1d47e797035
using SymPy
import_from(sympy)

# ╔═╡ cf087642-2c04-4e73-b0da-9dfbaadc9ba0
x,t,s=symbols("x, t, s")

# ╔═╡ f7370b2f-3093-4bb3-8938-53ef6abff119
u=inverse_laplace_transform(1/(1+s^2),s,t)

# ╔═╡ 2bd4bec2-8261-41ee-adc1-83fd4ed6f3dc
md"""
### Primjer difuzije


Neka $u(x,t)$ daje koncentraciju kemikalije na polu-beskonačnom prostoru $x>0$ koji je u početku bez kemikalije.
Neka tijekom vremena $t>0$ na rubu $x=0$ dajemo jediničnu koncentraciju kemikalije i želimo znati kako se kemikalija širi. 
Neka je difuzijska konstanata jednaka $1$.

Matematički model je

\begin{align*}
&u_t-u_{xx}=0, \quad x>0, t>0, \\
&u(x,0)=0, \quad x>0, \\
&u(0,t)=1, \quad t>0, \\
&u(x,t) \ \textrm{omeđena}.
\end{align*}

Laplaceova transformacija jednadžbe po vremenu $t$, pri čemu se prostorna varijabla $x$ ne transformira, daje
diferencijalnu jednadžbu po varijabli $x$:

$$
sU(x,s)-u(x,0)-U_{xx}(x,s)=0.
$$

Počeni uvjet daje jednadžbu

$$
sU(x,s)-U_{xx}(x,s)=0.
$$

Za detalje vidi [J. Logan, Applied Mathematics, 2nd ed., str. 226][JL97].

[JL97]: https://books.google.hr/books/about/Applied_Mathematics.html?id=tD0ZAQAAIAAJ&redir_esc=y "J. Logan, 'Applied Mathematics', 2nd ed., Wiley, New York, 1997"
"""

# ╔═╡ 9f6d2d4a-fc21-49da-815f-e87d9cc995a9
s=symbols("s",real=true,positive=true)

# ╔═╡ f7af0ec3-7a64-4827-895b-9dff33e55624
U = symbols("U", cls=SymFunction)
# U=SymFunction('U')
diffeq = Eq(s*U(x)-diff(U(x), x, 2), 0)

# ╔═╡ f08cf50f-d0ac-4303-aacf-424873897128
# Riješimo jednadžbu
dsolve(diffeq)

# ╔═╡ dfc8e2f9-ffd7-42d0-adc9-c8459b3ecb99
md"""
Rješili smo jednadžbu po $x$ pa je varijabla $s$ konstanta. Zato su $C_1$ i $C_2$ funkcije od $s$,

$$
C_1 \equiv a(s), \quad C_2\equiv b(s),
$$

odnosno,

$$
U(x,s)=a(s) e^{-\sqrt{s} x} + b(s)e^{\sqrt{s} x}.
$$

Zato što želimo omeđeno rješenje, mora biti $b(s)=0$ pa je

$$
U(x,s)=a(s) e^{-\sqrt{s} x}.
$$

Sada iskoristimo početni uvjet:

$$
U(0,s)=a(s)=\mathcal{L}(1)=\frac{1}{s}
$$

pa je 

$$
U(x,s)=\frac{1}{s} e^{-\sqrt{s} x}.
$$

Iz [tablice](http://integral-table.com/downloads/LaplaceTable.pdf) pod (33) slijedi

$$
u(x,s)=\mathop{\mathrm{erfc}} \left( \frac{x}{\sqrt{4t}}\right).
$$
"""

# ╔═╡ 8f62449c-2a77-4a95-91d5-0f031a824d54
f=laplace_transform(t^0,t,s)

# ╔═╡ f75b7b8e-7beb-4e62-8f7e-52a3ad86b1c7
inverse_laplace_transform(exp(-sqrt(s)*x)/s,s,t)

# ╔═╡ 9608ade0-bf77-443a-8978-cff8624eba00
md"""
__Napomena.__ Vrijedi (vidi [Error function](https://en.wikipedia.org/wiki/Error_function)):

\begin{align*}
\mathop{\mathrm{erf}}(x)&=\frac{2}{\sqrt{\pi}} \int\limits_{0}^x e^{-t^2}dt,\\
\mathop{\mathrm{erfc}}(x)&=1-\mathop{\mathrm{erf}}(x)=\frac{2}{\sqrt{\pi}} \int\limits_{x}^\infty e^{-t^2}dt.
\end{align*}
"""

# ╔═╡ 47c05b25-8da6-471e-b264-8b191367660a
md"""
## Fourierova transformacija

Za funkciju $u(x)$, $x\in\mathbb{R}$, definiramo __Fourierovu transformaciju__:

$$
(\mathcal{F} u)(\xi)\equiv \hat u(\xi) = \int\limits_{-\infty}^\infty u(x)\, \displaystyle e^{i\xi x} dx.
$$

Ovo je integralna transformacija s granicama $a=-\infty$ i $b=\infty$ i jezgrom
$K(\xi,x)=e^{i\xi x}$.

Fourierova transformacija postoji čim je $u$ apsolutno integrabilna funkcija, odnosno
$\int\limits_{-\infty}^\infty \left| u(x)\right| dx < \infty$.

Promatrat ćemo funkcije __Schwartzove klase__ koje, zajedno s derivacijama, 
opadaju brže od bilo koje potencije:

$$
\mathcal{S}=\left\{ u\in C^\infty : \left\|\displaystyle \frac{d^k u}{dx^k}\right\| = 
\mathcal{O} \left( \displaystyle \frac{1}{\left|x\right|^N} \right), \ \left|x\right|\to \infty,\ k=0,1,2,3,\ldots, \ 
\forall N\in\mathbb{N} \right\}.
$$

Inverzna Fourierova transformacija dana je formulom

$$
(\mathcal{F}^{-1}\hat u)(x)\equiv u(x) = \frac{1}{2\pi}\int\limits_{-\infty}^\infty \hat u(\xi) \,
\displaystyle e^{-i\xi x} d\xi.
$$

Fourierove transformacije i inverzne Fourierove transformacije možemo naći u [tablicama](http://uspas.fnal.gov/materials/11ODU/FourierTransformPairs.pdf).

Posebno, za transformacije derivacija vrijedi

$$
(\mathcal{F} u^{(k)})(\xi)=(-i\xi)^k\, \hat u(\xi), \quad u\in \mathcal{S}.
$$

__Konvolucija__ funkcija $u,v\in\mathcal{S}$ je funkcija

$$
(u\ast v)(x) =\int\limits_{-\infty}^\infty u(x-y)\,v(y)\, dy.
$$

Vrijedi

$$
\mathcal{F}(u\ast v)(\xi)=\hat u(\xi)\,\hat v(\xi).
$$

__Napomena.__ U navedenim tablicama korištena je jezgra $K(\xi,x)=e^{-i\xi x}$ pa parove treba prilagoditi tako da transformacije navedene u tablici daju $\hat u(-\xi)$.
"""

# ╔═╡ 2742f830-98f2-4b82-9b25-71c536ecb8cb
md"""
### Računanje Fourierove transformacije

Izračunajmo transformaciju funkcije $u(x)=e^{-ax^2}$, $a>0$. Vrijedi

$$
\hat u(\xi)=\int\limits_{-\infty}^\infty \displaystyle e^{-ax^2} \displaystyle e^{i\xi x} dx.
$$

Deriviranje od znakom integrala daje

$$
\hat u'(\xi)=i\int\limits_{-\infty}^\infty \displaystyle e^{-ax^2} x\, \displaystyle e^{i\xi x} dx.
$$

Primijenimo parcijalnu integraciju: neka je 

$$
u=e^{i\xi x}, \quad du=e^{i\xi x}\,i\xi\, dx,\\
dv=\int e^{-ax^2}x\, dx,\quad v=-\frac{1}{2a}e^{-ax^2}.
$$
"""

# ╔═╡ 308a9d74-cb69-49a4-9c5a-eb2696e0dba7
a=symbols("a", positive=true, real=true)

# ╔═╡ b937b75e-651e-4a2a-ba05-11e8f4c6bd7d
f=x*exp(-a*x^2)

# ╔═╡ fdf015d7-023a-4855-9b5d-654a5ac1b6d2
integrate(f,x)

# ╔═╡ 8ca0a65d-ea94-42d8-a73c-0a4ea2378014
md"""
Sada je

\begin{align*}
\hat u'(\xi)&= i\left[ -\frac{1}{2a} e^{-ax^2}e^{i\xi x}\bigg|_{-\infty}^\infty - \int\limits_{-\infty}^\infty  (-1)\frac{1}{2a} 
e^{-ax^2}e^{i\xi x} i\xi \, dx \right]\\
&= \frac{i^2\xi}{2a} \int\limits_{-\infty}^\infty e^{-ax^2}e^{i\xi x} \, dx = -\frac{\xi}{2a} \hat u(\xi).
\end{align*}

Dobili smo linearnu diferencijalnu jednadžbu prvog reda

$$
\hat u'(\xi)=\frac{d\hat u}{d\xi}=-\frac{\xi}{2a} \hat u(\xi).
$$

Separacija varijabli daje

$$
\frac{d\hat u}{\hat u}=-\frac{\xi}{2a} d\xi.
$$

Integriranje daje

$$
\ln |\hat u| =-\frac{1}{2a}\frac{\xi^2}{2}=-\frac{\xi^2}{4a}
$$

pa je 

$$
\hat u=C e^{-\xi^2/(4a)}.
$$

Vrijedi

$$
\hat u(0)=C=\int\limits_{-\infty}^\infty e^{-ax^2}\, dx
=\sqrt{\frac{\pi}{a}}
$$

pa je, konačno,

$$
\hat u(\xi)= \sqrt{\frac{\pi}{a}} e^{-\xi^2/(4a)}.
$$

Dakle, Fourierova transformacija Gaussove funkcije je Gaussova funkcija.
"""

# ╔═╡ 1a373079-9103-49e9-8c52-e0aac28f00db
g=exp(-a*x^2)
fourier_transform(g,x,s)

# ╔═╡ 325dd806-621d-490f-a8c8-3934750ced88
md"""
__Zašto se transformacije razlikuju?__

Paket `SymPy.jl` transformaciju definira koristeći jezgru $K(\xi,x)=e^{-2\pi i\xi x}$ pa parove treba prilagoditi.

Dokumentacija se nalazi na adresi https://docs.sympy.org/latest/search.html?q=fourier_transform.
"""

# ╔═╡ cfa135e0-3b93-4474-89ef-543ead2d1afc
md"""
### Problem rubnih vrijednosti

Za funkciju $f\in\mathcal{S}$ nađimo $u\in\mathcal{S}$ za koju je

$$
u''-u=f(x), \quad x\in\mathbb{R}.
$$

Transfomacije jednadžbe daje

$$
(-i\xi)^2\hat u-\hat u=\hat f
$$

pa je

$$
\hat u(\xi)=-\frac{1}{1+\xi^2} \,\hat f(\xi).
$$

Iz tablica vidimo da je 

$$
\mathcal{F}^{-1}\left(\frac{1}{1+\xi^2}\right) = \frac{1}{2} e^{-|x|}
$$

pa je po teoremu o konvoluciji

$$
u(x)=-\frac{1}{2} e^{-|x|} \ast f(x)=-\frac{1}{2} 
\int\limits_{-\infty}^\infty e^{-|x-y|}f(y)dy.
$$
"""

# ╔═╡ 33ce7e8a-0d0a-4aee-921c-3c0c379347d6
# Usporedi s tablicama
f=1/(1+s^2)
inverse_fourier_transform(f,s,x)

# ╔═╡ 17105dd2-613c-4c92-b8a2-b269e84051cf
md"""
### Jednadžba difuzije

Riješimo problem

$$
u_t-ku_{xx}=0,\quad u(x,0)=f(x), \quad x\in\mathbb{R},\quad t>0.
$$

Pretpostavljamo da je $f\in\mathcal{S}$. Fourierova transformacija jednadžbe po $x$ daje populacijsku 
jednadžbu

$$
\hat u_t=-\xi^2 k \hat u
$$

pa je

$$
\hat u(\xi, t)= C e^{-\xi^2k t}.
$$

Početni uvjet daje

$$
\hat u(\xi, 0)=C=\hat f(\xi)
$$

pa je 

$$
\hat u(\xi, t)= \hat f(\xi)\, e^{-\xi^2k t}.
$$

Iz tablica vidimo da je 

$$
\mathcal{F}^{-1}\left(e^{-\xi^2k t} \right) = \frac{1}{\sqrt{4\pi k t}}\, e^{-x^2/(4kt)}
$$

pa je po teoremu o konvoluciji

$$
u(x,t)=\frac{1}{\sqrt{4\pi k t}}\int\limits_{-\infty}^\infty e^{-(x-y)^2/(4kt)}\,f(y)\,dy.
$$

"""

# ╔═╡ 7466621e-3c11-44e1-bce8-68dafa5904a2
# Uz a=kt
f=exp(-s^2*a)
inverse_fourier_transform(f,s,x)

# ╔═╡ 86187835-80e6-453c-b071-04b0f457491b
md"""
### Laplaceova jednadžba

Riješimo problem

$$
u_{xx}+u_{yy}=0,\quad u(x,0)=f(x), \quad x\in\mathbb{R},\quad y>0,
$$

uz dodatni uvijet da je rješenje omeđeno kada $y\to\infty$.

Pretpostavljamo da je $f\in\mathcal{S}$. Fourier-ova transformacija jednadžbe po $x$ daje
jednadžbu

$$
\hat u_{yy}-\xi^2 \hat u=0,
$$

pa je

$$
\hat u(\xi, y)= a(\xi)\, e^{-\xi y} + b(\xi)\, e^{\xi y}.
$$

Dodatni uvijet omeđenosti rješenja povlači $b(\xi)=0$ pa je

$$
\hat u(\xi, y)= a(\xi)\, e^{-\xi y}.
$$

Međutim, i ovo rješenje će rasti kada je $\xi<0$ pa stoga uzimamo

$$
\hat u(\xi, y)= a(\xi)\, e^{-|\xi| y}.
$$

Rubni uvijet daje 
$$
\hat u(\xi, 0)=a(\xi)=\hat f(\xi)
$$

pa je rješenje problema u transformiranoj domeni

$$
\hat u(\xi, y)= \hat f(\xi)\, e^{-|\xi|y}.
$$

Iz Zadatka 4., str. 396, vidimo da je 

$$
\mathcal{F}^{-1}\left(e^{-|\xi|y} \right) = \frac{y}{\pi}\frac{1}{x^2+y^2}
$$

odakle, koristeći teorem o konvoluciji, slijedi

$$
u(x,y)=\frac{y}{\pi}\frac{1}{x^2+y^2} \ast f = 
\frac{y}{\pi}\int\limits_{-\infty}^\infty \frac{f(\tau)\,d\tau}{(x-\tau)^2+y^2}.
$$
"""

# ╔═╡ f7d68389-3690-4386-8023-113ba79d3a00
inverse_fourier_transform(exp(-abs(t)*s),t,x)

# ╔═╡ dba76381-ff96-4fca-91a0-057d3f83a4ac
md"""
Moramo provjeriti kako je definirana inverzna Fourierova transformacija u paketu `SymPy.jl`, vidi 
http://docs.sympy.org/latest/search.html?q=inverse_fourier_transform.
"""

# ╔═╡ 6f531aa2-ae12-4c06-a824-61cd3c5ee082
md"""
### Plancharelova jednakost

Vrijedi

$$
\int\limits_{-\infty}^\infty \big|u(x)\big|^2 dx =
\frac{1}{2\pi} \int\limits_{-\infty}^\infty \big|\hat u(\xi) \big|^2 d\xi.
$$
"""

# ╔═╡ edc86c0e-0a1b-4f7b-a177-97149d2a1931


# ╔═╡ Cell order:
# ╟─544b7449-12a0-48a1-ba7b-0dd429ee580b
# ╟─5bbb766f-dbf7-4555-834f-952c4c8aa180
# ╠═b8d575dd-985d-4d5b-90fa-c1d47e797035
# ╠═cf087642-2c04-4e73-b0da-9dfbaadc9ba0
# ╠═f7370b2f-3093-4bb3-8938-53ef6abff119
# ╟─2bd4bec2-8261-41ee-adc1-83fd4ed6f3dc
# ╠═9f6d2d4a-fc21-49da-815f-e87d9cc995a9
# ╠═f7af0ec3-7a64-4827-895b-9dff33e55624
# ╠═f08cf50f-d0ac-4303-aacf-424873897128
# ╟─dfc8e2f9-ffd7-42d0-adc9-c8459b3ecb99
# ╠═8f62449c-2a77-4a95-91d5-0f031a824d54
# ╠═f75b7b8e-7beb-4e62-8f7e-52a3ad86b1c7
# ╟─9608ade0-bf77-443a-8978-cff8624eba00
# ╟─47c05b25-8da6-471e-b264-8b191367660a
# ╟─2742f830-98f2-4b82-9b25-71c536ecb8cb
# ╠═308a9d74-cb69-49a4-9c5a-eb2696e0dba7
# ╠═b937b75e-651e-4a2a-ba05-11e8f4c6bd7d
# ╠═fdf015d7-023a-4855-9b5d-654a5ac1b6d2
# ╟─8ca0a65d-ea94-42d8-a73c-0a4ea2378014
# ╠═1a373079-9103-49e9-8c52-e0aac28f00db
# ╟─325dd806-621d-490f-a8c8-3934750ced88
# ╟─cfa135e0-3b93-4474-89ef-543ead2d1afc
# ╠═33ce7e8a-0d0a-4aee-921c-3c0c379347d6
# ╟─17105dd2-613c-4c92-b8a2-b269e84051cf
# ╠═7466621e-3c11-44e1-bce8-68dafa5904a2
# ╟─86187835-80e6-453c-b071-04b0f457491b
# ╠═f7d68389-3690-4386-8023-113ba79d3a00
# ╟─dba76381-ff96-4fca-91a0-057d3f83a4ac
# ╟─6f531aa2-ae12-4c06-a824-61cd3c5ee082
# ╠═edc86c0e-0a1b-4f7b-a177-97149d2a1931
