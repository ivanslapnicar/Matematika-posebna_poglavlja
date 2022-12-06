### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 93f68ec6-615e-4deb-9362-ffa26cc52664
using SymPy, PlutoUI

# ╔═╡ 82ae49b6-f696-433b-b5aa-626569856b2b
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ b8d575dd-985d-4d5b-90fa-c1d47e797035
import_from(sympy)
# import sympy.core.add

# ╔═╡ 544b7449-12a0-48a1-ba7b-0dd429ee580b
md"""
# Integralne transformacije


__Integralna transformacija__ funkcije $f(t)$ na intervalu $[a,b]$ je

$$
F(s)=\int_a^b K(s,t)\,f(t) \, dt.$$

Funkcija $K(s,t)$ je __jezgra__ transformacije.

# Laplaceova transformacija

Za $a=0$, $b=\infty$ i $K(s,t)=e^{-st}$ imamo __Laplaceovu transformaciju__:

$$
(\mathcal{L}u)(s)\equiv U(s)=\int_0^\infty u(t)\,e^{-st}\, dt.$$

Funkcije koje su 

* po djelovima neprekidne na svakom konačnom intervalu i 
* koje su __eksponencijalnog rasta__, odnosno za koje postoje konstante  $M>0$ i $a>0$ takve da je

$$
\big|f(t)\big|\leq Me^{at}$$

sigurno imaju Laplaceovu transformaciju. Ovo su __dovoljni uvjeti__, ali __ne i nužni__. 

Laplaceova transformacija je __linearni operator__.

Lapleaceova transformacija ima __inverz__:

$$
\mathcal{L}^{-1}U(s)=u(t)=\frac{1}{2\pi i}\int_{a-i\infty}^{a+i\infty} U(s)\,e^{st}\, ds.$$

Parovi transformacija i njihovih inverza se nalaze u [tablicama](http://integral-table.com/downloads/LaplaceTable.pdf).

Posebno su važne formule za deriviranje:

$$
\begin{aligned}
(\mathcal{L}u')(s)&=sU(s)-u(0), \\
(\mathcal{L}u'')(s)&=s^2U(s)-su(0)-u'(0).
\end{aligned}$$
"""

# ╔═╡ 5bbb766f-dbf7-4555-834f-952c4c8aa180
md"""
## Problem početnih vrijednosti

Riješimo problem

$$
u''+u=0, \quad t>0, \quad u(0)=0, \quad u'(0)=1.$$

Laplaceova transformacija cijele jednadžbe daje

$$
s^2U(s)-su(0)-u'(0)+U(s)=0.$$

Uvrštavanje početnih uvjeta daje

$$
s^2U(s)-1+U(s)=0$$

pa je 

$$
U(s)=\frac{1}{1+s^2}.$$

Primjena inverzne transformacije daje rješenje

$$
u(t)=\mathcal{L}^{-1}\bigg(\frac{1}{1+s^2}\bigg)= \sin t.$$
"""

# ╔═╡ cf087642-2c04-4e73-b0da-9dfbaadc9ba0
x,t,s=symbols("x, t, s")

# ╔═╡ f7370b2f-3093-4bb3-8938-53ef6abff119
# u(t)=
inverse_laplace_transform(1/(1+s^2),s,t)

# ╔═╡ 2bd4bec2-8261-41ee-adc1-83fd4ed6f3dc
md"""
## Difuzija na polu-beskonačnom prostoru


Neka $u(x,t)$ daje koncentraciju kemikalije na polu-beskonačnom prostoru $x>0$ koji je u početku bez kemikalije. Neka tijekom vremena $t>0$ na rubu $x=0$ dajemo jediničnu koncentraciju kemikalije i želimo znati kako se kemikalija širi. Neka je difuzijska konstanata jednaka $1$.

Matematički model je

$$
\begin{aligned}
&u_t-u_{xx}=0, \quad x>0, t>0, \\
&u(x,0)=0, \quad x>0, \\
&u(0,t)=1, \quad t>0, \\
&u(x,t) \ \textrm{omeđena}.
\end{aligned}$$

Laplaceova transformacija jednadžbe po vremenu $t$, pri čemu se prostorna varijabla $x$ ne transformira, daje diferencijalnu jednadžbu po varijabli $x$:

$$
sU(x,s)-u(x,0)-U_{xx}(x,s)=0.$$

Počeni uvjet daje jednadžbu

$$
sU(x,s)-U_{xx}(x,s)=0.$$

Za detalje vidi [J. Logan, Applied Mathematics, 2nd ed., str. 226](https://books.google.hr/books/about/Applied_Mathematics.html?id=tD0ZAQAAIAAJ&redir_esc=y).
"""

# ╔═╡ 9f6d2d4a-fc21-49da-815f-e87d9cc995a9
# s=symbols("s",real=true,positive=true)

# ╔═╡ cf90ecb2-1c72-477d-a389-4b401d5f04c9
U = symbols("U", cls=SymFunction)

# ╔═╡ 5c5f2b22-5063-4ce4-b2be-444c084d6761
# Definirajmo diferencijalnu jednadžbu
diffeq = Eq(s*U(x)-diff(U(x), x, 2), 0)

# ╔═╡ f08cf50f-d0ac-4303-aacf-424873897128
# Riješimo jednadžbu
dsolve(diffeq)

# ╔═╡ dfc8e2f9-ffd7-42d0-adc9-c8459b3ecb99
md"""
Rješili smo jednadžbu po $x$ pa je varijabla $s$ konstanta. Zato su $C_1$ i $C_2$ funkcije od $s$,

$$
C_1 \equiv a(s), \quad C_2\equiv b(s),$$

odnosno,

$$
U(x,s)=a(s) e^{-\sqrt{s} x} + b(s)e^{\sqrt{s} x}.$$

Zato što želimo omeđeno rješenje, mora biti $b(s)=0$ pa je

$$
U(x,s)=a(s) e^{-\sqrt{s} x}.$$

Sada iskoristimo početni uvjet:

$$
U(0,s)=a(s)=\mathcal{L}(1)=\frac{1}{s}$$

pa je 

$$
U(x,s)=\frac{1}{s} e^{-\sqrt{s} x}.$$

Iz [tablice](http://integral-table.com/downloads/LaplaceTable.pdf) pod (33) slijedi

$$
u(x,s)=\mathop{\mathrm{erfc}} \left( \frac{x}{\sqrt{4t}}\right).$$
"""

# ╔═╡ 8f62449c-2a77-4a95-91d5-0f031a824d54
# U(0,s)=
laplace_transform(t^0,t,s)

# ╔═╡ f75b7b8e-7beb-4e62-8f7e-52a3ad86b1c7
# u(x,s)=
inverse_laplace_transform(exp(-sqrt(s)*x)/s,s,t)

# ╔═╡ 9608ade0-bf77-443a-8978-cff8624eba00
md"""
__Napomena.__ Vrijedi (vidi [Error function](https://en.wikipedia.org/wiki/Error_function)):

$$
\begin{aligned}
\mathop{\mathrm{erf}}(x)&=\frac{2}{\sqrt{\pi}} \int\limits_{0}^x e^{-t^2}dt,\\
\mathop{\mathrm{erfc}}(x)&=1-\mathop{\mathrm{erf}}(x)=\frac{2}{\sqrt{\pi}} \int\limits_{x}^\infty e^{-t^2}dt.
\end{aligned}$$
"""

# ╔═╡ 47c05b25-8da6-471e-b264-8b191367660a
md"""
# Fourierova transformacija

Za funkciju $u(x)$, $x\in\mathbb{R}$, definiramo __Fourierovu transformaciju__:

$$
(\mathcal{F} u)(\xi)\equiv \hat u(\xi) = \int\limits_{-\infty}^\infty u(x)\, \displaystyle e^{i\xi x} dx.$$

Ovo je integralna transformacija s granicama $a=-\infty$ i $b=\infty$ i jezgrom $K(\xi,x)=e^{i\xi x}$.

Fourierova transformacija postoji čim je $u$ apsolutno integrabilna funkcija, odnosno $\int\limits_{-\infty}^\infty \left| u(x)\right| dx < \infty$.

Promatrat ćemo funkcije __Schwartzove klase__ koje, zajedno s derivacijama, opadaju brže od bilo koje potencije:

$$
\mathcal{S}=\left\{ u\in C^\infty : \left\|\displaystyle \frac{d^k u}{dx^k}\right\| = 
\mathcal{O} \left( \displaystyle \frac{1}{\left|x\right|^N} \right), \ \left|x\right|\to \infty,\ k=0,1,2,3,\ldots, \ 
\forall N\in\mathbb{N} \right\}.$$

Inverzna Fourierova transformacija dana je formulom

$$
(\mathcal{F}^{-1}\hat u)(x)\equiv u(x) = \frac{1}{2\pi}\int\limits_{-\infty}^\infty \hat u(\xi) \,
\displaystyle e^{-i\xi x} d\xi.$$

Fourierove transformacije i inverzne Fourierove transformacije možemo naći u [tablicama](http://uspas.fnal.gov/materials/11ODU/FourierTransformPairs.pdf).

Posebno, za transformacije derivacija vrijedi

$$
(\mathcal{F} u^{(k)})(\xi)=(-i\xi)^k\, \hat u(\xi), \quad u\in \mathcal{S}.$$

__Konvolucija__ funkcija $u,v\in\mathcal{S}$ je funkcija

$$
(u\ast v)(x) =\int\limits_{-\infty}^\infty u(x-y)\,v(y)\, dy.$$

Vrijedi

$$
\mathcal{F}(u\ast v)(\xi)=\hat u(\xi)\,\hat v(\xi).$$

__Napomena.__ U navedenim tablicama korištena je jezgra $K(\xi,x)=e^{-i\xi x}$ pa parove treba prilagoditi tako da transformacije navedene u tablici daju $\hat u(-\xi)$.
"""

# ╔═╡ 2742f830-98f2-4b82-9b25-71c536ecb8cb
md"""
## Računanje Fourierove transformacije

Izračunajmo transformaciju funkcije $u(x)=e^{-ax^2}$, $a>0$. Vrijedi

$$
\hat u(\xi)=\int\limits_{-\infty}^\infty \displaystyle e^{-ax^2} \displaystyle e^{i\xi x} dx.$$

Deriviranje od znakom integrala daje

$$
\hat u'(\xi)=i\int\limits_{-\infty}^\infty \displaystyle e^{-ax^2} x\, \displaystyle e^{i\xi x} dx.$$

Primijenimo parcijalnu integraciju: neka je 

$$
u=e^{i\xi x}, \quad du=e^{i\xi x}\,i\xi\, dx,\\
dv=\int e^{-ax^2}x\, dx,\quad v=-\frac{1}{2a}e^{-ax^2}.$$
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

$$
\begin{aligned}
\hat u'(\xi)&= i\left[ -\frac{1}{2a} e^{-ax^2}e^{i\xi x}\bigg|_{-\infty}^\infty - \int\limits_{-\infty}^\infty  (-1)\frac{1}{2a} 
e^{-ax^2}e^{i\xi x} i\xi \, dx \right]\\
&= \frac{i^2\xi}{2a} \int\limits_{-\infty}^\infty e^{-ax^2}e^{i\xi x} \, dx = -\frac{\xi}{2a} \hat u(\xi).
\end{aligned}$$

Dobili smo linearnu diferencijalnu jednadžbu prvog reda

$$
\hat u'(\xi)=\frac{d\hat u}{d\xi}=-\frac{\xi}{2a} \hat u(\xi).$$

Separacija varijabli daje

$$
\frac{d\hat u}{\hat u}=-\frac{\xi}{2a} d\xi.$$

Integriranje daje

$$
\ln |\hat u| =-\frac{1}{2a}\frac{\xi^2}{2}=-\frac{\xi^2}{4a}$$

pa je 

$$
\hat u=C e^{-\xi^2/(4a)}.$$

Vrijedi

$$
\hat u(0)=C=\int\limits_{-\infty}^\infty e^{-ax^2}\, dx
=\sqrt{\frac{\pi}{a}}$$

pa je, konačno,

$$
\hat u(\xi)= \sqrt{\frac{\pi}{a}} e^{-\xi^2/(4a)}.$$

Dakle, Fourierova transformacija Gaussove funkcije je Gaussova funkcija.
"""

# ╔═╡ 1a373079-9103-49e9-8c52-e0aac28f00db
begin
	g=exp(-a*x^2)
	fourier_transform(g,x,s)
end

# ╔═╡ 325dd806-621d-490f-a8c8-3934750ced88
md"""
__Zašto se transformacije razlikuju?__

Paket `SymPy.jl` transformaciju definira koristeći jezgru $K(\xi,x)=e^{-2\pi i\xi x}$ pa parove treba prilagoditi.

Dokumentacija se nalazi na [docs.sympy.org](https://docs.sympy.org/latest/search.html?q=fourier_transform). 
"""

# ╔═╡ cfa135e0-3b93-4474-89ef-543ead2d1afc
md"""
## Problem rubnih vrijednosti

Za funkciju $f\in\mathcal{S}$ nađimo $u\in\mathcal{S}$ za koju je

$$
u''-u=f(x), \quad x\in\mathbb{R}.$$

Transfomacije jednadžbe daje

$$
(-i\xi)^2\hat u-\hat u=\hat f$$

pa je

$$
\hat u(\xi)=-\frac{1}{1+\xi^2} \,\hat f(\xi).$$

Iz tablica vidimo da je 

$$
\mathcal{F}^{-1}\left(\frac{1}{1+\xi^2}\right) = \frac{1}{2} e^{-|x|}$$

pa je po teoremu o konvoluciji

$$
u(x)=-\frac{1}{2} e^{-|x|} \ast f(x)=-\frac{1}{2} 
\int\limits_{-\infty}^\infty e^{-|x-y|}f(y)dy.$$
"""

# ╔═╡ 33ce7e8a-0d0a-4aee-921c-3c0c379347d6
begin
	# Usporedi s tablicama
	f₁=1/(1+s^2)
	inverse_fourier_transform(f₁,s,x)
end

# ╔═╡ 17105dd2-613c-4c92-b8a2-b269e84051cf
md"""
## Jednadžba difuzije

Riješimo problem

$$
u_t-ku_{xx}=0,\quad u(x,0)=f(x), \quad x\in\mathbb{R},\quad t>0.$$

Pretpostavljamo da je $f\in\mathcal{S}$. Fourierova transformacija jednadžbe po $x$ daje populacijsku jednadžbu

$$
\hat u_t=-\xi^2 k \hat u$$

pa je

$$
\hat u(\xi, t)= C e^{-\xi^2k t}.$$

Početni uvjet daje

$$
\hat u(\xi, 0)=C=\hat f(\xi)$$

pa je 

$$
\hat u(\xi, t)= \hat f(\xi)\, e^{-\xi^2k t}.$$

Iz tablica vidimo da je 

$$
\mathcal{F}^{-1}\left(e^{-\xi^2k t} \right) = \frac{1}{\sqrt{4\pi k t}}\, e^{-x^2/(4kt)}$$

pa je po teoremu o konvoluciji

$$
u(x,t)=\frac{1}{\sqrt{4\pi k t}}\int\limits_{-\infty}^\infty e^{-(x-y)^2/(4kt)}\,f(y)\,dy.$$

"""

# ╔═╡ 7466621e-3c11-44e1-bce8-68dafa5904a2
begin
	# Uz a=kt
	f₂=exp(-s^2*a)
	inverse_fourier_transform(f₂,s,x)
end

# ╔═╡ 86187835-80e6-453c-b071-04b0f457491b
md"""
## Laplaceova jednadžba

Riješimo problem

$$
u_{xx}+u_{yy}=0,\quad u(x,0)=f(x), \quad x\in\mathbb{R},\quad y>0,$$

uz dodatni uvijet da je rješenje omeđeno kada $y\to\infty$.

Pretpostavljamo da je $f\in\mathcal{S}$. Fourier-ova transformacija jednadžbe po $x$ daje jednadžbu

$$
\hat u_{yy}-\xi^2 \hat u=0,$$

pa je

$$
\hat u(\xi, y)= a(\xi)\, e^{-\xi y} + b(\xi)\, e^{\xi y}.$$

Dodatni uvijet omeđenosti rješenja povlači $b(\xi)=0$ pa je

$$
\hat u(\xi, y)= a(\xi)\, e^{-\xi y}.$$

Međutim, i ovo rješenje će rasti kada je $\xi<0$ pa stoga uzimamo

$$
\hat u(\xi, y)= a(\xi)\, e^{-|\xi| y}.$$

Rubni uvijet daje 

$$
\hat u(\xi, 0)=a(\xi)=\hat f(\xi)$$

pa je rješenje problema u transformiranoj domeni

$$
\hat u(\xi, y)= \hat f(\xi)\, e^{-|\xi|y}.$$

Iz Zadatka 4., str. 396, vidimo da je 

$$
\mathcal{F}^{-1}\left(e^{-|\xi|y} \right) = \frac{y}{\pi}\frac{1}{x^2+y^2}$$

odakle, koristeći teorem o konvoluciji, slijedi

$$
u(x,y)=\frac{y}{\pi}\frac{1}{x^2+y^2} \ast f = 
\frac{y}{\pi}\int\limits_{-\infty}^\infty \frac{f(\tau)\,d\tau}{(x-\tau)^2+y^2}.$$
"""

# ╔═╡ f7d68389-3690-4386-8023-113ba79d3a00
inverse_fourier_transform(exp(-abs(t)*s),t,x)

# ╔═╡ dba76381-ff96-4fca-91a0-057d3f83a4ac
md"""
Moramo provjeriti kako je definirana inverzna Fourierova transformacija u paketu `SymPy.jl`, vidi [docs.sympy.org](https://docs.sympy.org/latest/search.html?q=inverse_fourier_transform).
"""

# ╔═╡ 6f531aa2-ae12-4c06-a824-61cd3c5ee082
md"""
## Plancharelova jednakost

Vrijedi

$$
\int\limits_{-\infty}^\infty \big|u(x)\big|^2 dx =
\frac{1}{2\pi} \int\limits_{-\infty}^\infty \big|\hat u(\xi) \big|^2 d\xi.$$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
PlutoUI = "~0.7.21"
SymPy = "~1.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "b68904528fd538f1cb6a3fbc44d2abdc498f9e8e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.21"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "8f8d948ed59ae681551d184b93a256d0d5dd4eae"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═93f68ec6-615e-4deb-9362-ffa26cc52664
# ╠═82ae49b6-f696-433b-b5aa-626569856b2b
# ╠═b8d575dd-985d-4d5b-90fa-c1d47e797035
# ╟─544b7449-12a0-48a1-ba7b-0dd429ee580b
# ╟─5bbb766f-dbf7-4555-834f-952c4c8aa180
# ╠═cf087642-2c04-4e73-b0da-9dfbaadc9ba0
# ╠═f7370b2f-3093-4bb3-8938-53ef6abff119
# ╟─2bd4bec2-8261-41ee-adc1-83fd4ed6f3dc
# ╠═9f6d2d4a-fc21-49da-815f-e87d9cc995a9
# ╠═cf90ecb2-1c72-477d-a389-4b401d5f04c9
# ╠═5c5f2b22-5063-4ce4-b2be-444c084d6761
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
