### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 93f68ec6-615e-4deb-9362-ffa26cc52664
using SymPyPythonCall, PlutoUI

# ╔═╡ 82ae49b6-f696-433b-b5aa-626569856b2b
TableOfContents(title="📚 Sadržaj", aside=true)

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
x=symbols("x")

# ╔═╡ 696ef8d6-b7d5-4fec-b924-d2011fc45974
t=symbols("t")

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
s=symbols("s",real=true,positive=true)

# ╔═╡ f7370b2f-3093-4bb3-8938-53ef6abff119
# u(t)=
sympy.inverse_laplace_transform(1/(1+s^2),s,t)

# ╔═╡ cf90ecb2-1c72-477d-a389-4b401d5f04c9
U = symbols("U",cls=sympy.Function);

# ╔═╡ 876316b4-de44-42a7-b895-ffd43cc16940
typeof(U)

# ╔═╡ ec5f5023-d7ea-42ad-b674-71f042bd2f80
# Ovo trenutno ne radi?
#=
begin
	U(x)
	U(x).diff(x)
	# Definirajmo diferencijalnu jednadžbu
	diffeq = Eq(s*U(x)-U(x).diff(x,x), 0)
	dsolve(diffeq)
end
=#

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
sympy.laplace_transform(t^0,t,s)

# ╔═╡ f75b7b8e-7beb-4e62-8f7e-52a3ad86b1c7
# u(x,s)=
factor(sympy.inverse_laplace_transform(exp(-√(s)*x)/s,s,t))

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
	sympy.fourier_transform(g,x,s)
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
	sympy.inverse_fourier_transform(f₁,s,x)
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
	sympy.inverse_fourier_transform(f₂,s,x)
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
sympy.inverse_fourier_transform(exp(-abs(t)*s),t,x)

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
SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[compat]
PlutoUI = "~0.7.60"
SymPyPythonCall = "~0.4.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "6ea967922a745d366e6df01d4ffc28ead5894e9b"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CommonEq]]
git-tree-sha1 = "6b0f0354b8eb954cdba708fb262ef00ee7274468"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "TOML"]
git-tree-sha1 = "8f7faef2ca039ee068cd971a80ccd710d23fb2eb"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.23"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "1d322381ef7b087548321d3f878cb4c9bd8f8f9b"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.1"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "011cab361eae7bcd7d278f0a7a00ff9c69000c51"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.14"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "REPL", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "06a778ec6d6e76b0c2fb661436a18bce853ec45f"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.23"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SymPyCore]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "bef92ec4c31804bdc9c44cb00eaf0348eac383fb"
uuid = "458b697b-88f0-4a86-b56b-78b75cfb3531"
version = "0.2.5"

    [deps.SymPyCore.extensions]
    SymPyCoreTermInterfaceExt = "TermInterface"

    [deps.SymPyCore.weakdeps]
    TermInterface = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"

[[deps.SymPyPythonCall]]
deps = ["CommonEq", "CommonSolve", "CondaPkg", "LinearAlgebra", "PythonCall", "SpecialFunctions", "SymPyCore"]
git-tree-sha1 = "a8e887c6a810ce59f68e7640bb12732d5c32f092"
uuid = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
version = "0.4.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "b4a5a3943078f9fd11ae0b5ab1bdbf7718617945"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.5.8+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═93f68ec6-615e-4deb-9362-ffa26cc52664
# ╠═82ae49b6-f696-433b-b5aa-626569856b2b
# ╟─544b7449-12a0-48a1-ba7b-0dd429ee580b
# ╟─5bbb766f-dbf7-4555-834f-952c4c8aa180
# ╠═cf087642-2c04-4e73-b0da-9dfbaadc9ba0
# ╠═696ef8d6-b7d5-4fec-b924-d2011fc45974
# ╠═f7370b2f-3093-4bb3-8938-53ef6abff119
# ╟─2bd4bec2-8261-41ee-adc1-83fd4ed6f3dc
# ╠═9f6d2d4a-fc21-49da-815f-e87d9cc995a9
# ╠═cf90ecb2-1c72-477d-a389-4b401d5f04c9
# ╠═876316b4-de44-42a7-b895-ffd43cc16940
# ╠═ec5f5023-d7ea-42ad-b674-71f042bd2f80
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
