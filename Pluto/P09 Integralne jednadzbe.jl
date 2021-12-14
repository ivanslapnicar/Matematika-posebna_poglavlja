### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ c2a942a6-79a8-42c9-a39a-7e6b6d5e9f63
using PlutoUI, SymPy, LinearAlgebra, Base.MathConstants

# ╔═╡ 421d8ff6-692a-4ce4-8caf-11baf4326183
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 1675701e-c972-4c76-a3ba-33fc3ff82ba2
md"""
# Integralne jednadžbe


Integralne jednadžbe su povezane s diferencijalnim jednadžbama. 

__Fredholmova jednadžba__ glasi

$$
\int\limits_a^b k(x,y)u(y)dy+\alpha(x)u(x)=f(x), \quad x\in[a,b].$$

__Volterra-ina jednadžba__ glasi

$$
\int\limits_a^x k(x,y)u(y)dy+\alpha(x)u(x)=f(x), \quad x\in[a,b].$$

Zadatak je naći funkciju $u$ koja zadovoljava jednadžbu.

Funkcija $k$ je __jezgra__. 

Ako je $k(x,y)=k(y,x)$, jezgra je __simetrična__. 

Ako je $f=0$, jednadžba je __homogena__.

Ako je $\alpha=0$, jednadžba je __prve vrste__, a inače je __druge vrste__.
"""

# ╔═╡ 10b7579b-23af-4e33-9a91-f3ceeb598ce7
md"""
Definirajmo __integralne operatore__:

__Fredholmov__: $Ku(x)=\int\limits_a^bk(x,y)u(y)dy$,

__Volerra-in__: $Ku(x)=\int\limits_a^x k(x,y)u(y)dy$.

Operatori su linearni.

U oba slučaja integralne jednadžbe možemo zapisati kao

$$
Ku+\alpha u=f$$

pa rješenje možemo tražiti i pomoću rješenje problema svojstvenih vrijednosti

$$
Ku=\lambda u.$$
"""

# ╔═╡ 42442f09-69a2-4c39-a9dd-69088bbf2a4f
md"""
__Skalarni produkt__ je

$$
(u,v)=\int\limits_a^b u(x)\cdot \overline{v(x)} dx.$$

Skalarni produkt je linearan i vrijedi:

*  $(u,v)=\overline{(v,u)}$,
*  $\|u\|=\sqrt{(u,u)}=\left(\int\limits_a^b \big| u(x)\big|^2dx\right)^{1/2}$,
*  $(u,u)=0 \Leftrightarrow u=0$,
*  ako je $(u,v)=0$, funkcije $u$ i $v$ su __ortogonalne__.
"""

# ╔═╡ 75f61491-9141-4440-bbd0-af2ffd450c71
md"""
## Kontrola inventara

Početna količina robe je $a$. Neka je $k(t)$ postotak robe koje nije prodana u trenutku $t$ nakon nabave.

Po kojoj stopi $u(t)$ treba naručivati robu ako želimo imati konstantnu zalihu?

Promotrimo vremenski interval $[\tau,\tau+\Delta \tau]$. Smatramo da su vrijednosti konstantne unutar kratkog intervala.
Ukupna količina nabavljene robe u tom intervalu je $u(\tau)\Delta \tau$.

Količina robe koja nije prodana u trenutku $t\in[\tau,\tau+\Delta \tau]$ je

$$
k(t-\tau)u(\tau)\Delta \tau.$$

Dakle, ukupna količina robe koja nije prodana u trenutku $t$ je zbroj dijela početne količine robe koji nije prodan u tenutku $t$ i neprodanog dijela robe koja je naručena do trenutka $t$:

$$
ak(t)+\int\limits_0^t k(t-\tau)u(\tau)d\tau.$$

Uvjet da je inventar konstantan daje Volterra-inu jednadžbu

$$
ak(t)+\int\limits_0^t k(t-\tau)u(\tau)d\tau=a.$$

Laplaceova transformacija i teorem o konvoluciji daju

$$
a\mathcal{L}(k) +\mathcal{L}(k) \cdot \mathcal{L}(u) =\frac{a}{s}$$

pa je 

$$
\mathcal{L}(u)=\bigg(\frac{a}{s}-a\mathcal{L}(k)\bigg) \frac{1}{\mathcal{L}(k)}=
a \bigg(\frac{1}{s \mathcal{L}(k)}-1\bigg).$$

Konačno,

$$
u(t)=a\mathcal{L}^{-1} \bigg(\frac{1}{s \mathcal{L}(k)} \bigg) -a\delta(t).$$

(Koristili smo $\mathcal{L}^{-1}(1)=\delta(t)$.)
"""

# ╔═╡ a8dbacd2-5ca2-4682-8007-c561548d0ef4
md"""
## Konvolucija

Riješimo jednadžbu

$$
u(x)=x-\int\limits_0^x (x-y) u(y) dy.$$

Laplaceova transformacija jednadžbe i teorem o konvoluciji daju

$$
\mathcal{L}(u)= \mathcal{L}(x)-\mathcal{L}(x)\cdot \mathcal{L}(u)=
\frac{1}{s^2}-\frac{1}{s^2}\mathcal{L}(u)$$

odnosno

$$
\mathcal{L}(u) \bigg(1+\frac{1}{s^2}\bigg)=\frac{1}{s^2}.$$

Dakle,

$$
\mathcal{L}(u)=\frac{1}{s^2}\cdot \frac{s^2}{1+s^2}=\frac{1}{1+s^2}$$

pa je 

$$
u(x)=\mathcal{L}^{-1} \bigg(\frac{1}{1+s^2}\bigg) =\sin x.$$
"""

# ╔═╡ 1cb24112-e22f-4519-9f5f-792cd774a6f6
x,y=symbols("x,y", real=true)

# ╔═╡ 1242db10-426b-46e7-98f7-0ee0a4ea3ce7
# Provjera
integrate((x-y)*sin(y),(y,0,x))

# ╔═╡ e6de4cf5-150c-4faa-b0d8-91414e79b0cc
md"""
# Prebacivanje DJ u IJ

Promotrimo problem početnih vrijednosti

$$
u'=f(x,u), \quad u(x_0)=u_0.$$

Ako je $u$ rješenje, onda (uz zamjenu varijabli) vrijedi

$$
u'(y)=f(y,u(y)), \quad \forall y.$$

Integriranje od $x_0$ do $x$ daje

$$
\int\limits_{x_0}^x u'(y) dy =\int\limits_{x_0}^x f(y,u(y)) dy.$$

Dakle, $u(x)$ je rješenje Volterra-ine jednadžbe

$$
u(x)=u(x_0)+\int\limits_{x_0}^x f(y,u(y)) dy.$$
"""

# ╔═╡ eb6624e0-b6da-4b3b-a638-5f77d58b5383
md"""
# Prebacivanje IJ u DJ

Neka je zadana integralna jednadžba

$$
u(x)=u_0+\int\limits_{x_0}^x f(y,u(y)) dy.$$

Uvrštavanje $x=x_0$ daje početni uvijet $u(x_0)=u_0$.

__Leibnitzova fomula__ glasi: ako je 

$$
I(\alpha)=\int\limits_{\phi(\alpha)}^{\psi(\alpha)} f(x,\alpha) dx.$$

onda je

$$
I'(\alpha)=\int\limits_{\phi(\alpha)}^{\psi(\alpha)} f_\alpha(x,\alpha) dx +f(\psi(\alpha),\alpha)\frac{d\psi}{d\alpha}-
f(\phi(\alpha),\alpha)\frac{d\phi}{d\alpha}.$$

Primjena Leibnitzove formule daje

$$
u'(x)=0+\int\limits_{x_0}^x f_x(y,u(y)) dy +f(x,u(x))\cdot x'-f(x_0,u(x_0))\cdot x_0'
= f(x,u(x)).$$

U zadnjoj jednakosti koristili smo činjenice:

$$
f_x(y,u(y))=0,\quad x'=1, \quad x_0'=0.$$
"""

# ╔═╡ 4613bdc3-2ffe-4e46-9798-37e6980fc711
md"""
# Prebacivanje DJ u IJ (II)

Promotrimo problem početnih vrijednosti

$$
u''+p(x)u'+q(x)u=f(x), \quad x>a, \quad u(a)=u_0, \ u'(a)=u_1.$$

Potreban nam je sljedeći rezultat:

__Lema.__
Neka je $f$ neprekidna funkcija za $x\geq a$. Onda je

$$
\int_a^x\int_a^s f(y)\,dy\,ds=\int_a^x f(y)(x-y)dy.$$

_Dokaz:_ Označimo $F(s)=\int_a^s f(y)\,dy$. Parcijalna integracija daje:

$$
\begin{aligned}
\int_a^x\int_a^s f(y)\,dy\,ds&=\int_a^x F(s)\, ds \quad (u=F(s), \ dv=ds) \\
&= sF(s)\bigg\vert_a^x -\int_a^x sF'(s)\, ds\\
&= xF(x)-aF(a)-\int_a^x sF'(s)\, ds.
\end{aligned}$$

Očito je $F(a)=0$, a Leibnitzova formula daje

$$
F'(s)=\int_a^s f_s(y)\,dy +f(s)\frac{ds}{ds}-f(a)\frac{da}{ds}=f(s).$$

Dakle,

$$
\int_a^x\int_a^s f(y)\,dy\,ds=xF(x)-\int_a^x sf(s)\,ds
=x\int_a^x f(y)\,dy-\int_a^x yf(y)\, dy.$$

_Q.E.D._

Vratimo se zadanom problemu početnih vrijednosti. 

Vrijedi (vidi [J. Logan, Applied Mathematics, str. 233](https://www.amazon.com/Applied-Mathematics-J-David-Logan/dp/0471746622))

$$
u''=-p(x)u'-q(x)u+f(x)$$

pa je

$$
\int_a^x u''(y)\, dy = -\int_a^x p(y)u'(y)\,dy -\int_a^x(q(y)u(y)-f(y))\, dy.$$

Integriranje lijeve strane i uvrštavanje početnih uvjeta te parcijalna integracija prvog integrala na desnoj strani daje

$$
u'(x)-u_1=-p(x)u(x)+p(a)u_0-\int_a^x [(q(y)-p'(y))u(y)-f(y)]\,dy$$


Ponovna integracija daje

$$
\begin{aligned}
u(x)-u_0=&-\int_a^x p(y)u(y)\, dy\\ & -\int_a^x \int_a^s [(q(y)-p'(y))u(y)-f(y)]\,dy\, ds
+(p(a) u_0+u_1)(x-a).
\end{aligned}$$

Iz leme konačno slijedi

$$
\begin{aligned}
u(x)=&-\int_a^x \{ p(y)+(x-y) [q(y)-p'(y)]\} u(y)\, dy\\ &-\int_a^x (x-y)f(y)\, dy
+(p(a) u_0+u_1)(x-a)+u_0,
\end{aligned}$$

što je Volterra-ina jednadžba oblika

$$
u(x)=\int_a^x k(x,y)u(y)\, dy +F(x).$$
"""

# ╔═╡ 3d5ffc95-c761-481d-87fe-0b694abf2920
md"""
# Fredholmova jednadžba

Promotrimo jednostavniju jednadžbu ($\alpha(x)=-\lambda$)

$$
Ku-\lambda u=f$$

i to za __separabilnu jezgru__,

$$
k(x,y)=\sum_{j=1}^{n}\alpha_j(x)\beta_j(y).$$

__Teorem.__ Za separabilnu jezgru $k$ vrijedi:

* za $\lambda=0$ vrijedi:
    * ako $f$ nije linearna kombinacija od $\alpha_i$, jednadžba nema rješenja,
    * ako je $f$ linearna kombinacija od $\alpha_i$, jednadžba ima beskonačno rješenja,
* za $\lambda\neq 0$ izračunajmo matricu skalarnih produkata $A_{ij}=(\beta_i,\alpha_j)$ i njene svojstvene vrijednosti. Vrijedi:
    * ako je $\lambda$ svojstvena vrijednost od $A$, jednadžba ili nema rješenje ili ima beskonačno rješenja,
    * ako $\lambda$ nije svojstvena vrijednost od $A$, jednadžba ima jedinstveno rješenje.
"""

# ╔═╡ 3fc91534-bcdb-4b5c-ac15-e79d4b524abf
md"""
## Primjer 1

Analizirajmo zadnji slučaj teorema.

Uvrštavanje jezgre u jednadžbu daje

$$
\int\limits_a^b \sum_{j=1}^{n}\alpha_j(x)\beta_j(y) u(y)dy-\lambda u(x)=f(x)$$

odnosno

$$
\sum_{j=1}^{n}\alpha_j(x) \int\limits_a^b \beta_j(y) u(y)dy -\lambda u(x)=f(x).$$

Uz oznake $c_j=(\beta_j,u)$ i $c=\begin{pmatrix} c_1 \\ \vdots \\ c_n\end{pmatrix}$ imamo

$$
\sum_{j=1}^{n}\alpha_j(x) c_j-\lambda u(x)=f(x) \tag{1}.$$

Pomnožimo ovu jednadžbu s $\beta_i(x)$ i integrirajmo od $a$ do $b$:

$$
\sum_{j=1}^{n}(\beta_i,\alpha_j) c_j -\lambda c_i=(\beta_i,f), \quad i=1,2,\ldots, n. \tag{2}$$

Uz oznaku $F=\begin{pmatrix} (\beta_1,f) \\ \vdots \\ (\beta_n,f) \end{pmatrix}$, dobili smo sustav linearnih jednadžbi

$$
(A-\lambda I)c=F.$$

Rješenje problema je

$$
u(x)=\frac{1}{\lambda} \bigg(-f(x)+\sum_{j=1}^{n}\alpha_j(x)c_j\bigg).
\tag{3}$$
"""

# ╔═╡ cfbf6735-06cc-4576-a652-ffba02f81b4e
md"""
## Primjer 2

Riješimo jednadžbu

$$
\int\limits_0^1 (1-3xy)u(y)dy-2 u(x) =e^x. \tag{4}$$

Pripadni operator je

$$
Ku(x)=\int\limits_0^1 (1-3xy)u(y)dy.$$

Jezgra je separabilna:

$$
k(x,y)=1-3xy=\alpha_1(x)\beta_1(y)+\alpha_2(x)\beta_2(y),$$

gdje je 

$$
\alpha_1(x)=1,\ \beta_1(y)=1,\ \alpha_2(x)=-3x, \ \beta_2(y)=y.$$

Matrica $A$ je 

$$
A=\begin{pmatrix} \int_0^1 1\cdot 1 \, dx & \int_0^1 1\cdot (-3x)\, dx\\
\int_0^1 x\cdot 1\, dx & \int_0^1 (-3x)\cdot x \,dx
\end{pmatrix} =\begin{pmatrix} 1 &-\frac{3}{2} \\ 
\frac{1}{2} &  -1 \end{pmatrix}.$$

Svojstvene vrijednosti matrice $A$ su rješenja jednadžbe

$$
\det(A-\lambda I)=0.$$
"""

# ╔═╡ 9278f978-90fd-4fae-8cbd-be9fa70b40c8
A=[1//1 -3//2; 1//2 -1]

# ╔═╡ 82da59e7-e988-448c-a101-8ee253f93f77
eye(n)=Matrix{Rational}(I,n,n)

# ╔═╡ 6b1b340e-fbd6-4826-9a1d-ced8996cee06
det(A-x*eye(2))

# ╔═╡ 88c36114-f86a-4c6f-a70c-0beac1df4553
md"""
Dakle, svojstvene vrijednosti matrice $A$ su 

$$
\lambda_1=\frac{1}{2}, \quad \lambda_2=-\frac{1}{2}.$$
"""

# ╔═╡ 7bc5579a-cccb-4ff5-b3ee-796c11b85770
eigen(A)

# ╔═╡ 968faa0c-c110-4040-aee6-c407a956458a
e

# ╔═╡ 5cb6c5d7-7f05-48ed-aa44-d9f56b417a92
md"""
U zadatku je $\lambda=2$ različita od svojstvenih vrijednosti matrice $A$ pa je prema teoremu rješenje jedinstveno.
"""

# ╔═╡ 348b6f15-8c8c-4281-8aae-889ce005941c
# Izračunajmo F
F=[integrate(e^x,(x,0,1)); integrate(e^x*x,(x,0,1))]

# ╔═╡ 7409348b-53c3-4e35-863e-cf5fa6c4503f
# Izračunajmo c kao rješenje sustava
c=(A-2*I)\F

# ╔═╡ f54a54b2-12bb-45f3-a7db-5251f3c355c4
# Rješenje polaznog problema prema (3)
u(x)=1/2*(-exp(x)+1*c[1]+(-3*x)*c[2])

# ╔═╡ d8b30280-0ec7-42dc-a503-a5842828f38f
# Na primjer
u(y)

# ╔═╡ 97c84300-4078-426f-a650-1935089a24ab
# Provjera - uvrštavanje u lijevu stranu (4) treba dati exp(x)
ex=integrate((1-3*x*y)*u(y),(y,0,1))-2*u(x)

# ╔═╡ a29e21b0-3a6e-4f22-ad04-213c3b246e99
# !!!!!!!!!!!!!
simplify(ex)

# ╔═╡ 62c3a9fe-291a-466c-a741-80207e5a465e
md"""
### Dodatak

Izračunajmo svojstvene vektore matrice $A$ i svojstvene funkcije operatora $K$

Svojstveni vektori matrice $A$ su rješenja jednadžbi

$$
A\begin{pmatrix}x \\ y\end{pmatrix}= \lambda \begin{pmatrix}x \\ y\end{pmatrix},$$

odnosno

$$
\begin{pmatrix} 1 &-\frac{3}{2} \\ 
\frac{1}{2} &  -1 \end{pmatrix}\begin{pmatrix}x \\ y\end{pmatrix}= \lambda \begin{pmatrix}x \\ y\end{pmatrix}.$$

Uvrštavanje $\lambda_1$ i $\lambda_2$ daje redom

$$
v_1=\begin{pmatrix}3 \\ 1\end{pmatrix}, \quad
v_2=\begin{pmatrix}1 \\ 1\end{pmatrix}.$$

Uvrštavanje $f=0$ u (1) daje jednadžbu $Ku=\lambda u$, odnosno

$$
\sum_{j=1}^{n}\alpha_j(x) c_j=\lambda u(x) \tag{5}$$

a uvrštavanje $f=0$ u (2) daje jednadžbu 

$$
\sum_{j=1}^{n}(\beta_i,\alpha_j) c_j -\lambda c_i=(\beta_i,f), \quad i=1,2,\ldots, n,$$

odnosno $Ac=\lambda c$.

Zadnja jednadžba će biti zadovoljena kada je $\lambda$ svojstvena vrijednost matrice $A$ ($\lambda_1$ ili $\lambda_2$), a $c$ pripadni svojstveni vektor, $v_1$ ili $v_2$. Iz (5) slijedi da su svojstvene vrijednosti matrice $A$ ujedno i svojstvene vrijednosti operatora $K$ i da su pripadne svojstvene funkcije dane formulom

$$
u(x)=\frac{1}{\lambda}\sum_{j=1}^{n}\alpha_j(x) c_j,$$

s time da u ovoj fromuli možemo uzeti i $\lambda=1$ jer svojstvene funkcije ne ovise o skaliranju. U našem slučaju su svojstvene funkcije opratora $K$ jednake 

$$
\begin{aligned}
u_1(x)&= \alpha_1(x)[v_1]_1+\alpha_2(x)[v_1]_2=1\cdot 3+(-3x)\cdot 1= 3(1-x), \\
u_2(x)&=\alpha_1(x)[v_2]_1 +\alpha_2(x)[v_2]_2=1\cdot 1+(-3x)\cdot 1= 1-3x.
\end{aligned}$$

Vrijedi $(u_1,u_2)=0$, odnosno $u_1\perp u_2$.
"""

# ╔═╡ 24d21607-ad66-4ce7-a30c-8a9d1f03d8ef
# Provjera Ku=λu za λ=1/2, u=3(1-x)
integrate((1-3*x*y)*3*(1-y),(y,0,1))

# ╔═╡ 3fccce23-64d2-4aa4-a700-c3674f177e06
# Provjera Ku=λu za λ=-1/2, u=1-3x 
integrate((1-3*x*y)*(1-3*y),(y,0,1))

# ╔═╡ 3e974d75-7adc-461e-82a8-304a2befb7df
# Ortogonalnost svojstvenih funkcija
integrate((1-x)*(1-3*x),(x,0,1))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
PlutoUI = "~0.7.23"
SymPy = "~1.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "8f8d948ed59ae681551d184b93a256d0d5dd4eae"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═c2a942a6-79a8-42c9-a39a-7e6b6d5e9f63
# ╠═421d8ff6-692a-4ce4-8caf-11baf4326183
# ╟─1675701e-c972-4c76-a3ba-33fc3ff82ba2
# ╟─10b7579b-23af-4e33-9a91-f3ceeb598ce7
# ╟─42442f09-69a2-4c39-a9dd-69088bbf2a4f
# ╟─75f61491-9141-4440-bbd0-af2ffd450c71
# ╟─a8dbacd2-5ca2-4682-8007-c561548d0ef4
# ╠═1cb24112-e22f-4519-9f5f-792cd774a6f6
# ╠═1242db10-426b-46e7-98f7-0ee0a4ea3ce7
# ╟─e6de4cf5-150c-4faa-b0d8-91414e79b0cc
# ╟─eb6624e0-b6da-4b3b-a638-5f77d58b5383
# ╟─4613bdc3-2ffe-4e46-9798-37e6980fc711
# ╟─3d5ffc95-c761-481d-87fe-0b694abf2920
# ╟─3fc91534-bcdb-4b5c-ac15-e79d4b524abf
# ╟─cfbf6735-06cc-4576-a652-ffba02f81b4e
# ╠═9278f978-90fd-4fae-8cbd-be9fa70b40c8
# ╠═82da59e7-e988-448c-a101-8ee253f93f77
# ╠═6b1b340e-fbd6-4826-9a1d-ced8996cee06
# ╟─88c36114-f86a-4c6f-a70c-0beac1df4553
# ╠═7bc5579a-cccb-4ff5-b3ee-796c11b85770
# ╠═968faa0c-c110-4040-aee6-c407a956458a
# ╟─5cb6c5d7-7f05-48ed-aa44-d9f56b417a92
# ╠═348b6f15-8c8c-4281-8aae-889ce005941c
# ╠═7409348b-53c3-4e35-863e-cf5fa6c4503f
# ╠═f54a54b2-12bb-45f3-a7db-5251f3c355c4
# ╠═d8b30280-0ec7-42dc-a503-a5842828f38f
# ╠═97c84300-4078-426f-a650-1935089a24ab
# ╠═a29e21b0-3a6e-4f22-ad04-213c3b246e99
# ╟─62c3a9fe-291a-466c-a741-80207e5a465e
# ╠═24d21607-ad66-4ce7-a30c-8a9d1f03d8ef
# ╠═3fccce23-64d2-4aa4-a700-c3674f177e06
# ╠═3e974d75-7adc-461e-82a8-304a2befb7df
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
