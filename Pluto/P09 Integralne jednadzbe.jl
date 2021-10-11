### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 1675701e-c972-4c76-a3ba-33fc3ff82ba2
md"""
# Integralne jednadžbe

---

Integralne jednadžbe su povezane s diferencijalnim jednadžbama. 

__Fredholmova jednadžba__ glasi

$$
\int\limits_a^b k(x,y)u(y)dy+\alpha(x)u(x)=f(x), \quad x\in[a,b].
$$

__Volterra-ina jednadžba__ glasi

$$
\int\limits_a^x k(x,y)u(y)dy+\alpha(x)u(x)=f(x), \quad x\in[a,b].
$$

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
Ku+\alpha u=f
$$

pa rješenje možemo tražiti i pomoću rješenje problema svojstvenih vrijednosti

$$
Ku=\lambda u.
$$
"""

# ╔═╡ 42442f09-69a2-4c39-a9dd-69088bbf2a4f
md"""
__Skalarni produkt__ je

$$
(u,v)=\int\limits_a^b u(x)\cdot \overline{v(x)} dx.
$$

Skalarni produkt je linearan i vrijedi:

* $(u,v)=\overline{(v,u)}$,
* $\|u\|=\sqrt{(u,u)}=\left(\int\limits_a^b \big| u(x)\big|^2dx\right)^{1/2}$,
* $(u,u)=0 \Leftrightarrow u=0$,
* ako je $(u,v)=0$, funkcije $u$ i $v$ su __ortogonalne__.
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
k(t-\tau)u(\tau)\Delta \tau.
$$

Dakle, ukupna količina robe koja nije prodana u trenutku $t$ je zbroj dijela početne količine robe koji nije prodan u tenutku $t$ i neprodanog dijela robe koja je naručena do trenutka $t$:

$$
ak(t)+\int\limits_0^t k(t-\tau)u(\tau)d\tau.
$$

Uvjet da je inventar konstantan daje Volterra-inu jednadžbu

$$
ak(t)+\int\limits_0^t k(t-\tau)u(\tau)d\tau=a.
$$

Laplaceova transformacija i teorem o konvoluciji daju

$$
a\mathcal{L}(k) +\mathcal{L}(k) \cdot \mathcal{L}(u) =\frac{a}{s}
$$

pa je 

$$
\mathcal{L}(u)=\bigg(\frac{a}{s}-a\mathcal{L}(k)\bigg) \frac{1}{\mathcal{L}(k)}=
a \bigg(\frac{1}{s \mathcal{L}(k)}-1\bigg).
$$

Konačno,

$$
u(t)=a\mathcal{L}^{-1} \bigg(\frac{1}{s \mathcal{L}(k)} \bigg) -a\delta(t).
$$

(Koristili smo $\mathcal{L}^{-1}(1)=\delta(t)$.)
"""

# ╔═╡ a8dbacd2-5ca2-4682-8007-c561548d0ef4
md"""
## Konvolucija

Riješimo jednadžbu

$$
u(x)=x-\int\limits_0^x (x-y) u(y) dy.
$$

Laplaceova transformacija jednadžbe i teorem o konvoluciji daju

$$
\mathcal{L}(u)= \mathcal{L}(x)-\mathcal{L}(x)\cdot \mathcal{L}(u)=
\frac{1}{s^2}-\frac{1}{s^2}\mathcal{L}(u)
$$

odnosno

$$
\mathcal{L}(u) \bigg(1+\frac{1}{s^2}\bigg)=\frac{1}{s^2}.
$$

Dakle,

$$
\mathcal{L}(u)=\frac{1}{s^2}\cdot \frac{s^2}{1+s^2}=\frac{1}{1+s^2}
$$

pa je 

$$
u(x)=\mathcal{L}^{-1} \bigg(\frac{1}{1+s^2}\bigg) =\sin x.
$$
"""

# ╔═╡ 670126d3-0f01-4100-979f-b3041eaa65f5
using SymPy

# ╔═╡ 1cb24112-e22f-4519-9f5f-792cd774a6f6
x,y=symbols("x,y", real=true)

# ╔═╡ 1242db10-426b-46e7-98f7-0ee0a4ea3ce7
# Provjera
integrate((x-y)*sin(y),(y,0,x))

# ╔═╡ e6de4cf5-150c-4faa-b0d8-91414e79b0cc
md"""
## Prebacivanje DJ u IJ

Promotrimo problem početnih vrijednosti

$$
u'=f(x,u), \quad u(x_0)=u_0.
$$

Ako je $u$ rješenje, onda (uz zamjenu varijabli) vrijedi

$$
u'(y)=f(y,u(y)), \quad \forall y.
$$

Integriranje od $x_0$ do $x$ daje

$$
\int\limits_{x_0}^x u'(y) dy =\int\limits_{x_0}^x f(y,u(y)) dy.
$$

Dakle, $u(x)$ je rješenje Volterra-ine jednadžbe

$$
u(x)=u(x_0)+\int\limits_{x_0}^x f(y,u(y)) dy.
$$
"""

# ╔═╡ eb6624e0-b6da-4b3b-a638-5f77d58b5383
md"""
## Prebacivanje IJ u DJ

Neka je zadana integralna jednadžba

$$
u(x)=u_0+\int\limits_{x_0}^x f(y,u(y)) dy.
$$

Uvrštavanje $x=x_0$ daje početni uvijet $u(x_0)=u_0$.

__Leibnitzova fomula__ glasi: ako je 

$$
I(\alpha)=\int\limits_{\phi(\alpha)}^{\psi(\alpha)} f(x,\alpha) dx.
$$

onda je

$$
I'(\alpha)=\int\limits_{\phi(\alpha)}^{\psi(\alpha)} f_\alpha(x,\alpha) dx +f(\psi(\alpha),\alpha)\frac{d\psi}{d\alpha}-
f(\phi(\alpha),\alpha)\frac{d\phi}{d\alpha}.
$$

Primjena Leibnitzove formule daje

$$
u'(x)=0+\int\limits_{x_0}^x f_x(y,u(y)) dy +f(x,u(x))\cdot x'-f(x_0,u(x_0))\cdot x_0'
= f(x,u(x)).
$$ 

U zadnjoj jednakosti koristili smo činjenice:

$$
f_x(y,u(y))=0,\quad x'=1, \quad x_0'=0.
$$
"""

# ╔═╡ 4613bdc3-2ffe-4e46-9798-37e6980fc711
md"""
## Primjer - Prebacivanje DJ u IJ (II)

Promotrimo problem početnih vrijednosti

$$
u''+p(x)u'+q(x)u=f(x), \quad x>a, \quad u(a)=u_0, \ u'(a)=u_1.
$$

Potreban nam je sljedeći rezultat:

__Lema.__
Neka je $f$ neprekidna funkcija za $x\geq a$. Onda je

$$
\int_a^x\int_a^s f(y)\,dy\,ds=\int_a^x f(y)(x-y)dy.
$$

_Dokaz:_ Označimo $F(s)=\int_a^s f(y)\,dy$. Parcijalna integracija daje:

\begin{align*}
\int_a^x\int_a^s f(y)\,dy\,ds&=\int_a^x F(s)\, ds \quad (u=F(s), \ dv=ds) \\
&= sF(s)\bigg\vert_a^x -\int_a^x sF'(s)\, ds\\
&= xF(x)-aF(a)-\int_a^x sF'(s)\, ds.
\end{align*}

Očito je $F(a)=0$, a Leibnitzova formula daje

$$
F'(s)=\int_a^s f_s(y)\,dy +f(s)\frac{ds}{ds}-f(a)\frac{da}{ds}=f(s).
$$

Dakle,

$$
\int_a^x\int_a^s f(y)\,dy\,ds=xF(x)-\int_a^x sf(s)\,ds
=x\int_a^x f(y)\,dy-\int_a^x yf(y)\, dy.
$$

_Q.E.D._

Vratimo se zadanom problemu početnih vrijednosti. Vrijedi

$$
u''=-p(x)u'-q(x)u+f(x)
$$

pa je

$$
\int_a^x u''(y)\, dy = -\int_a^x p(y)u'(y)\,dy -\int_a^x(q(y)u(y)-f(y))\, dy.
$$

Integriranje lijeve strane i uvrštavanje početnih uvjeta te parcijalna integracija prvog integrala na desnoj strani daje

$$
u'(x)-u_1=-p(x)u(x)+p(a)u_0-\int_a^x [(q(y)-p'(y))u(y)-f(y)]\,dy
$$


Ponovna integracija daje

\begin{align*}
u(x)-u_0=&-\int_a^x p(y)u(y)\, dy\\ & -\int_a^x \int_a^s [(q(y)-p'(y))u(y)-f(y)]\,dy\, ds
+(p(a) u_0+u_1)(x-a).
\end{align*}

Iz leme konačno slijedi

$$
u(x)=-\int_a^x \{ p(y)+(x-y) [q(y)-p'(y)]\} u(y)\, dy-\int_a^x (x-y)f(y)\, dy
+(p(a) u_0+u_1)(x-a)+u_0,
$$

što je Volterra-ina jednadžba oblika

$$
u(x)=\int_a^x k(x,y)u(y)\, dy +F(x).
$$


(Vidi [J. Logan, Applied Mathematics, str. 233][Log06].)

[Log06]: https://www.amazon.com/Applied-Mathematics-J-David-Logan/dp/0471746622 "J. Logan, 'Applied Mathematics', 3rd Edition, Wiley and Sons, New York, 2006"
"""

# ╔═╡ 3d5ffc95-c761-481d-87fe-0b694abf2920
md"""
## Fredholmova jednadžba

Promotrimo jednostavniju jednadžbu ($\alpha(x)=-\lambda$)

$$
Ku-\lambda u=f
$$

i to za __separabilnu jezgru__,

$$
k(x,y)=\sum_{j=1}^{n}\alpha_j(x)\beta_j(y).
$$

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
### Primjer

Analizirajmo zadnji slučaj teorema.

Uvrštavanje jezgre u jednadžbu daje

$$
\int\limits_a^b \sum_{j=1}^{n}\alpha_j(x)\beta_j(y) u(y)dy-\lambda u(x)=f(x)
$$

odnosno

$$
\sum_{j=1}^{n}\alpha_j(x) \int\limits_a^b \beta_j(y) u(y)dy -\lambda u(x)=f(x).
$$

Uz oznake $c_j=(\beta_j,u)$ i $c=\begin{pmatrix} c_1 \\ \vdots \\ c_n\end{pmatrix}$ imamo

$$
\sum_{j=1}^{n}\alpha_j(x) c_j-\lambda u(x)=f(x) \tag{1}.
$$

Pomnožimo ovu jednadžbu s $\beta_i(x)$ i integrirajmo od $a$ do $b$:

$$
\sum_{j=1}^{n}(\beta_i,\alpha_j) c_j -\lambda c_i=(\beta_i,f), \quad i=1,2,\ldots, n. \tag{2}
$$

Uz oznaku $F=\begin{pmatrix} (\beta_1,f) \\ \vdots \\ (\beta_n,f) \end{pmatrix}$, dobili smo sustav linearnih jednadžbi

$$
(A-\lambda I)c=F.
$$

Rješenje problema je

$$
u(x)=\frac{1}{\lambda} \bigg(-f(x)+\sum_{j=1}^{n}\alpha_j(x)c_j\bigg).
\tag{3}
$$
"""

# ╔═╡ cfbf6735-06cc-4576-a652-ffba02f81b4e
md"""
### Primjer

Riješimo jednadžbu

$$
\int\limits_0^1 (1-3xy)u(y)dy-2 u(x) =e^x. \tag{4}
$$

Pripadni operator je

$$
Ku(x)=\int\limits_0^1 (1-3xy)u(y)dy.
$$

Jezgra je separabilna:

$$
k(x,y)=1-3xy=\alpha_1(x)\beta_1(y)+\alpha_2(x)\beta_2(y),
$$

gdje je 

$$
\alpha_1(x)=1,\ \beta_1(y)=1,\ \alpha_2(x)=-3x, \ \beta_2(y)=y.
$$

Matrica $A$ je 

$$
A=\begin{pmatrix} \int_0^1 1\cdot 1 \, dx & \int_0^1 1\cdot (-3x)\, dx\\
\int_0^1 x\cdot 1\, dx & \int_0^1 (-3x)\cdot x \,dx
\end{pmatrix} =\begin{pmatrix} 1 &-\frac{3}{2} \\ 
\frac{1}{2} &  -1 \end{pmatrix}.
$$

Svojstvene vrijednosti matrice $A$ su rješenja jednadžbe

$$
\det(A-\lambda I)=0.
$$ 
"""

# ╔═╡ 9278f978-90fd-4fae-8cbd-be9fa70b40c8
A=[1//1 -3//2; 1//2 -1]

# ╔═╡ 82da59e7-e988-448c-a101-8ee253f93f77
using LinearAlgebra
eye(n)=Matrix{Rational}(I,n,n)

# ╔═╡ 6b1b340e-fbd6-4826-9a1d-ced8996cee06
det(A-x*eye(2))

# ╔═╡ 593f6546-d08e-487f-bf84-4badd352699e
simplify(det(A-x*eye(2)))

# ╔═╡ 88c36114-f86a-4c6f-a70c-0beac1df4553
md"""
Dakle, svojstvene vrijednosti matrice $A$ su 

$$
\lambda_1=\frac{1}{2}, \quad \lambda_2=-\frac{1}{2}.
$$
"""

# ╔═╡ 7bc5579a-cccb-4ff5-b3ee-796c11b85770
eigen(A)

# ╔═╡ c71ae1f2-26b5-48e3-bae4-456eb7ad8ecd
using Base.MathConstants

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
A\begin{pmatrix}x \\ y\end{pmatrix}= \lambda \begin{pmatrix}x \\ y\end{pmatrix},
$$

odnosno

$$
\begin{pmatrix} 1 &-\frac{3}{2} \\ 
\frac{1}{2} &  -1 \end{pmatrix}\begin{pmatrix}x \\ y\end{pmatrix}= \lambda \begin{pmatrix}x \\ y\end{pmatrix}.
$$

Uvrštavanje $\lambda_1$ i $\lambda_2$ daje redom

$$
v_1=\begin{pmatrix}3 \\ 1\end{pmatrix}, \quad
v_2=\begin{pmatrix}1 \\ 1\end{pmatrix}.
$$

Uvrštavanje $f=0$ u (1) daje jednadžbu $Ku=\lambda u$, odnosno

$$
\sum_{j=1}^{n}\alpha_j(x) c_j=\lambda u(x) \tag{5}
$$

a uvrštavanje $f=0$ u (2) daje jednadžbu 

$$
\sum_{j=1}^{n}(\beta_i,\alpha_j) c_j -\lambda c_i=(\beta_i,f), \quad i=1,2,\ldots, n,$$

odnosno

$$
Ac=\lambda c.
$$

Zadnja jednadžba će biti zadovoljena kada je $\lambda$ svojstvena vrijednost matrice $A$ ($\lambda_1$ ili $\lambda_2$), a $c$ pripadni svojstveni vektor, $v_1$ ili $v_2$. Iz (5) slijedi da su svojstvene vrijednosti matrice $A$ ujedno i svojstvene vrijednosti operatora $K$ i da su pripadne svojstvene funkcije dane formulom

$$
u(x)=\frac{1}{\lambda}\sum_{j=1}^{n}\alpha_j(x) c_j,
$$

s time da u ovoj fromuli možemo uzeti i $\lambda=1$ jer svojstvene funkcije ne ovise o skaliranju. 
U našem slučaju su svojstvene funkcije opratora $K$ jednake 

$$
u_1(x)= \alpha_1(x)[v_1]_1+\alpha_2(x)[v_1]_2=1\cdot 3+(-3x)\cdot 1= 3(1-x), \\
u_2(x)=\alpha_1(x)[v_2]_1 +\alpha_2(x)[v_2]_2=1\cdot 1+(-3x)\cdot 1= 1-3x.
$$

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

# ╔═╡ e85f1a18-811b-4a53-aa3c-ac88cb09583b


# ╔═╡ Cell order:
# ╟─1675701e-c972-4c76-a3ba-33fc3ff82ba2
# ╟─10b7579b-23af-4e33-9a91-f3ceeb598ce7
# ╟─42442f09-69a2-4c39-a9dd-69088bbf2a4f
# ╟─75f61491-9141-4440-bbd0-af2ffd450c71
# ╟─a8dbacd2-5ca2-4682-8007-c561548d0ef4
# ╠═670126d3-0f01-4100-979f-b3041eaa65f5
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
# ╠═593f6546-d08e-487f-bf84-4badd352699e
# ╟─88c36114-f86a-4c6f-a70c-0beac1df4553
# ╠═7bc5579a-cccb-4ff5-b3ee-796c11b85770
# ╠═c71ae1f2-26b5-48e3-bae4-456eb7ad8ecd
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
# ╠═e85f1a18-811b-4a53-aa3c-ac88cb09583b
