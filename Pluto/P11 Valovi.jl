### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ da811f07-9dda-4191-813e-4d11aa3b7357
md"""
# Valovi

---

Neka $x$ označava položaj, $t$ vrijeme, a $u$ veličinu poremećaja. Funkcija

$$
u(x,t)=f(x-ct)
$$

je __desni val__, odnosno val koji putuje udesno brzinom $c$. Vrijedi

$$
u_x=f'\cdot \frac{d(x-ct)}{dx}=f', \quad 
u_t=f'\cdot \frac{d(x-ct)}{dt}=-cf',
$$

što daje __advekcijsku__  parcijalnu diferencijalnu jednadžbu (_advekcija_ je transport tvari usmjerenim gibanjem)

$$
u_t+cu_x=0.
$$

Dakle, ako u trenutku $t=0$ poremećaj ima profil $u(x,0)=f(x)$ (početni uvjet), onda
se u trenutku $t>0$ poremećaj pomaknu udesno za $ct$ jedinica duljine.

Na sličan način __lijevi val__, $u(x,t)=f(x+ct)$, rješava jednadžbu $u_t-cu_x=0$.

Za lijevi val vrijedi

$$
u_x=f',\quad u_{xx}=f'',\quad u_t=cf',\quad u_{tt}=c^2f'',
$$

a za desni val vrijedi

$$
u_x=f',\quad u_{xx}=f'',\quad u_t=-cf',\quad u_{tt}=c^2f''.
$$

Dakle, i lijevi i desni val zadovoljavaju __valnu jednadžbu__

$$
u_{tt}-c^2u_{xx}=0.
$$

Obrnuto, zaključujemo da je rješenje valne jednadžbe proizlazi iz funkcije $f(x)$
kao superpozicija (zbroj) lijevog i desnog vala. 

"""

# ╔═╡ 7c761460-692c-4977-bdc6-9c9c0ade9ea3
md"""
__Primjer.__ Rješenje problema

$$
u_{tt}=u_{xx},\quad u(x,0)=\sin \pi x,
$$

je superpozicija lijevog vala 

$$
u_L(x,t)=\sin \pi (x+1\cdot t),
$$

i desnog vala 

$$
u_R(x,t)=\sin \pi (x-1\cdot t),
$$

odnosno,

$$
u(x,t)=\frac{1}{2}(u_L+u_R)=\frac{1}{2} [\sin \pi (x+t)+\sin \pi (x-t)]=
\sin \pi x \cos \pi t.
$$

"""

# ╔═╡ 3b596ad7-5489-4591-93ed-6613711d20c4
using Plots
#plotly()
gr()

# ╔═╡ 02409215-4359-4254-9ccc-b2c8b69bd20e
# Desni val - probajte razna vremena t od 0 do 3
t=1
u(x)=sin(π*(x-t))
s=1.7
v(x)=sin(π*(x-s))
plot(x->x,[u,v],0,2*π,xlabel="x",ylabel="u",aspect_ratio=1)

# ╔═╡ 9f51dc48-ebb2-43ee-bfc9-863e5e5e0014
# Rješenje - probajte razna vremena t od 0 do 3
t₁=1
u₁(x)=sin(pi*x)*cos(pi*t₁)
t₂=1.2
u₂(x)=sin(pi*x)*cos(pi*t₂)
t₃=1.6
u₃(x)=sin(pi*x)*cos(pi*t₃)
plot(x->x,[u₁,u₂,u₃],0,2*pi,xlabel="x",ylabel="u")

# ╔═╡ 83501b08-ce5b-4f1f-9e1d-6663190c4936
# Cijelo rješenje
X=range(0,stop=2*pi,length=101)
T=range(0,stop=1,length=101)
u(x,t)=sin(pi*x)*cos(pi*t)
surface(X,T,u,xlabel="x",ylabel="t",seriescolor=:blues)

# ╔═╡ 2f923b97-42c2-453d-9ff5-5f3ffa9f633f
md"""
## Sinusoidalni valovi

Definirajmo još neke pojmove vezane uz sinusoidalnu valnu funkciju

$$
u(x,t)=A\cos(kx-\omega t). 
$$

$A$ je __amplituda__. 

$k$ je __valni broj__, odnosno broj oscilacija u $2\pi$ jedinica prostora u danom trenutku $t$.

$\omega$ je __kutna frekvencija__, odnosno broj oscilacija u $2\pi$ jedinica vremena na zadanom mjestu $x$.

$\lambda=\displaystyle\frac{2\pi}{k}$ je __valna duljina__, odnosno udaljenost između susjednih kresta.

$P=\displaystyle\frac{2\pi}{\omega}$ je vremenski period, odnosno vrijeme nakon kojeg se na mjestu $x$ ponavlja isti obrazac poremećaja.

$c=\displaystyle\frac{\omega}{k}$ je __brzina__ kojom val putuje udesno.
"""

# ╔═╡ 95fb8bf4-ada6-4247-9a1e-bdabb67570c7
# Probajte za razne vrijednosti t, A, k i ω
t=0.5
A=2.0
k=0.6
ω=0.8
u(x)=A*cos(k*x-ω*t)
plot(u,0,4*pi,xlabel="x",ylabel="u")

# ╔═╡ 4e19a9f9-ee64-4490-9b1f-aa62067ed634
2*pi/k

# ╔═╡ eb656cbd-b950-4188-a829-9d00f42e37a9
# Probajte za razne vrijednosti x, A, k i ω
x=0.4
A=2.0
k=0.6
ω=0.8
u(t)=A*cos(k*x-ω*t)
plot(u,0,4*pi,xlabel="t",ylabel="u")

# ╔═╡ f255a3e4-0c73-42eb-a364-1fb401d04bca
2*pi/ω

# ╔═╡ de1707e0-f4b4-4477-92ec-b3a7abc6cbaa
md"""
## Relacija disperzije

Umjesto realnog rješenje $u(x,t)=A\cos(kx-\omega t)$, možemo koristiti __Eulerov oblik__ 

$$u(x,t)=Ae^{i(kx-\omega t)}$$ 

te promatrati realni ili imaginarni dio kako bi se dobilo realno rješenje.
Uvrštavanje Eulerovog oblika u parcijalnu diferencijalnu jednadžbu daje nam __relaciju disperzije__, odnosno vezu između valnog broja i kutne frekvencije,

$$
\omega=\omega(k).
$$

Vrijedi sljedeće: 

* ako je $\omega(k)$ kompleksna, kažemo da je jednadžba __difuzivna__,
* ako je $\omega(k)$ realna i brzina $c=\omega(k)/k$ ne ovisi o $k$, kažemo da je jednadžba __disperzivna__.
* ako je $\omega(k)$ realna i brzina $c=\omega(k)/k$ ovisi o $k$, kažemo da je jednadžba __hiperbolična__.
"""

# ╔═╡ 95c47f26-c1f3-4a7c-b964-f4200db5f945
md"""
### Advekcijsko-difuzijska jednadžba

Promotrimo __advekcijsko-difuzijsku jednadžbu__

$$
u_t+\gamma u_x=Du_{xx}. \tag{1}
$$

Uvrštavanje $u=e^{i(kx-\omega t)}$ daje 

$$
-i\omega +\gamma ki=D(ki)^2,
$$

odnosno (nakon množenja s $i$),

$$
\omega=\gamma k-iDk^2\equiv \omega(k).
$$

Dakle, relacija disperzije $\omega(k)$ je kompleksna, pa je jednadžba difuzivna. 

Rješenje jednadžbe (1) je 

$$ u(x,t)=e^{-Dk^2t}e^{ik(x-\gamma t)}.$$

Faktor $e^{ik(x-\gamma t)}$ je desni val s valnim brojem $k$, a faktor 
$\displaystyle e^{-Dk^2t}$ je opadajuća amplituda (gušenje). Zaključujemo sljedeće:

* za konstantni valni broj $k$ atenuacija (opadanje) vala se povećava s $D$ pa je $D$ mjera difuzije,
* za konstantnu vrijednost $\gamma$, atenuacija se povećava kada $k$ raste, što znači da valovi s manjim valnim duljinama opadaju brže od valova s većim valnim duljinama.

__Napomena.__ Advekcijsko-difuzijska jednadžba se još zove __Burgerova jednadžba__. Rješenje jednadžbe (1) izvest ćemo na drugi način u sljedećoj bilježnici. 
"""

# ╔═╡ 50c77c31-9d04-4d03-9e83-f93470da330b
md"""
### Poprečne vibracije grede

Promotrimo jednadžbu 

$$
u_{tt}+\gamma u_{xxxx}=0, \quad \gamma>0,
$$

koja se javlja prilikom analize poprečnih vibracija grede. 

Uvrštavanje $u=e^{i(kx-\omega t)}$ daje 

$$
(i\omega)^2+\gamma (ik)^4=0.
$$

Dakle, disperzijska relacija

$$
\omega=\pm \sqrt{\gamma}k^2
$$

je realna i brzina $c=\omega(k)/k=\pm \sqrt{\gamma}k$ nije konstantna pa je
jednadžba hiperbolična. Vidimo da vibracija s većim valnim brojem, odnosno 
manjom valnom duljinom $2\pi /k$ putuju brže.

"""

# ╔═╡ b6a15e67-d353-4065-94dc-1859298a6db1
md"""
## Linearni valovi

__Linearni valovi__ su rješenja valnih jednadžbi koje su linearne u rješenju $u$. Linearni val uvijek ima isti oblik. 

__Karakteristične krivulje__ su krivulje uzduž kojih je rješenje konstantno.

### Primjer

Rješenje problema početnih vrijednosti

\begin{align*}
&u_t+c\, u_x=0, \quad x\in\mathbb{R}, \quad t>0 \\
&u(x,0)=\phi(x), \quad x\in \mathbb{R},
\end{align*}

je 

$$
u(x,t)=\phi(x-ct).
$$

Karakteristične krivulje su pravci oblika

$$
x-ct=k, \quad k=konst.
$$

i uzduž njih je rješenje očito konstantno, odnosno početna vrijednost se propagira bez promjene. 

Alternativno, možemo provjeriti da je usmjerena derivacija uzduž karakteristične krivulje jednaka nuli: uz oznaku 

$$
x(t)=ct+k
$$

vrijedi

$$
\frac{du(x(t),t)}{dt}=u_x\frac{dx(t)}{dt}+u_t\frac{dt}{dt}
=u_x\cdot c+ u_t=0.
$$

Brzina kojom se val kreće jednaka je

$$
\frac{dx(t)}{dt}=c.
$$
"""

# ╔═╡ 214a1dfb-dde5-4ab1-92ea-d3ec99b813e9
md"""
### Primjer

Promotrimo slučaj kada je $c=c(x,t)$. Problem glasi

\begin{align*}
&u_t+c(x,t)\, u_x=0, \quad x\in\mathbb{R}, \quad t>0 \\
&u(x,0)=\phi(x), \quad x\in \mathbb{R}.
\end{align*}

Karakteristične krivulje su definirane s

$$
\frac{dx(t)}{dt}=c(x,t). \tag{1}
$$

Zaista, uzduž svake krivulje koja zadovoljava (1), usmjerena derivacije je jednaka nuli  pa je rješenje konstantno:

$$
\frac{du(x(t),t)}{dt}=u_x\frac{dx(t)}{dt}+u_t\frac{dt}{dt}
=u_x\cdot c(x,t)+ u_t=0.
$$

Brzina kojom se val kreće (kojom se početne vrijednosti propagiraju uzduž katakterističnih krivulja) dana je s (1).
"""

# ╔═╡ 6489afd0-46ab-4ecb-a132-ac16ebcd65bc
md"""
__Zadatak.__ Riješimo problem

\begin{align*}
&u_t+2t\, u_x=0, \quad x\in\mathbb{R}, \quad t>0 \\
&u(x,0)=e^{-x^2}, \quad x\in \mathbb{R}.
\end{align*}

Karakteristične krivulje su definirane jednadžbom

$$
\frac{dx(t)}{dt}=2t,
$$

čije rješenje je familija krivulja

$$
x(t)=t^2+k.
$$

Točka $(x,t)$ se nalazi na njoj pripadnoj karakterističnoj krivulji pa je vrijednost rješenje jednaka početnoj vrijednosti u točki $x=x(0)=k$. Dakle,

$$
u(x,t)=e^{-k^2}=e^{-(x-t^2)^2}.
$$

U točki $(x,t)$ val se kreće brzinom $2t$ te ubrzava s vremenom, ali uvijek ima isti oblik.
"""

# ╔═╡ b21a8fc3-d645-4ed1-9b77-7d9d8eef9ce3
# Rješenje
X=range(-2,stop=2,length=201)
T=range(0,stop=3,length=101)
u(x,t)=exp(-(x-t^2)^2)
surface(X,T,u,xlabel="x",ylabel="t",seriescolor=:blues)

# ╔═╡ d1f93b6b-186e-4b1d-9170-898969a966ab
plot(x->exp(-x^2),-5,5,label="Početni uvjet")

# ╔═╡ Cell order:
# ╟─da811f07-9dda-4191-813e-4d11aa3b7357
# ╟─7c761460-692c-4977-bdc6-9c9c0ade9ea3
# ╠═3b596ad7-5489-4591-93ed-6613711d20c4
# ╠═02409215-4359-4254-9ccc-b2c8b69bd20e
# ╠═9f51dc48-ebb2-43ee-bfc9-863e5e5e0014
# ╠═83501b08-ce5b-4f1f-9e1d-6663190c4936
# ╟─2f923b97-42c2-453d-9ff5-5f3ffa9f633f
# ╠═95fb8bf4-ada6-4247-9a1e-bdabb67570c7
# ╠═4e19a9f9-ee64-4490-9b1f-aa62067ed634
# ╠═eb656cbd-b950-4188-a829-9d00f42e37a9
# ╠═f255a3e4-0c73-42eb-a364-1fb401d04bca
# ╟─de1707e0-f4b4-4477-92ec-b3a7abc6cbaa
# ╟─95c47f26-c1f3-4a7c-b964-f4200db5f945
# ╟─50c77c31-9d04-4d03-9e83-f93470da330b
# ╟─b6a15e67-d353-4065-94dc-1859298a6db1
# ╟─214a1dfb-dde5-4ab1-92ea-d3ec99b813e9
# ╟─6489afd0-46ab-4ecb-a132-ac16ebcd65bc
# ╠═b21a8fc3-d645-4ed1-9b77-7d9d8eef9ce3
# ╠═d1f93b6b-186e-4b1d-9170-898969a966ab
