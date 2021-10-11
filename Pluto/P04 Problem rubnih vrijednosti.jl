### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 7d6aa743-e05c-4f4a-9c1e-5b29903a7386
using PlutoUI, SymPy, LinearAlgebra, Plots

# ╔═╡ 0c650d04-522a-4d3d-a3a7-eef4848452b8
plotly()

# ╔═╡ 6eff0f94-9807-43cf-97eb-a8b29d9e6f40
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 47bf3638-8134-496d-ba98-e4611aadb720
md"""
# Problem rubnih vrijednosti



# Klasifikacija

Neka je 

$$
a\cdot u_{xx}+b\cdot u_{xt}+c\cdot u_{tt}+d\cdot u_x+e\cdot u_t + f\cdot u + g = 0$$

i neka je

$$
D=b^2-4ac.$$

Vrijedi sljedeća klasifikacija:

| D               | D=0                        | D<0       | D>0            |  
| ---:            | :---                       | :---      | :--- |
| Vrsta           | parabolička                | eliptička | hiperbolička |
| Problem         | difuzija                   | ravnoteža                | valovi   |
| Domena / Metoda | omeđena / SLP              | omeđena / SLP              |        |
| Domena / Metoda | neomeđena / integr. trans. | neomeđena / integr. trans. |        |

Za neomeđeni interval $(0,\infty)$ koristi se Laplaceova transformacija, a za interval $(-\infty,\infty)$ koristi se Fourierova transformacija.
"""

# ╔═╡ 2d6cce29-1383-41e1-8e8c-653db856329b
md"""
# Jednadžba difuzije

Zadan je problem 

$$
\begin{aligned}
& u_t-u_{xx}=0  \\
& u(x,0)=|x|, \quad -2<x<2 \\
& u_x(-2,t)=0, \quad u_x(2,t)=0, \quad t>0 
\end{aligned}$$

Pretpostavimo __separaciju varijabli__ (rješenje je jedinstveno pa je svaka pretpostavka korektna ako daje rješenje):

$$
u(x,t)=X(x)T(t).$$ 

Uvrštavanje u jednadžbu daje

$$
XT'=X''T$$

odnosno (stavljamo $-\lambda$ po dogovoru)

$$
\frac{T'}{T}=\frac{X''}{X}=-\lambda,$$

za neki $\lambda \in\mathbb{R}$.

Dobili smo SLP i populacijsku jednadžbu:

1. SLP: $X''+\lambda X=0$ uz uvjete $X'(-2)=0$ i $X'(2)=0$
2. Populacijska jednadžba: $T'+\lambda T=0$

Analizirajmo tri slučaja.

**Slučaj 1**,  $\lambda=0$:

Rješenje jednadžbe $X''=0$ je $X=ax+b$ pa je $X'=a$. Rubni uvjet u lijevoj strani daje $X'(-2)=a=0$ pa je $X=b=konst.$ svojstvena funkcija, a $\lambda_0=0$ svojstvena vrijednost.

**Slučaj 2**, $\lambda<0$:

Rješenje jednadžbe $X''=-\lambda X$ je $X=a e^{\sqrt{-\lambda}x}+ b e^{-\sqrt{-\lambda}x}$. Vrijedi 

$$X'=a\sqrt{-\lambda} e^{\sqrt{-\lambda}x}- b\sqrt{-\lambda} e^{-\sqrt{-\lambda}x}.$$

Uvrštavanje rubnih uvjeta daje

$$
\begin{aligned}
X'(-2)& =\sqrt{-\lambda} \big(a e^{-2\sqrt{-\lambda}}-be^{2\sqrt{-\lambda}}\big)=0\\
X'(2) & = \sqrt{-\lambda} \big(a e^{2\sqrt{-\lambda}}-be^{-2\sqrt{-\lambda}}\big)=0.
\end{aligned}$$

Iz prve jednadžbe slijedi $a=be^{4\sqrt{-\lambda}}$ pa uvrštavanje u drugu jednadžbu daje $b\big(e^{6\sqrt{-\lambda}}-e^{-2\sqrt{-\lambda}}\big)=0$. Izraz u zagradi je nula samo za $\lambda=0$ pa je $b=0$. Iz prve jednadžbe onda slijedi $a=0$. Dakle, $X=0$ ne može biti svojstvena funkcija pa $\lambda<0$ nije svojstvena vrijednost.

**Slučaj 3**, $\lambda>0$:

Rješenje jednadžbe $X''=-\lambda X$ je $X=a \sin \sqrt{\lambda}x+ b \cos\sqrt{\lambda}x$. Početni uvjet $u(x,0)=|x|$ je parna funkcija pa možemo uzati koeficijent uz $\sin(\sqrt{\lambda}x)$ jednak nuli. Dakle, $X=b \cos\sqrt{\lambda}x$ pa je $X'=-b\sqrt{\lambda}\sin\sqrt{\lambda}x$. Rubni uvjeti daju

$$
\begin{aligned}
X'(-2)& =-b\sqrt{\lambda}\sin(-2\sqrt{\lambda})=b\sqrt{\lambda}\sin(2\sqrt{\lambda})=0 \\
X'(2) & = b\sqrt{\lambda}\sin(2\sqrt{\lambda})=0.
\end{aligned}$$

Kako je $b\neq 0$ i $\sqrt{\lambda}\neq 0$, vrijedi $2\sqrt{\lambda}=n\pi$, $n\in\mathbb{N}$. 

Dakle, za $\lambda \geq 0$ SLP ima beskonačno svojstvenih vrijednosti

$$
\lambda_n = \frac{n^2\pi^2}{4}, \quad n\in \mathbb{N}\cup \{0\},$$

i pripadne svojstvene funkcije 

$$
X_n(x)= A_n \cos \big(\frac{n\pi}{2}x\big).$$

Za svaki $\lambda_n$ rješenje populacijske jednadžbe glasi 

$$
T_n(t)=B_n e^{\displaystyle-\frac{n^2\pi^2}{4} t}$$

što zajedno daje 

$$
u_n(x,t)= C_n \cos \big(\frac{n\pi}{2}x\big ) e^{\displaystyle-\frac{n^2\pi^2}{4} t}.$$

Svaka funkcija $u_n$ zadovoljava jednadžbu i rubne uvjete. 

Prema __principu superpozicije__ i funkcija 

$$
u(x,t)=\sum_{n=0}^\infty u_n(x,t)= \sum_{n=0}^\infty C_n \cos \big(\frac{n\pi}{2}x\big ) e^{\displaystyle-\frac{n^2\pi^2}{4} t}$$

također zadovoljava jednadžbu i rubne uvjete. Treba još odabrati koeficijente $C_n$ tako da se zadovolji i početni uvijet - radi se o razvoju u __(generalizirani) Fourierov red__ :

$$
\begin{aligned}
& u(x,0)=\sum_{n=0}^\infty C_n \cos \big(\frac{n\pi}{2}x\big )=|x| \\
& C_n=  \displaystyle \frac{ \big(|x|, \cos \big(\frac{n\pi}{2}x\big ) \big)}
{\big (\cos \big(\frac{n\pi}{2}x\big ), \cos \big(\frac{n\pi}{2}x\big ) \big)}.
\end{aligned}$$
"""

# ╔═╡ 79e6ad88-8838-4029-bead-99de981da955
md"""
Probajmo simboličko računanje - treba nam paket `SymPy.jl`:
"""

# ╔═╡ fc7e565d-680b-476c-bf89-4091e9c7bc0b
begin
	# Definirajmo simbole
	n=symbols("n",integer=True,nonnegative=True)
	x=symbols("x",real=True)
end

# ╔═╡ 671eea4b-486a-4e77-98bb-1dca4c56873c
begin
	# Definirajmo skalarni produkt
	import LinearAlgebra.⋅
	⋅(f,g,a,b)=integrate(f*g,(x,a,b))
end

# ╔═╡ c4d89894-088f-446c-b4ce-d9c78af73d2c
# Umjesto 
g=abs(x)

# ╔═╡ 9b90ab1d-1843-4a3b-813e-eac93903bb05
md"""
Pogledajmo zadani početni uvjet:
"""

# ╔═╡ 407e4a6d-bb34-40ea-8694-1e4a21a15193
plot(g(x),-2,2)

# ╔═╡ b08a5be1-90b7-4c36-b4ac-83ca75b3d31f
f(n,x)=cos(n*PI*x/2)

# ╔═╡ 49de3080-a15b-4ea6-bfcb-d9fd4c5063fa
# Na primjer
f(2,x)

# ╔═╡ c66069da-f439-4f7a-b1a9-8970051b221f
f(0,x)

# ╔═╡ c8c5a472-c781-4894-b557-8ff121ba4077
md"""
Izračunajmo koeficijente $C_n$:
"""

# ╔═╡ 07016beb-48d0-41bd-9798-ba0cc9cf10fd
C(n)=⋅(g(x),f(n,x),-2,2)/⋅(f(n,x),f(n,x),-2,2)

# ╔═╡ 8a64a59f-b6c7-42e8-ace2-e47b27fc9d0e
C(0)

# ╔═╡ d70c37fa-6bfc-43a1-a3d1-e53d680bbd78
C(1)

# ╔═╡ 3d20eccd-637e-46f9-aa4e-f6f730e445ee
C(2)

# ╔═╡ 051777e6-45df-47df-81f7-46d07d679133
C(3)

# ╔═╡ f466ab29-7dc0-4435-8feb-f40549367b91
C(4)

# ╔═╡ 694f43db-af83-46aa-ba05-ae2c40e2f1b8
C(5)

# ╔═╡ b2b2fcef-69d6-4a6c-860d-92b7ab2cafca
# Opća formula
⋅(g(x),f(n,x),-2,2)/⋅(f(n,x),f(n,x),-2,2)

# ╔═╡ 9327f715-8ef6-4187-8e0f-674caba4ae41
md"""
Vidimo da je 

$$
\begin{aligned}
& C_0=1, \\
& C_{2k}=0, \\
& C_{2k-1}=\displaystyle\frac{-8}{(2k-1)^2\pi^2},
\end{aligned}$$

odnosno

$$
u(x,t)=1 - \sum_{k=1}^\infty \frac{8}{(2k-1)^2\pi^2}
\cos \bigg(\frac{(2k-1)\pi}{2}x\bigg ) e^{\displaystyle-\frac{(2k-1)^2\pi^2}{4} t}.$$
"""

# ╔═╡ 4c9eba19-f68e-4865-b8d2-e5b620212d4e
md"""
Definirajmo sumu prvih $n$ članova reda:
"""

# ╔═╡ f5edf223-60a7-48f1-8386-8dabaaba5790
begin
	k=symbols("k", integer=True, nonnegative=True)
	t=symbols("t", real=True, nonnegative=True)
end

# ╔═╡ 79b3caec-b5f5-4ba4-8550-cc202eb5d9d0
u(n,x,t)=C(n)*f(n,x)*exp(-(n^2*PI^2/4)*t)

# ╔═╡ 9a95ed5f-94fd-487a-bfd2-fbea45124db3
# Na primjer
u(0,x,t)

# ╔═╡ f96fdfc8-6bf7-4a67-b18e-0771efe377eb
u(3,x,t)

# ╔═╡ 1f7ef259-df1e-4e1b-8bff-5f03ae2252ed
# u(3,x,t) u zadanoj točki
u(3,0.5,0.5)

# ╔═╡ 05419561-4764-418c-8184-c6be1ea6bfcc
# Numerička vrijednost (BigFloat)
N(u(3,0.5,0.5))

# ╔═╡ 6f71559b-356a-414d-b226-a70fd50db330
# Suma prvih n članova reda
U(n,x,t)=sum([u(k,x,t) for k=0:n])

# ╔═╡ cb246811-3d02-461e-8700-47dc8f8f347a
# Na primjer
U(5,0.5,0.5)

# ╔═╡ 42f86e27-6812-4c3c-979d-860c40728a67
# Numerička vrijednost
N(U(5,0.5,0.5))

# ╔═╡ 77b2e2ef-3456-4e53-b71f-0a49d2e0478c
# Za t=0 ovo mora konvergirati u |x|
@time N(U(11,0.5,0.0))

# ╔═╡ a6b00d81-41b9-4880-961a-5a63f02cd9bd
md"""
__Napomena:__ Radi se o simboličkom računanju pa ne treba pretjerivati s $n$.

## Crtanje
"""

# ╔═╡ d9d895a8-7b12-4415-8902-805ce2471f1c
begin
	m=17
	X=range(-2,2, length=m)
	T=range(0,5,length=m)
end

# ╔═╡ 73b2542c-1e8a-48ce-badb-e46baf35a22d
X

# ╔═╡ 5b99b88d-095f-4b43-b103-1604e0a2aa43
collect(X)

# ╔═╡ 29280217-8af2-4679-a9a8-26b83c65b7c2
# Radi brzine pripremimo U(9) unaprijed
U9=U(9,x,t)

# ╔═╡ 4b8df6fe-6a9b-4cd1-b25b-fc8c74759d64
# Sada je puno brže jer se samo uvrštava.
@time N(U9(0,2.0))

# ╔═╡ 3ed3b591-f6fb-4294-b9a2-0fe2dc6ce1f2
begin
	# Ovo je sada brzo (3-4 sek)
	z = zeros(m,m)
	for i in 1:m
	    for j in 1:m
	        z[i,j] = Float64(U9(X[i],T[j]))
	    end
	end
end

# ╔═╡ ba8fcfb5-029a-4b40-a0b4-f4c9d0c9c067
z

# ╔═╡ 7888a4ef-6761-40ce-b5b6-bae3381d8264
surface(X,T,z',xlabel="x",ylabel="t")

# ╔═╡ 510c9f4c-bcf1-4a19-ae12-2211d01bf192
md"""
# Numeričko računanje i crtanje
"""

# ╔═╡ 52f2f326-8fb0-465b-93be-e019afe93061
begin
	Xₙ=range(-2,2,length=51)
	Tₙ=range(0,5,length=51)
	XT=collect(Iterators.product(Xₙ,Tₙ))
end

# ╔═╡ 4b811a75-9b72-45ed-bbaa-8e323c588a36
begin
	# Probajmo l od 1 do 10
	l=5
	h(xt)=1-8*sum([cos.((2*k-1)*pi*xt[1]/2).*exp.(-(2*k-1)^2*pi^2*xt[2]/4)/((2*k-1)^2*π^2) 
	        for k=1:l])
	surface(Xₙ,Tₙ,map(h,XT)',xlabel="x",ylabel="t")
end

# ╔═╡ de02f5f9-9559-445c-b883-96aad3292dc8
md"""
## Primjer

$$
\begin{aligned}
& u_t-u_{xx}=-u \\
& u(x,0)=f(x)=\begin{cases}0, \quad -1<x<0 \\ x,\quad 0<x<1 \end{cases} \\
& u(-1,t)=0,\quad u(1,t)=0 
\end{aligned}$$

Za detalje o simboličkom računanju pogledajte [SymPy Tutorial](https://github.com/jverzani/SymPy.jl/blob/master/examples/tutorial.md).

Uvrštavanjem

$$
u(x,t)=X(x)T(t)$$

jednadžba prelazi u jednadžbu

$$
T'X-TX''=-TX,$$

što daje dvije jednadžbe:

$$
\frac{X''}{X}=\frac{T'+T}{T}=-\lambda.$$

Jednadžba po $T$ je populacijska jednadžba koja glasi

$$
T'=-(\lambda+1)T$$

i čije rješenje je

$$
T=Ce^{-(\lambda+1)t}.$$

Riješimo SLP po $X$:

$$
X''=-\lambda X, \quad X(-1)=0, \quad X(1)=0.$$
"""

# ╔═╡ 80cb356e-065c-43ac-85a5-0d9c64cd9176
F = SymFunction("F")

# ╔═╡ e3c68a65-82ef-492e-9283-88edd8cc26e2
begin
	l₁=symbols("l",real=true,positive=true)
	diffeq = Eq(diff(F(x), x, x) +l₁*F(x), 0)
end

# ╔═╡ ab5389ea-0126-4974-8e57-7d2782ca8b76
ex = dsolve(diffeq)

# ╔═╡ 3e56b49d-9462-4426-a88d-d6a45ce2bf9c
ex1 = rhs(ex)

# ╔═╡ 79794948-8b69-480d-a34c-a40ab0b13300
md"""
Uvrstimo rubne uvjete:
"""

# ╔═╡ 75b46887-f4c5-4514-aec9-e94a54477ba1
ex1a=subs(ex1,x,-1)

# ╔═╡ e6abef29-40bb-4a09-9163-39947e5364e3
ex1b=subs(ex1,x,1)

# ╔═╡ cf271e99-4be3-4fac-bcde-0988323000ce
solve(cos(sqrt(l₁)),l₁)

# ╔═╡ 0a6d295e-68be-41a9-9312-1c8fd8549cd3
md"""
Sustav jednadžbi je homogen i glasi

$$
\begin{bmatrix} -C_1 & C_2 \\ C_1 & C_2 \end{bmatrix} \begin{bmatrix}\sin \sqrt{\lambda} \\ \cos\sqrt{\lambda} \end{bmatrix} = \begin{bmatrix} 0\\ 0\end{bmatrix}.$$

Trivijalno rješenje je u ovom slučaju očito nemoguće, a netrivijalna rješenje postoje kada je matrica sustava singularna, odnosno kada je $C_1=0$ ili $C_2=0$.

Kada je $C_1=0$ onda je $\cos\sqrt{\lambda}=0$ pa je 

$$
\sqrt{\lambda}=\frac{2n+1}{2}\pi, \quad n=0,1,2,3,\ldots$$

Kada je $C_2=0$ onda je $\sin\sqrt{\lambda}=0$ pa je 

$$
\sqrt{\lambda}=n\pi, \quad n=0,1,2,3,\ldots$$

Dakle, rješenje problema koje zadovoljava jednadžbu i rubne uvjete ima oblik:

$$
u(x,t)=\sum_{n=0}^\infty a_n \cos \bigg(\frac{2n+1}{2}\pi x\bigg)
e^{-\big(\big[\frac{2n+1}{2}\pi\big]^2+1\big)t}+b_n \sin (n\pi x)\,e^{-([n\pi]^2+1)t}.$$

Potrebno je zadovoljiti još početni uvjet:

$$
u(x,0)=\sum_{n=0}^\infty a_n \cos \bigg(\frac{2n+1}{2}\pi x\bigg)+b_n \sin (n\pi x)=f(x).$$
 
Radi se o razvoju u generalizirani Fourierov red funkcije f(x): 
"""

# ╔═╡ 788235d8-c271-4d8f-abd6-9e3316c01298
# p=piecewise((0,Lt(x,0)),(x,Ge(x,0)))
p(x)=x*Heaviside(x)

# ╔═╡ 8610ef28-dd89-4cb1-9717-b237a20eeede
md"""
Provjerimo ortonormiranost sustava funkcija.
"""

# ╔═╡ e02874e6-d676-420b-9e6a-aa7cfd5b6ce2
⋅(cos((2*n+1)*PI*x/2),cos((2*n+1)*PI*x/2),-1,1)

# ╔═╡ e0a56710-d29a-4e43-ac26-afb0d64056f0
⋅(sin(n*PI*x),sin(n*PI*x),-1,1)

# ╔═╡ e6fb6e70-5621-4f12-876d-a7fc1c083f91
md"""
Norme svih funkcija su jednake $1$ pa ne trebamo računati nazivnike.
"""

# ╔═╡ 329fc0df-884a-4842-b27f-c6115e42c0bc
a(n)=⋅(p(x),cos((2*n+1)*PI*x/2),-1,1)

# ╔═╡ 0c2f0034-26a2-4d03-bd70-478ebe5110c3
a(0)

# ╔═╡ fd086bd3-1a9e-4f10-9707-f9a175ae4518
N(a(0))

# ╔═╡ 5d0b9b19-e55b-46e2-a12c-83c730828a83
b(n)=⋅(p(x),sin(n*PI*x),-1,1)

# ╔═╡ d02415b3-4e57-44a6-8be3-09c7abbdbcc4
b(0)

# ╔═╡ bcd1fe66-6389-463e-aaa0-a9764338b3d9
b(17)

# ╔═╡ 60b3aaeb-cb9c-4a6d-837e-bd2e20737292
# Opća formule za a(n)
⋅(p(x),cos((2*n+1)*PI*x/2),-1,1)

# ╔═╡ 24e26b9e-0894-4e53-9cb3-fb1b79b2ec6f
# Opća formule za b(n)
⋅(p(x),sin(n*PI*x),-1,1)

# ╔═╡ 816ebb55-4808-4829-b232-5b7f670ae0c7
md"""
Pripremimo se za brže računanje tako da ćemo unaprijed izračunati numeričke
vrijednosti koeficijenata $a_n$ i $b_n$. 
"""

# ╔═╡ 67d96ec1-688c-4fbd-a046-989c4001c089
A=[Float64(a(n)) for n=0:20]

# ╔═╡ aa81d6a7-9103-4d13-8071-ec1355d2d1a2
B=[Float64(b(n)) for n=0:20]

# ╔═╡ bab59fda-bcc6-4094-b3d4-a471b4c305ee
begin
	X₁=range(-1,1,length=51)
	T₁=range(0,5,length=51)
	XT₁=collect(Iterators.product(X₁,T₁))
end;

# ╔═╡ a21d0a52-ba97-4d31-be17-c0ce97cee4fe
begin
	lₙ=20
	h₁(xt)=sum([A[k]*cos.((2*k-1)*pi*xt[1]/2).*exp.(-(((2*k-1)*pi/2)^2/4+1)*xt[2])+
	        B[k]*sin.((k-1)*pi*xt[1]).*exp.(-(((k-1)*pi)^2+1)*xt[2]) for k=1:lₙ]) 
	surface(X₁,T₁,map(h₁,XT₁)',xlabel="x",ylabel="t")
end

# ╔═╡ a336d3c7-66fe-4ce4-b318-3736b705fd13
md"""
# Homogenizacija

U oba prethodna primjera zadani su homogeni rubni uvjeti. Ukoliko rubni uvjeti nisu homogeni, zadani problem je potrebno __homogenizirati__ kako bi mogli dobiti regularni SLP.

Navedimo primjer. Neka je zadan problem

$$
\begin{aligned}
& u_t -u_{xx}=0,\quad 0<x<l,\quad t>0 \\
& u(x,0)=f(x),\quad 0<x<l \\
& u(0,t)=g(t),\quad u(l,t)=h(t),\quad t>0.
\end{aligned}$$
 
Nađimo rješenje u obliku

$$
u(x,t)=v(x,t)+U(x,t),$$

gdje je $v$ rješenje problema sa homogenim rubnim uvjetima. Vrijedi

$$
\begin{aligned}
& u=v+U\\
& u_t=v_t+U_t\\
& u_{xx}=v_{xx}+U_{xx}
\end{aligned}$$

pa zadana PDJ prelazi u 

$$
v_t+U_t=v_{xx}+U_{xx}.$$

Početni uvjet za $v$ glasi

$$
v(x,0)=u(x,0)-U(x,0)=f(x)-U(x,0),$$

a rubni uvjeti glase

$$
\begin{aligned}
& v(0,t)=u(0,t)-U(0,t)=g(t)-U(0,t)=0\quad  \textrm{(želimo homogeni uvjet)}\\
& v(l,t)=u(l,t)-U(l,t)=h(t)-U(l,t)=0 \quad  \textrm{(želimo homogeni uvjet)}
\end{aligned}$$

Zaključujemo da će $v$ zadovoljavati homogene rubne uvjete ako je 

$$
U(x,t)=g(t)+\displaystyle\frac{x}{l}[h(t)-g(t)],\quad 0<x<l.$$

Za ovako definiranu funkciju $U$ vrijedi

$$
\begin{aligned}
& U_t=g'(t)+\displaystyle\frac{x}{l}[h'(t)-g'(t)]\\
& U_{xx}=0.
\end{aligned}$$

Uvrštavanjem slijedi da je $v$ rješenje __homogenog__ reakcijsko-difuzijskog problema

$$
\begin{aligned}
&v_t=v_{xx}-g'(t)-\displaystyle\frac{x}{l}[h'(t)-g'(t)], \quad 0<x<l,\quad t>0
\\
&v(x,0)=f(x)-g(0)-\displaystyle\frac{x}{l}[h(0)-g(0)], \quad 0<x<l
\\
& v(0,t)=0,\quad v(l,t)=0,\quad t>0,
\end{aligned}$$

dok je rješenje polaznog problema

$$
u(x,t)=v(x,t)+g(t)+\displaystyle\frac{x}{l}[h(t)-g(t)].$$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
Plots = "~1.22.4"
PlutoUI = "~0.7.16"
SymPy = "~1.0.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "a325370b9dd0e6bf5656a6f1a7ae80755f8ccc46"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.7.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

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
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cafe0823979a5c9bff86224b3b8de29ea5a44b2e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.61.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "f6532909bf3d40b308a0f360b6a0e626c0e263a8"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.1"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "34dc30f868e368f8a17b728a1238f3fcda43931a"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.3"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "6841db754bd01a91d281370d9a0f8787e220ae08"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.4"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

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
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

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

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "793793f1df98e3d7d554b65a107e9c9a6399a6ed"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.7.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "65fb73045d0e9aaa39ea9a29a5e7506d9ef6511f"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.11"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SymPy]]
deps = ["CommonEq", "CommonSolve", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "1ef257ecbcab8058595a68ca36a6844b41babcbd"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.0.52"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═7d6aa743-e05c-4f4a-9c1e-5b29903a7386
# ╠═0c650d04-522a-4d3d-a3a7-eef4848452b8
# ╠═6eff0f94-9807-43cf-97eb-a8b29d9e6f40
# ╟─47bf3638-8134-496d-ba98-e4611aadb720
# ╟─2d6cce29-1383-41e1-8e8c-653db856329b
# ╟─79e6ad88-8838-4029-bead-99de981da955
# ╠═fc7e565d-680b-476c-bf89-4091e9c7bc0b
# ╠═671eea4b-486a-4e77-98bb-1dca4c56873c
# ╠═c4d89894-088f-446c-b4ce-d9c78af73d2c
# ╟─9b90ab1d-1843-4a3b-813e-eac93903bb05
# ╠═407e4a6d-bb34-40ea-8694-1e4a21a15193
# ╠═b08a5be1-90b7-4c36-b4ac-83ca75b3d31f
# ╠═49de3080-a15b-4ea6-bfcb-d9fd4c5063fa
# ╠═c66069da-f439-4f7a-b1a9-8970051b221f
# ╟─c8c5a472-c781-4894-b557-8ff121ba4077
# ╠═07016beb-48d0-41bd-9798-ba0cc9cf10fd
# ╠═8a64a59f-b6c7-42e8-ace2-e47b27fc9d0e
# ╠═d70c37fa-6bfc-43a1-a3d1-e53d680bbd78
# ╠═3d20eccd-637e-46f9-aa4e-f6f730e445ee
# ╠═051777e6-45df-47df-81f7-46d07d679133
# ╠═f466ab29-7dc0-4435-8feb-f40549367b91
# ╠═694f43db-af83-46aa-ba05-ae2c40e2f1b8
# ╠═b2b2fcef-69d6-4a6c-860d-92b7ab2cafca
# ╟─9327f715-8ef6-4187-8e0f-674caba4ae41
# ╟─4c9eba19-f68e-4865-b8d2-e5b620212d4e
# ╠═f5edf223-60a7-48f1-8386-8dabaaba5790
# ╠═79b3caec-b5f5-4ba4-8550-cc202eb5d9d0
# ╠═9a95ed5f-94fd-487a-bfd2-fbea45124db3
# ╠═f96fdfc8-6bf7-4a67-b18e-0771efe377eb
# ╠═1f7ef259-df1e-4e1b-8bff-5f03ae2252ed
# ╠═05419561-4764-418c-8184-c6be1ea6bfcc
# ╠═6f71559b-356a-414d-b226-a70fd50db330
# ╠═cb246811-3d02-461e-8700-47dc8f8f347a
# ╠═42f86e27-6812-4c3c-979d-860c40728a67
# ╠═77b2e2ef-3456-4e53-b71f-0a49d2e0478c
# ╟─a6b00d81-41b9-4880-961a-5a63f02cd9bd
# ╠═d9d895a8-7b12-4415-8902-805ce2471f1c
# ╠═73b2542c-1e8a-48ce-badb-e46baf35a22d
# ╠═5b99b88d-095f-4b43-b103-1604e0a2aa43
# ╠═29280217-8af2-4679-a9a8-26b83c65b7c2
# ╠═4b8df6fe-6a9b-4cd1-b25b-fc8c74759d64
# ╠═3ed3b591-f6fb-4294-b9a2-0fe2dc6ce1f2
# ╠═ba8fcfb5-029a-4b40-a0b4-f4c9d0c9c067
# ╠═7888a4ef-6761-40ce-b5b6-bae3381d8264
# ╠═510c9f4c-bcf1-4a19-ae12-2211d01bf192
# ╠═52f2f326-8fb0-465b-93be-e019afe93061
# ╠═4b811a75-9b72-45ed-bbaa-8e323c588a36
# ╟─de02f5f9-9559-445c-b883-96aad3292dc8
# ╠═80cb356e-065c-43ac-85a5-0d9c64cd9176
# ╠═e3c68a65-82ef-492e-9283-88edd8cc26e2
# ╠═ab5389ea-0126-4974-8e57-7d2782ca8b76
# ╠═3e56b49d-9462-4426-a88d-d6a45ce2bf9c
# ╟─79794948-8b69-480d-a34c-a40ab0b13300
# ╠═75b46887-f4c5-4514-aec9-e94a54477ba1
# ╠═e6abef29-40bb-4a09-9163-39947e5364e3
# ╠═cf271e99-4be3-4fac-bcde-0988323000ce
# ╟─0a6d295e-68be-41a9-9312-1c8fd8549cd3
# ╠═788235d8-c271-4d8f-abd6-9e3316c01298
# ╟─8610ef28-dd89-4cb1-9717-b237a20eeede
# ╠═e02874e6-d676-420b-9e6a-aa7cfd5b6ce2
# ╠═e0a56710-d29a-4e43-ac26-afb0d64056f0
# ╟─e6fb6e70-5621-4f12-876d-a7fc1c083f91
# ╠═329fc0df-884a-4842-b27f-c6115e42c0bc
# ╠═0c2f0034-26a2-4d03-bd70-478ebe5110c3
# ╠═fd086bd3-1a9e-4f10-9707-f9a175ae4518
# ╠═5d0b9b19-e55b-46e2-a12c-83c730828a83
# ╠═d02415b3-4e57-44a6-8be3-09c7abbdbcc4
# ╠═bcd1fe66-6389-463e-aaa0-a9764338b3d9
# ╠═60b3aaeb-cb9c-4a6d-837e-bd2e20737292
# ╠═24e26b9e-0894-4e53-9cb3-fb1b79b2ec6f
# ╟─816ebb55-4808-4829-b232-5b7f670ae0c7
# ╠═67d96ec1-688c-4fbd-a046-989c4001c089
# ╠═aa81d6a7-9103-4d13-8071-ec1355d2d1a2
# ╠═bab59fda-bcc6-4094-b3d4-a471b4c305ee
# ╠═a21d0a52-ba97-4d31-be17-c0ce97cee4fe
# ╟─a336d3c7-66fe-4ce4-b318-3736b705fd13
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002