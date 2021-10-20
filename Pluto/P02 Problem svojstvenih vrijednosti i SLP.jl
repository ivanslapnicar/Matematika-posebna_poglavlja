### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ 55e86888-4b1b-4f6d-9bd8-c839e76fe381
using PlutoUI, Random, LinearAlgebra, SymPy, Plots

# ╔═╡ c5e5c059-ba8d-4d1d-b171-9c7673761320
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 68499249-33bd-4c9b-9ba4-072e096c4820
md"""
# Problem svojstvenih vrijednosti i SLP


# Matrični problem svojstvenih vrijednosti

Neka je $A\in\mathbb{R}^{n\times n}$ kvadratna realna matrica. 

Tražimo __svojstvene vrijednosti__  $\lambda\in\mathbb{R}$ i __svojstvene vektore__ $x\in\mathbb{R}^{n}\neq 0$, takve da je 

$$Ax=\lambda x.$$

Dakle, $A$ djeluje na vektor $x$ tako da ga produži ili skrati, eventualno promijeni orijentaciju, dok smjer ostaje isti. 

Vrijedi 

$$Ax-\lambda I x=(A-\lambda I)x=0.$$

Ovo je homogeni sustav linearnih jednadžbi koji ima netrivijalna rješenja ($x\neq 0$) ako i samo ako je matrica sustava $(A-\lambda I)$ singularna, odnosno ako i samo ako je 

$$
\det(A-\lambda I)=0.$$

Izraz $\det(A-\lambda I)$ je polinom stupnja $n$ u varijabli $\lambda$ s realnim koeficijentima, koji, prema osnovnom teoremu algebre, ima 
$n$ nul-točaka koje su ili realne ili dolaze u konjugirano kompleksnim parovima. 

__Teorem__. Svojstveni vektori koji odgovaraju različitim svojstvenim vrijednostima su linearno nezavisni.

_Dokaz_ : Neka je

$$Ax=\lambda x, \quad Ay=\mu y, \quad x,y\neq 0, \quad \lambda\neq\mu.$$

Pretpostavimo da su $x$ i $y$ linearno zavisni, odnosno, 

$$
\alpha x+\beta y=0, \quad |\alpha|+|\beta|>0.$$

Vrijedi

$$
A\cdot (\alpha x+\beta y)=\alpha\lambda x+\beta\mu y=A\cdot 0=0.$$

Množenje prve jednakosti s $\lambda$ daje sustav

$$
\begin{aligned}
\lambda\alpha x+\lambda\beta y&=0 \\
\alpha\lambda x+\beta\mu y&=0.
\end{aligned}$$

Oduzmanje prve jednadžbu od druge daje

$$
\beta(\mu-\lambda)y=0$$

pa je, zbog $\mu-\lambda\neq 0$ i $y\neq 0$, nužno $\beta=0$. 
Uvrštavanjem u originalnu lineranu kombinaciju, zbog $x\neq 0$ 
slijedi $\alpha=0$ pa su $x$ i $y$ linearno nezavisni.

__Teorem__.  Ako je $A$ simetrična matrica, $A=A^T$, tada su sve svojstvene vrijednosti realne i imaju ortogonalni skup svojstvenih vektora, osnosno postoji matrica $U$ takva da je 

$$U^TU=UU^T=I,\quad 
A=U\Lambda U^T, \quad AU=U\Lambda, \quad A=\sum_{i=1}^n \lambda_i u_i u_i^T.$$
"""

# ╔═╡ 5bad68c8-ccbf-4fdd-b22b-8fa78a3b95bb
md"""
## Geometrijska interpretacija_

Zadana je matrica $A=\begin{pmatrix} 1 & 2\cr 2 & 3 \end{pmatrix}$. Vrijedi

$$
\det (A-\lambda I) = \det \begin{pmatrix} 1-\lambda & 2\cr 2 & 3-\lambda \end{pmatrix}
=\lambda^2-4\lambda-1.$$

Nul-točke su 

$$
\lambda_1=2-\sqrt{5}, \qquad \lambda_2=2+\sqrt{5}.$$

"""

# ╔═╡ ab8033eb-60b2-4cbf-9dae-15c6753bba76
A=[1 2;2 3]

# ╔═╡ 21d8db39-7058-49d7-9174-f1b88ff7c062
eigvals(A)

# ╔═╡ d6a6f2b1-d64f-4d7b-94d2-f98fa36c0996
md"""
Engleski pojam __eigenvalue__ je prijevod njemačkog pojma __Eigenwert__.
Riječca _Eigen_ se u njemačkom jeziku nalazi u riječima _Eigentum_ - _vlasništvo_ i 
_Eigenshaft_ - _svojstvo_ pa se u hrvatskoj literaturi koriste pojmovi
__svojstvena vrijednost__ ili __vlastita vrijednost__.

U literaturi se još može naići na termin __principal values__ , odnosno __glavne vrijednosti__. 

"""

# ╔═╡ b2deff23-48d4-4495-bb73-35dfa0f55487
2-sqrt(5),2+sqrt(5)

# ╔═╡ a46ec17d-6c24-455a-b049-e521f79b0667
begin
	λ₁=2-sqrt(5)
	λ₂=2+sqrt(5)
end

# ╔═╡ ac54aa9d-b480-434d-83fc-69f2f71044bd
md"""
Svojstveni vektori su netrivijalna rješenja homogenog sustava

$$
\begin{pmatrix} 1 & 2\cr 2 & 3 \end{pmatrix}\begin{pmatrix} x \cr y\end{pmatrix}
=\lambda \begin{pmatrix} x \cr y\end{pmatrix}.$$

U ovom slučaju možemo staviti $x=1$ pa je $1+2y=\lambda$, odnosno $y=(\lambda-1)/2$.
Dakle, svojstveni vektori su

$$
x_1=\begin{pmatrix} 1 \cr \displaystyle\frac{1-\sqrt{5}}{2} \end{pmatrix},
\quad
x_2=\begin{pmatrix} 1 \cr \displaystyle\frac{1+\sqrt{5}}{2} \end{pmatrix}.$$

"""

# ╔═╡ 52597f5d-faed-438f-ae5a-554b63a88e78
begin
	x₁=[1; (1-sqrt(5))/2]
	x₂=[1;(1+sqrt(5))/2]
end

# ╔═╡ 3f53d041-f8b1-4483-92de-93e918a3c2e3
A*x₁-λ₁*x₁

# ╔═╡ 1f90b141-4403-477d-b57d-56493559b275
A*x₂-λ₂*x₂

# ╔═╡ b1309118-5f38-412b-8e51-832feacfb0e2
md"""
Funkcija `eigen()` računa i svojstvene vrijednosti i svojstvene vektore.
Odgovarajuća Matlab naredba je `eig`.
"""

# ╔═╡ 7156ab8c-5ec6-46d0-ba7f-65cc231f2fe7
B=eigen(A)

# ╔═╡ 7f0b797f-57a0-4287-abf7-028c357e7da2
B.vectors

# ╔═╡ 299ff3f2-24bb-46d8-bbea-c2e050282779
x₁/norm(x₁)

# ╔═╡ 52faa36a-b161-4d4e-a474-c49ec7f9f245
B.vectors'*B.vectors

# ╔═╡ cdbfa73e-9783-4256-9035-4afe0ae54cbb
B.vectors*B.vectors'

# ╔═╡ 1ab3bd33-27b6-4f95-ae41-6e2689b6ddc0
A*B.vectors-B.vectors*Diagonal(B.values)

# ╔═╡ ab619c12-f92f-42d2-877a-407e17c3d671
λ₁*B.vectors[:,1]*B.vectors[:,1]'+λ₂*B.vectors[:,2]*B.vectors[:,2]'

# ╔═╡ 79f3c8f8-30e6-423e-8a2e-a402f9f8b972
# Kutevi od 0 do 2Π
ϕ=range(0,stop=2*pi,length=200)

# ╔═╡ 5df96e85-43f3-46e2-b1b0-739637782191
# Stupci od X su točke na jediničnoj kružnici
X=[sin.(ϕ)'; cos.(ϕ)']

# ╔═╡ 73f65db7-9958-4044-990f-db91e0c4ea1e
scatter(X[1,:],X[2,:],aspect_ratio=1)

# ╔═╡ d2038ec5-8957-4c9f-b3eb-89ed0a9fbc16
# Stupci od Y su točke u koje se preslikavaju točke na jediničnoj kružnici
Y=A*X

# ╔═╡ 3c3bf03a-1473-45a7-9aae-f4b33f6268b0
scatter!(Y[1,:],Y[2,:])

# ╔═╡ 4789f475-d239-4488-93d4-0989e2489835
md"""
Jedinična kružnica se preslikava u elipsu. 

Duljina manje poluosi elipse je $|\lambda_1|$, duljina veće poluosi je $|\lambda_2|$.

Smjerovi poluosiju su pripadni svojstveni vektori.
"""

# ╔═╡ 8b032f4c-7ef2-4582-9e1e-ec9ec4e45d80
Random.seed!(123)

# ╔═╡ 80a7e0ba-87a6-48ae-84b1-9aeb354f9a6a
md"
## Primjer
"

# ╔═╡ ed37fc89-01a7-497d-a4ca-ba584ef21781
A₁=Symmetric(rand(-8:8,6,6))

# ╔═╡ 12964840-d516-423d-b0fd-17baea64fc58
B₁=eigen(A₁);

# ╔═╡ d4aba6a0-3012-4e40-8495-9957dce33a1c
B₁.values

# ╔═╡ a9160ec4-50ff-4067-bfde-61e2234af4d6
B₁.vectors

# ╔═╡ 27423aed-47a5-473c-8016-98a39b8fd7e2
# Ortogonalnost matrice svojstvenih vektora
norm(B₁.vectors'*B₁.vectors-I)

# ╔═╡ 89fe6b2f-6106-4695-92fb-62e57d1bae0b
norm(B₁.vectors*B₁.vectors'-I)

# ╔═╡ d477b32b-f85c-4f3c-9f96-8726725a7b3a
# Provjerimo točnost rastava
B₁.vectors*Diagonal(B₁.values)*B₁.vectors'

# ╔═╡ 6bf187b7-0b0b-43ae-80b3-3bc886b94ce0
sum([B₁.values[i]*B₁.vectors[:,i]*B₁.vectors[:,i]' for i=1:size(A₁,1)])

# ╔═╡ 5ab73935-218b-4840-9f76-1c580ac9ac7d
md"""
## Rešavanje algebarskih problema pomoću svojstvenih vrijednosti i vektora

Riješimo problem (prema [Logan, Applied Mathematics, str. 205][Log97])

$$
Ax=\mu x + f.$$

Neka je $A$ simetrična, $A=U\Lambda U^T$ i $\mu\neq \lambda_i$. Stupci matrice $U$ su ortogonalni i tvore bazu $n$-dimenzionalnog prostora, odnosno svaki vektor se može prikazati kao njihova linearna kombinacija:

$$
x=\sum_{i=1}^n c_i u_i, \quad  f=\sum_{i=1}^n f_i u_i.$$

Imamo

$$
A\cdot \big(\sum c_i u_i\big)=\mu \big(\sum c_i u_i\big) + \sum f_i u_i,$$

odnosno,

$$
\sum c_i \lambda_i u_i=\mu \big(\sum c_i u_i\big) + \sum f_i u_i.$$

Izjednačavanje koeficijenata daje

$$
c_i\lambda_i= \mu c_i + f_i$$

pa je 

$$
c_i=\frac{f_i}{\lambda_i-\mu}.$$

(Vidi J. David Logan, 'Applied Mathematics', 2nd Edition, Wiley, New York, 1997.)
"""

# ╔═╡ 809c18d2-07b4-4787-ad07-7199a36b3b83
md"""
Riješimo problem $Ax=\mu x+f$ za $\mu=2$:

$$
A=\begin{bmatrix} 1&2 &3&4 \\ 2 &5&6&7 \\ 3&6&8&9 \\ 4&7&9&10\end{bmatrix},\quad
f=\begin{bmatrix}1 \\ 2\\ 1\\ 2\\\end{bmatrix}.$$
"""

# ╔═╡ 84695256-71d9-49ac-9736-92ec7c87797f
begin
	μ=2
	A₃=[1 2 3 4;2 5 6 7;3 6 8 9;4 7 9 10]
	f₃=[1;2;1;2]
end

# ╔═╡ 2b0735aa-6178-4fee-b21a-615a9a7a4d04
A₃

# ╔═╡ 4cd59738-1f83-49da-860a-ebce6a37ac94
B₃=eigen(A₃);

# ╔═╡ a2dd8de3-18a0-4c5f-9950-985eede74fb7
U₃=B₃.vectors

# ╔═╡ 70d9ff94-71d6-41af-a952-3b3bad42ad6f
md"""
Izračunajmo koeficijente vektora $f$ u bazi $U$
"""

# ╔═╡ 3f5a2e21-d3ae-4be9-9227-b55a1debe7f1
md"""
Izračunajmo koeficijente $c$ rješenja $x$ u bazi $U$ 
"""

# ╔═╡ c7c533e8-c191-418f-91aa-205eceb7da23
md"""
# Linearni operatori

__Operator__ je preslikavnje $L:X\to X$ gdje je $X$ vektorski prostor. 


Neka su 
$x,y\in X$ i $\alpha, \beta \in \mathbb{R}$.

Operator je __linearan__ ako je __aditivan__,

$$
L(x+y)=L(x)+L(y),$$

i __homogen__,

$$
L(\alpha x)=\alpha L(x).$$

Oba svojstva zajedno možemo pisati kao 

$$
L(\alpha x+\beta y)=\alpha L(x) + \beta L(y).$$


### Matrica je linearni operator na skupu vektora

Uz definiciju

$$
A(x)\equiv A\cdot x,$$

vrijedi

$$
A(x+y)=A(x)+A(y),\quad A(\alpha x)=\alpha A(x), $$

odnosno

$$
A(\alpha x+\beta y)=\alpha A(x)+\beta A(y).$$
"""

# ╔═╡ 389409ff-62eb-4375-bf51-00368191069e
md"""
Detalji:

$$A(\alpha x+\beta y) =A\cdot (\alpha x+\beta y)=\alpha A\cdot x+\beta A\cdot y=\alpha A(x)+\beta A(y).$$
"""

# ╔═╡ 9d5e1a3f-1003-4e83-9d8b-e7be5ed1ddc9
md"""
# Skalarni produkt, norma, ortogonalnost i baza

Općenito, __norma__ na vektorskom prostoru $X$ je svaka funkcija $\| \phantom{x} \| : X\to \mathbb{R}$ sa sljedećim svojstvima:

1.  $\| x\|=0 \Leftrightarrow x=0$
2.  $\| \lambda x\|=|\lambda| \|x\|$
3.  $\| x+y\| \leq \|x\|+\|y\|$ (nejednakost trokuta)

__Skalarni produkt__ na vektorskom prostoru $X$ je svako preslikavanje $\cdot : X\times X \to \mathbb{R}$ sa sljedećim svojstvima:

1.  $x\cdot x\geq 0$
1.  $x\cdot x=0 \Leftrightarrow x=0$
2.  $x\cdot y=y\cdot x$
3.  $(\alpha x)\cdot y =\alpha (x\cdot y)$
3.  $(x+y)\cdot z=x\cdot z+y \cdot z$

Ukoliko je na vektorskom prostoru definiran skalarni produkt, normu možemo definirati kao

$$
\|x\|=\sqrt{x\cdot x}.$$

Također, ako je $x \cdot y=0$ kažemo da su vektori $x$ i $y$ __međusobno ortogonalni (okomiti)__.


## Konačno-dimenzionalni prostor

Neka su zadani vektori $x,y\in \mathbb{R}^n$. Definiramo sljedeće:

__skalarni produkt__: $(x,y)=x\cdot y=\displaystyle\sum_{i=1}^n x_i y_i$,

__norma__: $\|x \|=\sqrt{(x,x)}=\sqrt{\displaystyle\sum_{i=1}^n x_i x_i}=\sqrt{\displaystyle\sum_{i=1}^n x_i^2}$,

__ortogonalnost__: $x\perp y \Leftrightarrow (x,y)=0$,

__baza__: Skup od $n$ vektora, $x_1,x_2,\ldots, x_n$ je __potpun__ (baza) ako za svaki vektor $y$ vrijedi

$$
y=\displaystyle\sum_{i=1}^n \xi_i x_i.$$

Ako su, dodatno, vektori $x_i$ međusobno ortogonalni, onda je 

$$
\xi_j=\frac{(y,x_j)}{(x_j,x_j)}\equiv \frac{(y,x_j)}{\|x_j\|^2}.$$
"""

# ╔═╡ cb1db8f1-a43a-4ef4-a521-1df9031c33ce
U₃

# ╔═╡ aa87e87a-53ca-4caa-af7d-82d8a92c8763
md"""
# Vektorski prostor funkcija

Neka su zadane funkcije $f,g\in C[a,b]$, gdje je $C[a,b]$ skup svih funkcija neprekidnih na intervalu $[a,b]$. 

__Napomena__. Umjesto skupa $C[a,b]$ može se uzeti i neki skup, na primjer, skup svih kvadratno integrabilnih funkcija na intervalu $[a,b]$ kojeg označavamo s  $L^2[a,b]$.

Definirajmo __skalarni produkt__: 

$$(f,g)=f\cdot g=\displaystyle \int\limits_{a}^b f(x) g(x) \, dx.$$

Odavde slijede definicije:

__norma__: $\|f \|=\sqrt{(f,f)}=\sqrt{\displaystyle\int\limits_{a}^b f(x)\cdot f(x)\, dx}=
\sqrt{\displaystyle\int\limits_{a}^b f^2(x)\, dx}$

__ortogonalnost__: $f\perp g \Leftrightarrow (f,g)=0$

__baza__: Skup od $\infty$ funkcija, $f_1,f_2,\ldots$, je __potpun__ (baza) ako za svaku funkciju $g$ vrijedi

$$
g(x)=\displaystyle\sum_{i=1}^\infty \xi_i f_i(x).$$

Ukoliko su, dodatno, funkcije $f_i$ međusobno ortogonalne, tada je 

$$
\xi_j=\frac{(g,f_j)}{(f_j,f_j)}\equiv \frac{(g,f_j)}{\|f_j\|^2}.$$
"""

# ╔═╡ 3a2c17c9-e058-472b-97b7-e91d0046e5aa
md"""
## Numeričko i simboličko računanje

Julia ima više paketa pomoću kojih možemo računati određene integrale. 

Najjednostavnija je funkcija `quadgk()` iz paketa 
[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) (numeričko računanje).

Ovdje ćemo skalarni produkt definirati pomoću funkcije
`integrate()` iz paketa `SymPy.jl` (simboličko računanje).
"""

# ╔═╡ 89b5d46f-203c-43db-890e-04c04261c64d
#?integrate

# ╔═╡ 03f67024-69a0-49bd-9633-276b9faa3a82
md"""
## Fourierov red

Promotrimo periodične funkcije s periodom $2\pi$ na intervalu $[-\pi,\pi]$.

Funkcije

$$
1, \sin x, \cos x, \sin(2x), \cos(2x), \sin(3x),\cos(3x), \ldots$$

su međusobno ortogonalne, Vrijedi $\|1\|=\sqrt{2\pi}$, a norma svih ostalih funkcija je $\sqrt{\pi}$. Skup je potpun, odnosno svaka periodična funkcija $f$ se može prikazati kao 

$$
f(x)=\sum_{i=0}^\infty \xi_i f_i(x), \quad \xi_i=\frac{(f,f_i)}{(f_i,f_i)},$$

u smislu teorema o konvergenciji Fourierovog reda. 
Ovo su standardne formule za razvoj funkcije u Fourierov red.
"""

# ╔═╡ 3a116293-f9ca-46bd-aab7-4a2bc3e8d692
# Definirajmo simboličku varijablu x
x=Sym("x")

# ╔═╡ 88d2761c-5f71-4aba-92db-5dd65ec4bf73
begin
	# Definirajmo skalarni produkt
	import LinearAlgebra.⋅
	⋅(f,g,a,b)=integrate(f*g,(x,a,b))
end

# ╔═╡ 38b4daf6-4c03-4914-a08b-a7ce1c9baf88
fU=[f₃⋅U₃[:,i] for i=1:4]

# ╔═╡ 62254437-6a35-4813-b5dc-758619199ac3
# Provjera
U₃*fU

# ╔═╡ 9421c609-c322-4076-b0c0-f3545bcbd373
c=fU./(B₃.values.-μ)

# ╔═╡ 2c7c297d-3362-4901-b137-a5dfd707ab8b
# Rješenje
x₃=sum([c[i]*U₃[:,i] for i=1:4])

# ╔═╡ 11c369af-58df-44a4-ba69-64019714e17d
# Provjera
A₃*x₃.-μ*x₃.-f₃

# ╔═╡ 1db60195-4666-4efd-a9c1-ae66594c89f2
U₃*c

# ╔═╡ e2413b82-40f9-4c82-9584-43d5f642c47f
# Primjer za vektore - ortogonalnost i norma
U₃[:,2]⋅U₃[:,3], U₃[:,3]⋅U₃[:,3]

# ╔═╡ 66e0cb40-9a63-427a-af09-0b410f99784d
begin
	# Baza
	n₃=size(A₃,1)
	z=rand(n₃)
	# Računamo koeficijente po bazi stupaca od U
	ξ=Array{Float64}(undef,n₃)
	for i=1:n₃
	    ξ[i]=z⋅U₃[:,i]
	end
	# Provjera
	y=sum([ξ[i]*U₃[:,i] for i=1:n₃])
	[z y]
end

# ╔═╡ 93e28bed-a897-468d-9d88-ee9ce563a33b
# Provjerimo ortogonalnost funkcija
⋅(sin(x),sin(x),-π,π), ⋅(sin(2*x),cos(3*x),-π,π)

# ╔═╡ 82472731-ba8a-4c8f-8ebd-01f5d45a7f6e
m,n = symbols("m,n", integer=true, positive=true)

# ╔═╡ 7bbee2f4-2abd-4c20-a84e-22b73ec53c09
⋅(cos(m*x),cos(n*x),-PI,PI)

# ╔═╡ 2cc6f672-e827-46d1-b08e-8de0c0cbf71b
md"""
__Primjer:__ Razvijmo funkciju definiranu formulom $f(x)=x^2$ na intervalu $[-1,1]$ u Fourier-ov red.
"""

# ╔═╡ f96d7049-d0f6-4176-a238-ffba3a9e08fa
f=x^2

# ╔═╡ c742c464-b181-4886-86b0-8a028d22048d
f₀=x/x

# ╔═╡ 7277a5ab-1ab0-430c-9f29-90945c1330c8
⋅(f₀,f₀,-1,1)

# ╔═╡ b1477d47-a9d8-4041-9ebe-5cb33bfe074c
a₀=⋅(f,f₀,-1,1)/⋅(f₀,f₀,-1,1)

# ╔═╡ e5b3b137-dcda-4590-be90-ba8373fbc313
# Funkcija je parna pa imamo samo članove uz kosinuse
aₙ=⋅(f,cos(n*PI*x),-1,1)/⋅(cos(n*PI*x),cos(n*PI*x),-1,1)

# ╔═╡ 5ba75774-f60c-4292-8dca-0bab561ec9ce
bₙ=⋅(f,sin(n*PI*x),-1,1)/⋅(sin(n*PI*x),sin(n*PI*x),-1,1)

# ╔═╡ b4d6d3ed-ff18-4b19-9bed-c7847947cc98
# Na primjer
aₙ(2)

# ╔═╡ 4218c90b-00b6-47b8-8805-5e7c8b9c5d64
aₙ(3)

# ╔═╡ 64096fe3-6ec0-4f03-9c1a-5e758b2e43fb
plot([x->x^2,x->a₀+sum([aₙ(i)*cos(i*PI*x) for i=1:10])],-1,1,label=["Funkcija" "Fourierov red"])

# ╔═╡ 5ff59714-d411-426a-a995-c1416f8c6055
md"""
# Diferencijalni problem svojstvenih vrijednosti

Skup $C^2[a,b]$ je skup svih funkcija koje na intervalu $[a,b]$ imaju dvije neprekidne derivacije. 


## Operator druge derivacije

Operator druge derivacije $A\equiv\displaystyle\frac{d^2}{dx^2}$ je linearan operator. 

Riješimo problem svojstvenih vrijednosti

$$
\frac{d^2}{dx^2} \Phi=\lambda \Phi, \quad 0<x<l,\quad \Phi(0)=\Phi(l)=0,\quad \Phi\neq 0.$$

Razlikujemo slučajeve $\lambda=0$, $\lambda<0$ i $\lambda>0$.

__Slučaj__ $\lambda=0$. 

Vrijedi $\Phi(x)=a x+b$. Iz rubnog uvjeta $\Phi(0)=0$ slijedi $b=0$ pa je $\Phi(x)=ax$. Iz rubnog uvjeta $\Phi(l)=0$ slijedi $al=0$ pa je i $a=0$. Dakle, $\Phi(x)=0$, što ne može biti svojstvena funkcija, pa $\lambda=0$ nije svojstvena vrijednost.

__Slučaj__ $\lambda>0$. 

Vrijedi (vidi [Linearne diferencijalne jednadžbe drugog reda s konstantnim koeficijentima][Mat2])

$$
\Phi(x)=a e^{\displaystyle\sqrt{\lambda}x}+ b e^{-\displaystyle\sqrt{\lambda}x}.$$

Iz rubnog uvjeta $\Phi(0)=0$ slijedi $a+b=0$ pa je $b=-a$. 

Iz rubnog uvjeta $\Phi(l)=0$ slijedi 

$$
a\big(e^{\displaystyle\sqrt{\lambda}l}-e^{-\displaystyle\sqrt{\lambda}l}\big)=0$$

pa je $a=0$. Dakle, $\Phi(x)=0$, što ne može biti svojstvena funkcija, pa niti jedna $\lambda>0$ nije svojstvena vrijednost.

__Slučaj__ $\lambda<0$.

Vrijedi

$$
\Phi(x)=a \sin (\sqrt{-\lambda}x)+b \cos (\sqrt{-\lambda}x).$$

Iz rubnog uvjeta $\Phi(0)=0$ slijedi $b=0$ pa je $\Phi(x)=a\sin(\sqrt{-\lambda}x)$. 

Iz rubnog uvjeta $\Phi(l)=0$ slijedi 

$$
a \sin(\sqrt{-\lambda}l)=0$$

pa je ili $a=0$, što opet ne daje svojstvenu funkciju, ili 

$$
\sqrt{-\lambda}l=n\pi, \quad n\in\mathbb{N}.$$

Dakle, svojstvene vrijednosti su 

$$
\lambda_n=-\frac{n^2\pi^2}{l^2}, \quad n\in\mathbb{N},$$

a pripadne svojstvene funkcije su

$$
\Phi_n(x)=\sin \big(\frac{n\pi}{l}x\big).$$


Funkcije $\Phi_n(x)$ su međusobno ortogonalne i čine bazu promatranog prostora.


Ponovite postupak rješavanja linearna diferencijalne jednadžbe drugog reda s konstantnim koeficijentima ➡ [Matematika 2](http://lavica.fesb.unist.hr/mat2/predavanja/node95.html).

"""

# ╔═╡ c9931d6a-6ea8-4fc3-a9cd-e1629b89993b
l = symbols("l", real=true, positive=true)

# ╔═╡ 58e3a214-a053-4379-8f92-6f4197dc7a92
# Provjerimo ortogonalnost
⋅(sin(n*PI*x/l),sin(m*PI*x/l),0,l)

# ╔═╡ 31fb102e-d8e3-4e17-b0fe-2d1d3dd99011
md"""
# Regularni Sturm-Liouvilleov problem (SLP)

Problem glasi:

$$
\begin{aligned}
&A(\Phi) \equiv -(p(x)\,\Phi')'+q(x)\,\Phi = \lambda\, \Phi, \quad a\leq x\leq b,\\
&\alpha_1\Phi(a)+\alpha_2\Phi'(a)=0,\\
&\beta_1 \Phi(b)+\beta_2\Phi'(b)=0,
\end{aligned}$$

gdje je 

$$\Phi\in C^2[a,b],\quad p\in C^1[a,b],\quad q\in C^0[a,b],\quad 
\alpha_i,\beta_i\in\mathbb{R}.$$

Operator $A$ je linearan (provjerite!).

__Teorem__. Za regularni SLP vrijedi:

1. Postoji beskonačno mnogo svojstvenih vrijednosti $\lambda_n$, $n=1,2,3,\ldots$, koje su sve realne i vrijedi

$$\lim\limits_{n\to\infty} |\lambda_n|=\infty.$$

2. Svojstvene funkcije koje odgovaraju različitim svojstvenim vrijednostima su ortogonalne. 

3. Skup svih svojstvenih funkcija $\Phi_1,\Phi_2,\Phi_3,\ldots$ je potpun u smislu da se svaka funkcija $f\in L^2[a,b]$ može razviti u red

$$f(x)=\sum_{n=1}^\infty \xi_n \Phi_n(x), \quad \xi_n=\frac{(f,\Phi_n)}{(\Phi_n,\Phi_n)}$$ 

koji konvergira u $L^2[a,b]$. 

Konvergencija u $L^2[a,b]$ znači

$$
\big\|f-\sum_{n=1}^N \xi_n\Phi_n\big\|^2\equiv \int\limits_a^b 
\big(f-\sum_{n=1}^N \xi_n\Phi_n\big)^2 dx \to 0 \quad \textrm{kada} \quad  N\to\infty.$$

Na primjer, teorem vrijedi za regularni SLP iz prethodnog primjera, gdje je 

$$p(x)=-1,\quad q(x)=0,\quad a=0,\quad b=l, \\  
\alpha_1=1,\quad \alpha_2=0,\quad \beta_1=1,\quad \beta_2=0.$$

_Dokaz 2. tvrdnje_ (prema Logan, Applied Mathematics, 2. izdanje, str. 209): Neka su $\lambda$ i $\mu$ dvije različite svojstvene vrijednosti sa svojstvenim funkcijama $\phi$ i $\psi$, redom. Tada vrijedi

$$
\begin{aligned}
-(p\phi')'+q\phi&=\lambda\phi, \\
-(p\psi')'+q\psi&=\mu\psi.
\end{aligned}$$

Pomnožimo prvu jednadžbu sa $\psi$ i drugu sa $\phi$ te ih oduzmimo:

$$
\phi(p\psi')'-\psi(p\phi')'=(\lambda-\mu)\phi\psi.
$$

Integriranje od $a$ do $b$ daje

$$
\int\limits_a^b (\phi(p\psi')'-\psi(p\phi')')\, dx=(\lambda-\mu) (\phi,\psi).
$$

Parcijalna integracija daje

$$
\int\limits_a^b \phi(p\psi')'\, dx = \left\{ {u=\phi, \quad du=\phi'\, dx  \atop dv=(p\psi')'\, 
dx, \quad v=p\psi' } \right\}
=\phi(p\psi')\big|_a^b -\int\limits_a^b p\psi'\phi' \, dx,
$$

i, slično,

$$
\int\limits_a^b \psi(p\phi')'\, dx = \psi(p\phi')\big|_a^b -\int\limits_a^b p\psi'\phi' \, dx.
$$

Dakle, 

$$
p(\phi\psi'-\psi\phi')\big|_a^b=(\lambda-\mu)(\phi,\psi).
$$

Iz rubnih uvjeta slijedi da je lijeva strana jednaka nuli: na primjer, 
ako su sva dijeljenja definirana, onda je 

$$
\frac{\phi(a)}{\phi'(a)}=-\frac{\alpha_2}{\alpha_1}=\frac{\psi(a)}{\psi'(a)}
$$

pa je $\phi(a)\psi'(a)-\phi'(a)\psi(a)=0$. Slično se analiziraju i ostali slučajevi.

Dakle, 

$$0=(\lambda-\mu)(\phi,\psi).$$

Ako je $\lambda\neq \mu$, onda je $(\phi,\psi)=0$, odnosno,
$\phi\perp\psi$. 

Primjere rješavanja regularnog SLP dati ćemo kasnije.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
Plots = "~1.22.4"
PlutoUI = "~0.7.15"
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
git-tree-sha1 = "633f8a37c47982bff23461db0076a33787b17ecd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.15"

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
# ╠═55e86888-4b1b-4f6d-9bd8-c839e76fe381
# ╠═c5e5c059-ba8d-4d1d-b171-9c7673761320
# ╟─68499249-33bd-4c9b-9ba4-072e096c4820
# ╟─5bad68c8-ccbf-4fdd-b22b-8fa78a3b95bb
# ╠═ab8033eb-60b2-4cbf-9dae-15c6753bba76
# ╠═21d8db39-7058-49d7-9174-f1b88ff7c062
# ╟─d6a6f2b1-d64f-4d7b-94d2-f98fa36c0996
# ╠═b2deff23-48d4-4495-bb73-35dfa0f55487
# ╠═a46ec17d-6c24-455a-b049-e521f79b0667
# ╟─ac54aa9d-b480-434d-83fc-69f2f71044bd
# ╠═52597f5d-faed-438f-ae5a-554b63a88e78
# ╠═3f53d041-f8b1-4483-92de-93e918a3c2e3
# ╠═1f90b141-4403-477d-b57d-56493559b275
# ╟─b1309118-5f38-412b-8e51-832feacfb0e2
# ╠═7156ab8c-5ec6-46d0-ba7f-65cc231f2fe7
# ╠═7f0b797f-57a0-4287-abf7-028c357e7da2
# ╠═299ff3f2-24bb-46d8-bbea-c2e050282779
# ╠═52faa36a-b161-4d4e-a474-c49ec7f9f245
# ╠═cdbfa73e-9783-4256-9035-4afe0ae54cbb
# ╠═1ab3bd33-27b6-4f95-ae41-6e2689b6ddc0
# ╠═ab619c12-f92f-42d2-877a-407e17c3d671
# ╠═79f3c8f8-30e6-423e-8a2e-a402f9f8b972
# ╠═5df96e85-43f3-46e2-b1b0-739637782191
# ╠═73f65db7-9958-4044-990f-db91e0c4ea1e
# ╠═d2038ec5-8957-4c9f-b3eb-89ed0a9fbc16
# ╠═3c3bf03a-1473-45a7-9aae-f4b33f6268b0
# ╟─4789f475-d239-4488-93d4-0989e2489835
# ╠═8b032f4c-7ef2-4582-9e1e-ec9ec4e45d80
# ╟─80a7e0ba-87a6-48ae-84b1-9aeb354f9a6a
# ╠═ed37fc89-01a7-497d-a4ca-ba584ef21781
# ╠═12964840-d516-423d-b0fd-17baea64fc58
# ╠═d4aba6a0-3012-4e40-8495-9957dce33a1c
# ╠═a9160ec4-50ff-4067-bfde-61e2234af4d6
# ╠═27423aed-47a5-473c-8016-98a39b8fd7e2
# ╠═89fe6b2f-6106-4695-92fb-62e57d1bae0b
# ╠═d477b32b-f85c-4f3c-9f96-8726725a7b3a
# ╠═6bf187b7-0b0b-43ae-80b3-3bc886b94ce0
# ╟─5ab73935-218b-4840-9f76-1c580ac9ac7d
# ╟─809c18d2-07b4-4787-ad07-7199a36b3b83
# ╠═84695256-71d9-49ac-9736-92ec7c87797f
# ╠═2b0735aa-6178-4fee-b21a-615a9a7a4d04
# ╠═4cd59738-1f83-49da-860a-ebce6a37ac94
# ╠═a2dd8de3-18a0-4c5f-9950-985eede74fb7
# ╟─70d9ff94-71d6-41af-a952-3b3bad42ad6f
# ╠═38b4daf6-4c03-4914-a08b-a7ce1c9baf88
# ╠═62254437-6a35-4813-b5dc-758619199ac3
# ╟─3f5a2e21-d3ae-4be9-9227-b55a1debe7f1
# ╠═9421c609-c322-4076-b0c0-f3545bcbd373
# ╠═2c7c297d-3362-4901-b137-a5dfd707ab8b
# ╠═11c369af-58df-44a4-ba69-64019714e17d
# ╠═1db60195-4666-4efd-a9c1-ae66594c89f2
# ╟─c7c533e8-c191-418f-91aa-205eceb7da23
# ╟─389409ff-62eb-4375-bf51-00368191069e
# ╟─9d5e1a3f-1003-4e83-9d8b-e7be5ed1ddc9
# ╠═cb1db8f1-a43a-4ef4-a521-1df9031c33ce
# ╠═e2413b82-40f9-4c82-9584-43d5f642c47f
# ╠═66e0cb40-9a63-427a-af09-0b410f99784d
# ╟─aa87e87a-53ca-4caa-af7d-82d8a92c8763
# ╟─3a2c17c9-e058-472b-97b7-e91d0046e5aa
# ╠═89b5d46f-203c-43db-890e-04c04261c64d
# ╠═88d2761c-5f71-4aba-92db-5dd65ec4bf73
# ╟─03f67024-69a0-49bd-9633-276b9faa3a82
# ╠═3a116293-f9ca-46bd-aab7-4a2bc3e8d692
# ╠═93e28bed-a897-468d-9d88-ee9ce563a33b
# ╠═82472731-ba8a-4c8f-8ebd-01f5d45a7f6e
# ╠═7bbee2f4-2abd-4c20-a84e-22b73ec53c09
# ╟─2cc6f672-e827-46d1-b08e-8de0c0cbf71b
# ╠═f96d7049-d0f6-4176-a238-ffba3a9e08fa
# ╠═c742c464-b181-4886-86b0-8a028d22048d
# ╠═7277a5ab-1ab0-430c-9f29-90945c1330c8
# ╠═b1477d47-a9d8-4041-9ebe-5cb33bfe074c
# ╠═e5b3b137-dcda-4590-be90-ba8373fbc313
# ╠═5ba75774-f60c-4292-8dca-0bab561ec9ce
# ╠═b4d6d3ed-ff18-4b19-9bed-c7847947cc98
# ╠═4218c90b-00b6-47b8-8805-5e7c8b9c5d64
# ╠═64096fe3-6ec0-4f03-9c1a-5e758b2e43fb
# ╟─5ff59714-d411-426a-a995-c1416f8c6055
# ╠═c9931d6a-6ea8-4fc3-a9cd-e1629b89993b
# ╠═58e3a214-a053-4379-8f92-6f4197dc7a92
# ╟─31fb102e-d8e3-4e17-b0fe-2d1d3dd99011
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
