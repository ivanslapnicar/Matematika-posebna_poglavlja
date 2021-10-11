### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ a20ab82c-881d-4a01-97cc-b18fccb76af0
using PlutoUI

# ╔═╡ 16e6736a-084f-4a86-a0e4-988c22778768
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 4c99ed0d-147a-4a4a-a634-6d0debca4c5b
md"""
# Jednadžba difuzije


> zakon očuvanja + jednadžba stanja (konstitutivna jednadžba) + rubni uvjeti


# Zakon očuvanja 

Neka je $I$ dio cilindrične cijevi od $a$ do $b$, neka je $A$ površina presjeka i neka je

$$
u(x,t)$$

gustoća ili koncentracija _nečega/tvari/energije_ na mjestu $x$ u trenutku $t$.

![cijev](https://ivanslapnicar.github.io/Matematika-posebna_poglavlja/files/fig.png)

Količina tvari $u$ unutar $I$ je gustoća $\times$ volumen (površina baze $\times$ visina):

$$
\int\limits_a^b u(x,t) A \, dx.$$

Neka se tvar kreće i neka je 

$$
\phi(x,t)$$

tok kroz presjek $A$ na mjestu $x$ u trenutku $t$ po jedinici površine. Ako je $\phi(x,t)>0$, tok je prema pozitivnom smjeru $x$-osi.

Stopa po kojoj tvar ulazi u $I$ jednaka je razlici količine koja ulazi na mjestu $a$ i količini koja izlazi na mjestu $b$:

$$
A\phi(a,t)-A\phi(b,t).$$

Na primjer, ako se radi o vodi, ova razlika je nula.

Neka tvar nastaje ili nestaje pomoću izvora ili uvira po stopi

$$
f(x,t,u)$$

po jedinici površine. Stopa po kojoj tvar nastaje/nestaje unutar $I$ je

$$
\int\limits_a^b f(x,t,u)A\, dx.$$

Vremenska promjena količine tvari unutar $I$ jednaka je zbroju stope po kojoj tvar ulazi u $I$ i stope po kojoj tvar nastaje/nestaje u $I$:

$$
\frac{d}{dt} \int\limits_a^b u(x,t) A \, dx =A\phi(a,t)-A\phi(b,t) + \int\limits_a^b f(x,t,u)A\, dx.
$$

Možemo skratiti s $A$ i dobili smo __zakon očuvanja u integralnom obliku__:

$$
\frac{d}{dt} \int\limits_a^b u(x,t)\, dx =\phi(a,t)-\phi(b,t) + \int\limits_a^b f(x,t,u)\, dx.
$$

Pretpostavimo da su $u$ i $\phi$ neprekidno diferencijabilne (glatke) funkcije. Na lijevu stranu jednadžbe primijenimo [postupak deriviranja pod znakom integrala (Leibnitz-ovu formulu)](http://lavica.fesb.hr/mat2/predavanja/node81.html), a desnu stranu jednadžbe primijenimo [Newton-Leibnitzovu formulu](http://lavica.fesb.unist.hr/mat2/predavanja/node27.html) u obrnutom smislu, pa imamo

$$
\int\limits_a^bu_t(x,t) \, dx = -\int\limits_a^b \phi_x(x,t)\, dx + \int\limits_a^b f(x,t,u)\, dx,$$

odnosno

$$
\int\limits_a^b [u_t(x,t)+\phi_x(x,t)-f(x,t,u)]\, dx=0.$$

Jednakost vrijedi za sve $a$ i $b$ pa je podintegralna funkcija jednaka nuli, odnosno vrijedi __zakon očuvanja u diferencijalnom obliku__:

$$
u_t(x,t)+\phi_x(x,t)-f(x,t,u)=0.$$
"""

# ╔═╡ bf50b209-2e29-48ac-82a7-56dd5bdde7ad
md"""
# Jednadžba stanja

__Jednadžba stanja__ ili __konstitutivna jednadžba__ je empirijska. 

__Fickov zakon__ kaže da je tok proporcionalan promjeni koncentracije:

$$
\phi(x,t)=-D u_x(x,t),$$

gdje je $D$ __konstanta difuzije__. Na primjer, toplina ide s mjesta veće koncentracije na mjesto manje koncentracije. 

Dakle,

$$
\phi_x(x,t)=-D u_{xx}(x,t).$$


## Primjeri

__Difuzijska jednadžba__ (bez izvora) glasi:

$$
u_t(x,t)-Du_{xx}=0,$$

a __reakcijska-difuzijska jednadžba__ glasi:

$$
u_t(x,t)-Du_{xx}=f(x,t,u).$$

Ponekad i $\phi$ ovisi o $u$, $\phi(x,t,u)$.

Primjer za oba slučaja je __toplinska jednadžba__. 

__Fisherova jednadžba__ nastaje iz [logističke jednadžbe](http://lavica.fesb.unist.hr/mat2/predavanja/node88.html):

$$
u'(t)=k u(t)\left( 1-\frac{u(t)}{K}\right).$$

Ako populacija ovisi i o mjestu te ako imamo tok (npr. kukci se sele i nema ih svugdje jednako), jednadžba glasi

$$
u_t(x,t)+\phi_x=ku(x,t)\bigg( 1-\frac{u(x,t)}{K}\bigg).$$

Ako pretpostavimo da vrijedi Fickov zakon (npr. kukci se sele u manje naseljena područja), imamo

$$
u_t-Du_{xx}=ku\big( 1-\frac{u}{K}\big).$$
"""

# ╔═╡ da5b4f30-cf98-4f10-bd32-6fb863191456
md"""
# Rubni uvjeti

Prirodno je zadati gustoću u početnom trenutku $t=0$,

$$
u(x,0)=f(x), \quad a\leq x\leq b,$$

te gustoću ili njenu promjenu na rubovima.

__Dirichletovi__ ili __geometrijski uvjeti__ glase

$$
u(a,t)=g(t), \quad u(b,t)=h(t),\quad t>0.$$

Na primjer, $u(a,t)=1$ znači da je temperatura fiksirana na rubu.

__Neumannovi__ ili __prirodni uvjeti__ glase:

$$
u_x(a,t)=g(t), \quad u_x(b,t)=h(t),\quad t>0.$$

Na primjer, $u_x(a,t)=0$ znači da je rub izoliran. 

__Robinovi__ ili __mješoviti uvjeti__ glase:

$$
\alpha_1 u(a,t)+\alpha_2 u_x(a,t)=g(t), \quad \beta_1 u(b,t)+\beta_2 u_x(b,t)=h(t), \quad t>0.$$

Na primjer, prema [Newtonovom zakon hlađenja](http://lavica.fesb.unist.hr/mat2/predavanja/node89.html) je tok topline kroz rub proporcionalan razlici temperature ruba i okolne temperature:

$$
u_x(a,t)=\alpha(u(a,t)-T_{ambienta}(t)).$$

Moguće su i kombinacije uvjeta u različitim rubovima.
"""

# ╔═╡ 97219b4d-9430-48e1-8f74-d9252c61a283
md"""
# Jedinstvenost rješenja

> Ako dokažemo da problem rubnih vrijednosti ima jedinstveno rješenje, onda znamo da je rješenje koje smo dobili bilo kojom metodom upravo rješenje koje tražimo. 

Na primjer:

__Teorem__ Rješenje problema rubnih vrijednosti

$$
\begin{aligned}
& u_t-Du_{xx}=0,\quad D>0, \quad 0<x<l, \quad 0<t<T,\\
&u(x,0)=f(x),\quad 0<x<l,\\
&u(0,t)=g(t),\quad u(l,t)=h(t),\quad 0<t<T,
\end{aligned}$$

pri čemu su funkcije $f$, $g$ i $h$ neprekidne, je _jedinstveno_ za svaki $T>0$.

_Dokaz_ : Dokaz je kontradikcijom. Pretpostavimo da postoje dva rješenje, $u_1$ i $u_2$, i definirajmo

$$
w(x,t)=u_1(x,t)-u_2(x,t).$$

Tada je $w$ rješenje problema rubnih vrijednosti

$$
\begin{aligned}
& w_t-Dw_{xx}=0,\\
&w(x,0)=f(x)-f(x)=0,\\
&w(0,t)=0,\quad w(l,t)=0.
\end{aligned}$$

Definirajmo __integral energije__:

$$
E(t)=\int\limits_0^l w^2(x,t)\, dx.$$

Vrijedi $E(t)\geq 0$ i $E(0)=\int_0^0 w^2 dx=0$. Deriviranje pod znakom integrala, uvrštavanje jednadžbe i parcijalna integracija daju:

$$
\begin{aligned}
E'(t)&=\int\limits_0^l 2w(x,t) w_t(x,t)\, dx\\
&=2D\int\limits_0^l w(x,t) w_{xx}(x,t)\, dx\\
&=2D\bigg[ w(x,t)w_x(x,t)\big|_0^l - \int\limits_0^l w_{x}^2(x,t)\, dx\bigg]\\
&=2D\bigg[ w(l,t)w_x(l,t)-w(0,t)w_x(0,t) - \int\limits_0^l w_{x}^2(x,t)\, dx\bigg].
\end{aligned}$$

Prva dva izraza u zagradi su jednaka nuli zbog rubnih uvjeta, a integral nenegativan, pa je 

$$
E'(t)\leq 0,$$

odnosno $E(t)$ pada. Kako je $E(0)=0$ i $E(t)\geq 0$, zaključujemo da je $E(t)=0$ za svaki $t$. 

No, onda je i $w(x,t)=0$ pa je $u_1=u_2$ i rješenje je jedinstveno.


__Zadatak__: Dokažite jedinstvenost u slučaju Neumannovih rubnih uvjeta.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═a20ab82c-881d-4a01-97cc-b18fccb76af0
# ╠═16e6736a-084f-4a86-a0e4-988c22778768
# ╟─4c99ed0d-147a-4a4a-a634-6d0debca4c5b
# ╟─bf50b209-2e29-48ac-82a7-56dd5bdde7ad
# ╟─da5b4f30-cf98-4f10-bd32-6fb863191456
# ╟─97219b4d-9430-48e1-8f74-d9252c61a283
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
