### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 20dc641a-1cd4-4806-92e4-0fb85767754e
using PlutoUI

# ╔═╡ 31e0ee64-1b7c-4f26-beae-b0ce127aaf72
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ a1a4c928-af9a-40bb-9d2b-9458461e364f
md"""
# Zakon očuvanja u 3D


__zakon očuvanja + jednadžba stanja (konstitutivna jednadžba) + rubni uvjeti__


# Zakon očuvanja 

Izvod zakona očuvanja u trodimenzionalnom prostoru $\mathbb{R}^3$ vrlo je sličan izvodu jednodimenzionalnog slučaja.

Neka je $x\equiv(x,y,z)$ točka u $\mathbb{R}^3$. 
Neka je

$$
u(x,t)$$

skalarna gustoća ili __koncentracija__ _nečega/tvari/energije_ u točki $x$ u trenutku $t$ (količina po jedinici volumena).

Neka je $V\subset \mathbb{R}^3$ područje i neka je $\partial V$ rub područja koji je ili glatak ili se sastoji od konačno po djelovima glatkih ploha:

![područje](https://ivanslapnicar.github.io/Matematika-posebna_poglavlja/files/V.png)

Količina tvari unutar $V$ u trenutku $t$ jednaka je trostrukom integralu (gustoća $\times$ volumen):

$$
\int\limits_V u(x,t) \, dx,$$

pri čemu je $dx\equiv dx dy dz$ element volumena.

Neka se tvar kreće. U tri dimenzije tok može biti u bilo kojem smjeru pa je zadan vektorskim poljem

$$
\vec \phi(x,t).$$

Neka je $\vec n(x)$ jedinični vektor vanjske normale na područje $V$ u točki $x$.
Ukupan tok prema vani kroz rub $\partial V$ jednak je [plošnom integralu vektorskog polja](http://lavica.fesb.hr/mat3/predavanja/node20.html):

$$
\int\limits_{\partial V} \vec \phi(x,t)\cdot \vec n(x)\, dS,$$

gdje je $dS$ element površine $\partial V$.


Ako tvar nastaje ili nestaje pomoću izvora ili uvira po stopi

$$
f(x,t,u),$$

onda je stopa po kojoj tvar nastaje/nestaje unutar $V$ jednaka

$$
\int\limits_V f(x,t,u)\, dx.$$


__Zakon očuvanja u integralnom obliku__ glasi:

$$
\frac{d}{dt} \int\limits_V u(x,t)\, dx =-\int\limits_{\partial V} \vec \phi \cdot \vec n \, dS + \int\limits_V f(x,t,u)\, dx.$$

Pretpostavimo da su $u$ i $\vec \phi$ neprekidno diferencijabilne (glatke) funkcije. Na lijevu stranu jednadžbe primijenimo [postupak deriviranja pod znakom integrala (Leibnitz-ovu formulu)](http://lavica.fesb.hr/mat2/predavanja/node81.html), a na desnu stranu jednadžbe primijenimo [Teorem o divergenciji](http://lavica.fesb.hr/mat3/predavanja/node21.html), pa imamo

$$
\int\limits_V u_t(x,t) \, dx = -\int\limits_V \mathop{\mathrm{div}} \vec \phi(x,t)\, dx + \int\limits_V f(x,t,u)\, dx,$$

odnosno

$$
\int\limits_V [u_t(x,t)+\mathop{\mathrm{div}}\vec \phi(x,t)-f(x,t,u)]\, dx=0.$$

Jednakost vrijedi za proizvoljno područje $V$ pa je podintegralna funkcija jednaka nuli, odnosno vrijedi __zakon očuvanja u diferencijalnom obliku__:

$$
u_t(x,t)+\mathop{\mathrm{div}}\vec \phi(x,t)=f(x,t,u),\quad x\in V,\quad t>0.$$
"""

# ╔═╡ 3d4ac4e3-3bfa-4c19-834a-f91db0490ad4
md"""
# Jednadžba stanja

Kao i u jednodimenzionalnom slučaju, __jednadžba stanja__ ili __konstitutivna jednadžba__ je empirijska. 

__Fickov zakon__ kaže da je tok proporcionalan promjeni koncentracije, a koncentracija najbrže pada u smjeru suprotnom od smjera gradijenta, odnosno 

$$
\vec \phi(x,t)=-D \mathop{\mathrm{grad}} u(x,t),$$

gdje je $D$ __konstanta difuzije__. Dakle,

$$
\mathop{\mathrm{div}} \vec \phi(x,t)=-D \mathop{\mathrm{div}}\mathop{\mathrm{grad}}u(x,t)\equiv
-D \Delta u(x,t),$$

gdje je 

$$
\Delta=\nabla\cdot\nabla=\nabla^2=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}+\frac{\partial^2}{\partial z^2}$$

__Laplaceov operator__.

## Primjeri

__Reakcijsko-difuzijska jednadžba__ u 3D glasi:

$$
u_t-D\Delta u=f(x,t,u), \quad x\in V,\quad t>0.$$

Ako nema izvora, $f\equiv 0$, tada možemo tražiti stabilno stanje (__steady-state solution__) ili rješenje $u\equiv u(x)$ koje ovisi samo o položaju, a ne o vremenu. Takvo rješenje zadovoljava __Laplaceovu jednadžbu__:

$$
\Delta u=0, \quad x\in V.$$

Ako izvori ovise samo o položaju, $f\equiv f(x)$, stabilno stanje zadovoljava __Poissonovu jednadžbu__:

$$
\Delta u =-\frac{f}{D}, \quad x\in V.$$

Primjeri za oba slučaja su __statička električna polja__ određena nabojima koje je nalaze izvan $V$, odnosno unutar $V$. 
"""

# ╔═╡ 31f19a11-c08e-4d5a-a702-1621ba693e18
md"""
# Rubni uvjeti

Zadavanje početnih uvjeta nije prirodno jer mala promjena početnih uvjeta dovodi do velike promjene u rješenju. 


Možemo zadati __Dirichletov__ ili __geometrijski uvjet__, to jest gustoću na rubu:

$$
u=g(x), \quad x\in \partial V,$$

__Neumannov__ ili __prirodni uvjet__, to jest promjenu gustoće na rubu u smjeru vanjske normale:

$$
\frac{\partial u}{\partial \vec n}=g(x), \quad x\in \partial V,$$

ili __mješoviti uvjet__:

$$
\alpha(x) u  +\beta(x)\frac{\partial u}{\partial \vec n} =g(x) , \quad x\in \partial V.$$
"""

# ╔═╡ 70259200-0fc2-42b8-a471-d941347adb2f
md"""
# Jedinstvenost rješenja

Ako dokažemo da problem rubnih vrijednosti ima jedinstveno rješenje, onda znamo da je rješenje koje smo dobili bilo kojom metodom upravo rješenje koje tražimo. Navodimo tri primjera.


__Teorem.__ Neka je $g$ neprekidna na $\partial V$ i $f$ neprekidna na $V$. Rješenje problema

$$
\Delta u=f, \quad x\in V;\quad u=g,\quad x\in \partial V,$$

je jedinstveno.

_Dokaz_ : Pretpostavimo da postoje dva rješenje, $u_1$ i $u_2$ i definirajmo $w=u_1-u_2$.

Tada je 

$$
\Delta w=0, \quad x\in V;\quad w=0,\quad x\in \partial V.$$

Uvrštavanje $\phi=\psi=w$ u [prvi Greenov identitet](https://en.wikipedia.org/wiki/Green%27s_identities) daje 

$$
\int\limits_V \mathop{\mathrm{grad}} w \cdot \mathop{\mathrm{grad}} w \, dx =0.$$

Dakle, $\mathop{\mathrm{grad}} w=0$ na području $V$ pa je $w$ konstantno polje na području $V$. Zbog $w=0$ na rubu $\partial V$, na poručju $V$ vrijedi  $w=0$, odnosno $u_1=u_2$.


__Teorem.__ Neka je $g$ neprekidna na $\partial V$ i $f$ neprekidna na $V$. Rješenje problema

$$
\Delta u=f, \quad x\in V;\quad \frac{\partial u}{\partial \vec n}=g,\quad x\in \partial V,$$

zadovoljava

$$
\int\limits_Vf\, dx=\int\limits_{\partial V} g\, dS.$$

(U stabilnom stanju je tok tvari kroz rub jednak količni tvari koja nastaje unutar područja.)

_Dokaz_ : Tvrdnja slijedi uvrštavanjem  $\phi=u$ i $\psi=1$ u [drugi Greenov identitet](https://en.wikipedia.org/wiki/Green%27s_identities).


__Teorem.__ Rješenje problema rubnih vrijednosti

$$
\begin{aligned}
& u_t-D\Delta u=f, \quad x\in V,\quad t>0,\\
&u(x,0)=g(x),\quad  x\in V, \\
&u(x,t)=h(x,t),\quad x\in \partial V,\quad t>0,
\end{aligned}$$

pri čemu su funkcije $f$, $g$ i $h$ neprekidne, je jedinstveno.

__Dokaz__: Dokaz je kontradikcijom koristeći integral energije. 
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

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "49e3a36e9dabc387a7878ff958c2158d7abd723c"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "f6532909bf3d40b308a0f360b6a0e626c0e263a8"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.1"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "98f59ff3639b3d9485a03a72f3ab35bab9465720"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.6"

[[deps.PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═20dc641a-1cd4-4806-92e4-0fb85767754e
# ╠═31e0ee64-1b7c-4f26-beae-b0ce127aaf72
# ╟─a1a4c928-af9a-40bb-9d2b-9458461e364f
# ╟─3d4ac4e3-3bfa-4c19-834a-f91db0490ad4
# ╟─31f19a11-c08e-4d5a-a702-1621ba693e18
# ╟─70259200-0fc2-42b8-a471-d941347adb2f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
