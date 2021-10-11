# Matematika - posebna poglavlja

Bilježnica za predmet __[Matematika - posebna poglavlja](https://nastava.fesb.unist.hr/nastava/predmeti/14077)__ koji se predaje na diplomskim studijima Konstrukcijsko-enegetskog strojarstva (261) i Računalnog projektiranja i inženjerstva (262) na [FESB-u](https://www.fesb.unist.hr/).

Bilježnice su pisane u programskom jeziku [Julia](https://julialang.org) koristeći paket
[Pluto.jl](https://github.com/fonsp/Pluto.jl). Prilagođene su direktnoj i on-line nastavi.

# Pregledavanje bilježnica

Unutar svojeg preglednika, bilježnice možete pregledati na poveznici
[https://ivanslapnicar.github.io/Matematika-posebna_poglavlja/](https://ivanslapnicar.github.io/Matematika-posebna_poglavlja/)

#  Izvršavanje bilježnica

## Binder

1. Idite na adresu  https://ivanslapnicar.github.io/Matematika-posebna_poglavlja/ i odaberite željenu bilježnicu.
2. Pritisnite `Edit or run this notebook` i odaberite `binder`. Učitat će se svi poptrebni paketi i pokrenuti bilježnica (kroz nekoliko minuta).

## Računalo

1. Klonirajte cijeli repozitorij koristeći `git` naredbu:
```
git clone https://github.com/ivanslapnicar/Matematika-posebna_poglavlja.git
```
Ako niste upoznati s `git` alatom možete pogledati GitHubove [stranice za pomoć](https://help.github.com/articles/set-up-git/) ili direktno preuzeti bilježnice (repozitorij) kao zip datoteku. Možete koristiti i GitHub Dekstop.

2. Instalirajte [Julia-u](https://julialang.org/downloads/). U Julia terminalu izvedite naredbe
```
> using Pkg
> Pkg.add("Pluto")
```
Prethodne naredbe je potrebno izvršiti samo jednom.

3. Server Pluto bilježnica se pokreće pomoću naredbi
```
> using Pluto
> Pluto.run()
```

Sada možete izvršavati bilježnice koje se nalaze u direktoriju `Matematika-posebna_poglavlja/Pluto/`.
