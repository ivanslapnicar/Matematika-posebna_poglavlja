# Matematika - posebna poglavlja

Jupyter bilježnice predmeta _[Matematika - posebna poglavlja](https://nastava.fesb.unist.hr/nastava/predmeti/8237)_ koji se predaje na diplomskim studijima Konstrukcijsko-enegetskog strojarstva (261) i Računalnog projektiranja i inženjerstva (262) na [FESB-u](https://www.fesb.unist.hr/).

## Korištenje

Materijali su pisani kao [Jupyter](http://jupyter.org/) ([IJulia](https://github.com/JuliaLang/IJulia.jl)) bilježnice (engl. _notebooks_). Bilježnice možete koristiti na sljedeće načine.

### Korištenjem preglednika
Unutar svojeg preglednika bilježnice možete pregledati pomoću [Jupyter notebook viewera](http://nbviewer.jupyter.org/) na sljedećoj [poveznici](http://nbviewer.ipython.org/url/github.com/ivanslapnicar/Matematika-posebna_poglavlja/tree/master/src/).

###  Lokalno preuzimanje i pregled na vlastitom računalu

* Preuzmite bilježnice (repozitorij) korištenjem `git` naredbe:
```
git clone https://github.com/ivanslapnicar/Matematika-posebna_poglavlja.git
```
Ako niste upoznati s `git` alatom možete pogledati GitHubove [stranice za pomoć](https://help.github.com/articles/set-up-git/) ili direktno preuzeti bilježnice (repozitorij) kao zip datoteku.

* Instalirajte [Julia-u](https://julialang.org/downloads/). U Julia terminalu izvedite naredbe
```
> using Pkg
> Pkg.add("IJulia")
```
Prethodne naredbe je potrebno izvršiti samo jednom.

* Server bilježnica se pokreće u web pregledniku  pomoću naredbi
```
> using IJulia
> notebook(detached=true)
```
Sada možete izvršavati bilježnice koje se nalaze u direktoriju `Matematika-posebna_poglavlja/src`.
