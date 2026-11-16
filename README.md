# Filogenia gronkowców na podstawie MLST – potok Snakemake

Potok [Snakemake](https://snakemake.readthedocs.io/) typizujący genomy gronkowców metodą MLST (*Multilocus Sequence Typing*) i budujący na tej podstawie drzewo filogenetyczne metodą łączenia sąsiadów (*neighbour-joining*, NJ).

> Projekt wykonany w ramach przedmiotu *Potoki obliczeniowe w Bioinformatyce*, 2026/27.
> Autorzy: Agata Zając, Karina Dal, Hanna Milnikel

## Spis treści

1. [Opis działania](#opis-działania)
2. [Struktura repozytorium](#struktura-repozytorium)
3. [Wymagania i środowiska](#wymagania-i-środowiska)
4. [Dane wejściowe](#dane-wejściowe)
5. [Konfiguracja](#konfiguracja)
6. [Uruchomienie](#uruchomienie)
7. [Wyniki](#wyniki)
8. [Opis reguł](#opis-reguł)

## Opis działania

Dla każdego genomu (asemblacji) z katalogu wskazanego w konfiguracji potok:

1. wyszukuje programem `blastn` fragmenty genów MLST (sekwencje zapytań) w kontigach genomu,
2. zachowuje wyłącznie trafienia obejmujące 100% sekwencji zapytania, bez przerw w dopasowaniu (pełnej długości kopie alleli),
3. łączy trafienia dla wszystkich loci danego genomu,
4. określa numery alleli i typ sekwencyjny (ST) genomu,
5. łączy wyniki wszystkich genomów w jedną tabelę,
6. wybiera po jednej (połączonej) sekwencji loci dla każdego ST,
7. oblicza drzewo NJ i je wizualizuje.

```
blastn → filter → merge_flt → mlst → merge_mlst → unique → tree
 (locus × genom)              (genom)    (wszystkie genomy)
```

## Struktura repozytorium

```
.
├── Snakefile              # definicja potoku (reguły, zależności)
├── config.yml             # konfiguracja (m.in. katalog z genomami)
├── envs/                  # opisy środowisk conda
│   ├── powb-snake.yml     #   Snakemake
│   ├── powb-blast.yml     #   BLAST+
│   └── powb-python.yml    #   Python i biblioteki dla skryptów
├── scripts/               # skrypty Pythona wywoływane przez reguły
│   ├── filter.py          #   filtrowanie wyników BLAST
│   ├── mlst.py            #   przypisanie alleli i typu sekwencyjnego (ST)
│   ├── unique.py          #   wybór unikatowych sekwencji (po jednej na ST)
│   └── tree.py            #   budowa i rysowanie drzewa NJ
├── input/
│   ├── queries/           #   sekwencje zapytań: {seqid}.fna (po jednym locus)
│   ├── alleles.fna        #   baza alleli MLST (FASTA)
│   └── bigsdb.tsv         #   schemat MLST: profile alleli i odpowiadające im ST
├── genomesdir/            # genomy do analizy: {asmid}.fna.gz (nazwa katalogu = klucz genomesdir w config.yml)
├── output/                # wyniki (tworzone przez potok, ignorowane przez git)
├── log/                   # logi poszczególnych reguł (ignorowane przez git)
├── .gitignore
└── README.md
```

## Wymagania i środowiska

Wymagany jest menedżer pakietów [conda](https://docs.conda.io/) (lub mamba/miniforge). Potok korzysta z trzech środowisk:

| Środowisko    | Plik                  | Zawartość                                         | Używane przez                  |
|---------------|-----------------------|---------------------------------------------------|--------------------------------|
| `powb-snake`  | `envs/powb-snake.yml` | Snakemake 9.27                                    | uruchamianie potoku            |
| `powb-blast`  | `envs/powb-blast.yml` | BLAST+ 2.17                                       | reguła `blastn`                |
| `powb-python` | `envs/powb-python.yml`| Python 3.15, polars, biopython, toytree, packaging| reguły `filter`, `mlst`, `unique`, `tree` |

Ręcznie trzeba utworzyć tylko środowisko `powb-snake`; pozostałe Snakemake tworzy automatycznie (opcja `--use-conda`).

```bash
conda env create -f envs/powb-snake.yml
conda activate powb-snake
```

## Dane wejściowe

Dane wejściowe nie są częścią repozytorium. Należy je pobrać z serwera `celsius`
do głównego katalogu projektu. Potrzebne są dwa katalogi:

- `input/` – sekwencje zapytań i schemat MLST,
- katalog z genomami do analizy, o nazwie zgodnej z kluczem `genomesdir` w `config.yml`.

    scp -r <użytkownik>@celsius:<ścieżka_do_input> input
    scp -r <użytkownik>@celsius:<ścieżka_do_genomów> genomesdir

Po pobraniu struktura powinna wyglądać następująco:

    input/
    ├── queries/          # sekwencje zapytań {seqid}.fna
    ├── alleles.fna       # baza alleli MLST
    └── bigsdb.tsv        # schemat MLST
    genomesdir/           # genomy do analizy {asmid}.fna.gz

| Ścieżka                      | Opis                                                                      |
|------------------------------|---------------------------------------------------------------------------|
| `input/queries/{seqid}.fna`  | po jednym pliku FASTA na locus MLST (fragment genu, np. `gmk`)            |
| `<genomesdir>/{asmid}.fna.gz`| złożenia genomów gronkowców (FASTA skompresowany gzip); `{asmid}` to numer dostępowy asemblacji |
| `input/alleles.fna`          | baza alleli MLST (FASTA)                                                  |
| `input/bigsdb.tsv`           | schemat MLST (profile alleli → ST) w formacie tabelarycznym               |

Nazwy zapytań (`{seqid}`) i asemblacji (`{asmid}`) są wykrywane automatycznie na podstawie nazw plików (`glob_wildcards`), więc dodanie nowego genomu lub locus nie wymaga zmian w kodzie.


## Konfiguracja

Plik `config.yml` zawiera parametry potoku. Wykorzystywany jest klucz:

| Klucz       | Opis                                               | Przykład            |
|-------------|----------------------------------------------------|---------------------|
| `genomesdir`| ścieżka do katalogu z plikami `{asmid}.fna.gz`     | `genomesdir`        |

Pozostałe parametry (próg e-value `0.01`, minimalne pokrycie zapytania `100`%) są zapisane w sekcji `params` odpowiednich reguł w `Snakefile`.

## Uruchomienie

Z katalogu głównego repozytorium, po aktywacji środowiska `powb-snake`:

```bash
# podgląd planu wykonania (bez uruchamiania)
snakemake -n

# uruchomienie potoku na 4 rdzeniach (środowiska reguł tworzone automatycznie)
snakemake --use-conda --cores 4

# opcjonalnie: wykres grafu zależności (wymaga Graphviz)
snakemake --dag | dot -Tpng > dag.png
```

Skrypty z katalogu `scripts/` są wywoływane bezpośrednio, więc muszą mieć uprawnienie do wykonywania (`chmod +x scripts/*.py`) i nagłówek `#!/usr/bin/env python`.

## Wyniki

| Plik                                  | Opis                                                          |
|---------------------------------------|---------------------------------------------------------------|
| `output/tree.newick`                  | drzewo NJ w formacie Newick                                   |
| `output/tree.png`                     | wizualizacja drzewa                                           |
| `output/merge_mlst.sh`                | zbiorcza tabela MLST (plik tabelaryczny TSV mimo rozszerzenia `.sh`) |
| `output/unique.fna`                   | sekwencje loci, po jednej na ST (FASTA)                       |
| `output/blastn/`, `output/filter/`, `output/merge_flt/`, `output/mlts/` | wyniki pośrednie (na genom / locus) |
| `log/`                                | logi każdej reguły, osobno dla każdego genomu i locus         |

## Opis reguł

| Reguła       | Wejście                                   | Wyjście                                    | Opis |
|--------------|-------------------------------------------|--------------------------------------------|------|
| `blastn`     | zapytanie `{seqid}`, genom `{asmid}`      | `output/blastn/blastn_{asmid}_{seqid}.tsv` | Wyszukiwanie `blastn` (e-value ≤ 0.01); genom jest dekompresowany i przekazywany strumieniem (`-subject -`). |
| `filter`     | wynik `blastn`                            | `output/filter/filtered_{asmid}_{seqid}.tsv` | Zachowanie trafień z 100% pokrycia zapytania i długością równą długości zapytania (dopasowanie bez przerw). |
| `merge_flt`  | wyniki `filter` dla wszystkich loci       | `output/merge_flt/merged_{asmid}.tsv`      | Połączenie przefiltrowanych trafień wszystkich loci jednego genomu. |
| `mlst`       | `merge_flt`, `alleles.fna`, `bigsdb.tsv`  | `output/mlts/mlst_{asmid}.tsv`             | Przypisanie alleli i typu sekwencyjnego (ST) dla jednego genomu. |
| `merge_mlst` | wyniki `mlst` wszystkich genomów          | `output/merge_mlst.sh`                     | Jedna tabela MLST dla wszystkich genomów. |
| `unique`     | `merge_mlst`                              | `output/unique.fna`                        | Po jednej połączonej sekwencji loci dla każdego ST (FASTA). |
| `tree`       | `unique`, `merge_mlst`                    | `output/tree.newick`, `output/tree.png`    | Obliczenie drzewa NJ i jego wizualizacja. |
| `all`        | –                                         | `output/tree.newick`, `output/tree.png`    | Reguła docelowa: definiuje końcowe wyniki potoku. |
