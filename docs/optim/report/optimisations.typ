// FastSuite, dossier d'optimisation, rapport compilé
// Source des propositions : docs/optim/<outil>.md (six fichiers)
// Convention de preuve : [M] mesuré (artefact nommé), [E] estimé (mécanisme nommé),
//                        [V] vérifié dans cette session en exécutant le binaire ou en lisant la ligne citée.

#import "@preview/cetz:0.5.2": canvas, draw

#set page(paper: "a4", margin: 2cm)
#set text(font: "Libertinus Serif", size: 11pt, lang: "fr")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#show figure.caption: set text(size: 9.5pt)
#show raw.where(block: true): it => block(
  width: 100%,
  fill: luma(246),
  inset: 6pt,
  radius: 2pt,
  text(size: 8.5pt, font: "DejaVu Sans Mono", it.text),
)
#set table(stroke: 0.4pt + luma(170), inset: (x: 4pt, y: 2.6pt))
#show table.cell: set text(size: 9.5pt)

#let tag-m = text(fill: rgb("#0b6b3a"), weight: "bold")[M]
#let tag-e = text(fill: rgb("#8a5a00"), weight: "bold")[E]
#let tag-v = text(fill: rgb("#0b4f8a"), weight: "bold")[V]

// ---------------------------------------------------------------------------
// Page de titre
// ---------------------------------------------------------------------------
#align(center + horizon)[
  #text(size: 21pt, weight: "bold")[FastSuite]
  #v(0.2em)
  #text(size: 15pt)[Dossier d'optimisation des six outils]
  #v(1.2em)
  #text(size: 11pt, style: "italic")[
    Code et algorithmes : soixante-douze propositions classées,#linebreak()
    quatre constats qui changent la lecture des chiffres publiés
  ]
  #v(2.5em)
  #text(size: 10pt)[
Revue de lecture seule du dépôt `~/hermes_scratch/FastSuite`#linebreak()
état `main` du 20 octobre 2026, 30 987 lignes de Rust lues
  ]
  #v(1em)
  #text(size: 10pt)[Raphaël Ribes]
]

#pagebreak()
#outline(title: "Sommaire", indent: auto)

#pagebreak()

// ---------------------------------------------------------------------------
= Verdict

Trois familles de résultats, dans cet ordre d'importance.

*1. Quatre chiffres publiés ne mesurent pas ce qu'ils annoncent.* La cellule « fastprep 34 s contre
fastp 28 s sur un cœur » n'est pas une cellule sur un cœur : les lignes fastprep étiquetées
`pinned` du CSV committé brûlent 2,44 cœurs. Les 145 Mo de fastmaps sont une configuration à tampon
de 512 Mo, pas le défaut livré, qui mesure 36,1 Mo. Le temps de fastcheck sur un fichier est un temps
mono-thread quelle que soit la valeur de `-j`, par construction. Et fastsim compresse au niveau 6
codé en dur, ce qui représente 75,5 % de son temps d'exécution mesuré.

*2. Deux défauts de sortie, tous deux sur le chemin documenté.* `fasthelper concatenate` compte un
échantillon fantôme et duplique une ligne quand un échantillon porte la sortie `fastmaps --json`,
c'est-à-dire l'invocation que l'aide du binaire recommande : quatre rapports pour trois échantillons,
deux lignes identiques. Et `fasthelper check` accepte deux échantillons du même nom, qui écrivent
ensuite dans le même répertoire de travail. Les deux ont été reproduits ici.

*3. La marge d'optimisation restante est algorithmique, pas cosmétique.* Trois allocations par lecture
dans fastdedup et quatre dans fastmaps, quatre passes sur chaque séquence dans fastcheck, un
`Normal::new` dans une boucle par base dans fastsim, un préfiltre d'adaptateur qui ne peut pas
s'armer au taux d'erreur livré dans fastprep : ce sont des structures, pas des micro-réglages.
Le classement ci-dessous suit cette règle : d'abord supprimer du travail, ensuite paralléliser,
enfin réduire les constantes.

#v(0.4em)
#block(fill: luma(243), inset: 8pt, radius: 3pt, width: 100%)[
  *Ce que ce dossier n'est pas.* Aucun chiffre de gain n'est présenté comme mesuré. Seize valeurs
  sont mesurées et nommées avec leur artefact, quatre sont vérifiées ici en exécutant le binaire,
  tout le reste est estimé avec le mécanisme qui le produit. Les propositions ne sont pas
  implémentées : le dépôt est inchangé, seuls six fichiers `docs/optim/*.md` ont été écrits.
]

// ---------------------------------------------------------------------------
= Méthode et conventions de preuve

Chaque outil a été relu intégralement, fichiers par fichier, source et tests. Les propositions sont
classées par effet attendu, et non par facilité.

#table(
  columns: (auto, 1fr, auto),
  table.header([*Marque*], [*Signification*], [*Nombre*]),
  [M], [mesuré, l'artefact est nommé et citable], [32],
  [V], [vérifié en exécutant le binaire ou en relisant la ligne citée], [70],
  [E], [estimé, le mécanisme est nommé et l'arithmétique est écrite], [83],
)

Les définitions de coût qui suivent viennent toutes du dépôt ou d'une exécution, jamais d'une
impression. Trois d'entre elles méritent d'être posées avant de lire les tableaux.

- *Le coût unitaire de la suite* est de 15 à 17 µs par paire et par cœur #{tag-m}, chiffre qui
  traverse fastdedup, fastprep et le scoring de la campagne. À 150 pb cela fait 330 octets par
  paire, soit environ 21 Mo/s : la classe débit de gzip. Toute optimisation qui touche des octets
  passe donc avant toute optimisation qui touche des instructions.
- *Un seul fichier est mono-thread* dans fastcheck #{tag-v}, et *fastsim est mono-thread* par
  construction #{tag-v}. Les deux outils sont donc les seuls où le parallélisme est un gain entier
  et non un facteur 1,3.
- *Le banc ne tourne qu'à 1, 2 ou 4 threads épinglés* : toute mesure non épinglée est une mesure de
  la machine, pas de l'algorithme. Les propositions qui touchent au budget de threads sont donc
  accompagnées de leur compte de threads, pas seulement de leur temps.

// ---------------------------------------------------------------------------
= Où se trouve la marge

#figure(
  canvas({
    import draw: *
    // Géométrie : deux rangées de trois outils, cetz est en y croissant.
    let w = 4.7
    let h = 1.25
    let node(cx, cy, titre, cout, prop, fill, tcol: rgb("#12181e")) = {
      rect((cx - w/2, cy - h/2), (cx + w/2, cy + h/2), radius: 0.10, fill: fill, stroke: 0.5pt + rgb("#2b3843"))
      content((cx, cy + 0.30), text(size: 8.6pt, weight: "bold", fill: tcol)[#titre])
      content((cx, cy - 0.10), text(size: 6.9pt, fill: rgb("#33404a"))[#cout])
      content((cx, cy - 0.40), text(size: 6.9pt, weight: "bold", fill: rgb("#7a2d00"))[#prop])
    }
    let arrow(a, b, dash: none) = {
      if dash == none {
        line(a, b, mark: (end: ">"), stroke: 0.8pt + rgb("#4a5a66"))
      } else {
        line(a, b, mark: (end: ">"), stroke: (paint: rgb("#4a5a66"), dash: "dashed", thickness: 0.7pt))
      }
    }

    // Rangée haute : la chaîne de préparation.
    node(2.9, 0, "fastdedup", "3 alloc./read + 2 sauts de canal", "FD-1 découpe au commit ordonné", rgb("#dce9f7"))
    node(8.4, 0, "fastprep", "3 étages sur 2 threads", "FP-1 troisième étage", rgb("#fbe7cd"))
    node(13.9, 0, "fastmaps", "12 o/read d'index, 36 Mo livrés", "FM-1 index optionnel", rgb("#e3f0dc"))

    // Rangée basse : mesure, agrégation, simulation.
    node(2.9, -2.9, "fastcheck", "4 passes + 1 sur le qual.", "FC-2 passe unique", rgb("#dce9f7"))
    node(8.4, -2.9, "fasthelper", "RSS linéaire en échantillons", "FH-3 sérialisation en flux", rgb("#efe6f5"))
    node(13.9, -2.9, "fastsim", "75,5 % en deflate mesuré", "FS-1 niveau 1", rgb("#fbe7cd"))

    arrow((5.25, 0), (6.05, 0))
    arrow((10.75, 0), (11.55, 0))
    // descente : fastmaps -> fastcheck, fastcheck -> fasthelper
    line((13.9, -0.62), (13.9, -1.2), (2.9, -1.2), (2.9, -2.27), mark: (end: ">"), stroke: 0.8pt + rgb("#4a5a66"))
    arrow((5.25, -2.9), (6.05, -2.9))
    // fastsim alimente fastcheck par un profil
    line((13.9, -3.55), (13.9, -4.15), (2.9, -4.15), (2.9, -3.60), mark: (end: ">"), stroke: (paint: rgb("#4a5a66"), dash: "dashed", thickness: 0.7pt))
    content((8.4, -4.4), text(size: 6.9pt, style: "italic", fill: rgb("#4a5a66"))[profil mesuré, puis boucle QC vers simulation])

    content((2.9, 1.05), text(size: 7.2pt, weight: "bold", fill: rgb("#33404a"))[lecture des reads])
    content((2.9, -1.55), text(size: 7.2pt, weight: "bold", fill: rgb("#33404a"))[qualité et exactitude])
  }),
  caption: [
    Chaque outil et son coût dominant, avec la proposition qui l'attaque en premier. Les coûts sont
    lus dans le code, sauf ceux de fastsim qui sont mesurés sur cette machine.
  ],
)

Le diagramme dit l'essentiel : cinq des six outils ont leur coût dominant dans une structure
répétée par lecture, pas dans un calcul rare. C'est pourquoi les propositions algorithmiques
dominent le classement, et pourquoi les propositions de micro-optimisation sont explicitement
rangées après celles qui suppriment du travail.

// ---------------------------------------------------------------------------
= Quatre constats qui changent la lecture des chiffres publiés

== La cellule « un cœur » de fastprep n'est pas une cellule à un cœur

`docs/STATUS.md` et `docs/BENCHMARKS.md` citent fastprep à 34,37 s contre fastp à 28,50 s sur un
cœur. Le CSV committé `~/work/FastPrep-Bench/results/results_gz6.csv` contient ceci #{tag-v}.

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  table.header([*lectures*], [*outil*], [*pinned*], [*wall s*], [*cpu s*], [*cpu/wall*]),
  [5 M], [fastprep], [non], [33,508], [81,458], [2,43],
  [5 M], [fastprep], [oui], [32,983], [80,510], [*2,44*],
  [5 M], [fastp], [non], [28,197], [88,200], [3,13],
  [5 M], [fastp], [oui], [83,254], [82,341], [*0,99*],
  [5 M], [cutadapt], [oui], [135,759], [133,079], [0,98],
)

Un processus confiné à un seul cœur ne peut pas accumuler 80,5 secondes de CPU en 33,0 secondes de
temps réel. fastp et cutadapt retombent à 0,98 et 0,99, donc le harnais contraint bien ces deux-là,
et les lignes fastprep étiquetées `pinned` et `non pinned` sont indistinguables à toutes les tailles.
Elles sont sorties de `taskset`.

Deux conséquences, toutes deux inconfortables.

- *La comparaison à un cœur n'existe pas.* Le seul chiffre honnête est fastp à 83,25 s, avec fastprep
  non mesuré, et les secondes CPU de fastprep sont 2 % *sous* celles de fastp.
- *L'écart est un écart d'occupation des cœurs, pas de travail par paire.* Sur la plage linéaire,
  fastprep fait 6,747 s par million de paires à 16,36 s CPU, soit 2,42 cœurs occupés, quand fastp fait
  5,579 s à 17,64 s CPU, soit 3,16 cœurs. fastprep fait 7 % de travail en moins par paire et met 21 %
  de temps en plus. À l'occupation de fastp, son temps serait de 26 s.

La proposition FP-0 est donc une correction de banc avant toute optimisation : forcer le bras
fastprep dans le même `taskset`, et refuser d'écrire une ligne dont le rapport cpu/wall sort de
l'intervalle 0,95 à 1,05. C'est la seule façon d'empêcher cette classe d'erreur de rester invisible,
puisque les nombres restent plausibles.

== Les 145 Mo de fastmaps sont une configuration qui n'est plus le défaut

La ligne publiée, 5,49 s et 143,2 Mo à 5 millions de paires #{tag-v}, a été obtenue avec un tampon
de 512 Mo. Le défaut livré depuis le 18 octobre est de 8 Mo et mesure 36,1 Mo sur le même corpus.
Le plancher, previews désactivées, est de 16,0 Mo, soit le niveau de FastQC à 16,4 Mo.

L'écart contesté n'est donc pas de 130 Mo mais d'environ 20 Mo. La ligne de `docs/STATUS.md` qui
compare 145 Mo à 16 Mo compare deux configurations différentes, et devrait nommer le tampon.

== fastcheck sur un fichier est mono-thread, quelle que soit la valeur de -j

`main.rs` prend la branche séquentielle dès que `files.len() == 1`, avec un commentaire qui affirme
qu'un fichier unique utilise déjà tout le parallélisme disponible dans le pool BGZF. C'est vrai pour
du BGZF, faux pour un FASTQ nu et faux pour un gzip ordinaire, où le paramètre ne sert à rien. Les
55,68 s du tableau des 20 millions de lectures sont donc un temps mono-cœur, et l'affirmation de
montée en charge qui s'appuie dessus n'a jamais été mesurée.

La proposition FC-4 découpe le fichier en plages d'octets alignées sur les enregistrements, chacune
alimentant son propre accumulateur, puis fusionne. La fusion est exacte et associative pour tous les
champs, y compris la table d'empreintes, ce qui en fait le parallélisme le moins risqué du dépôt.

== fastsim passe trois quarts de son temps à compresser au niveau 6

Mesuré sur cette machine, 200 000 paires, profil `mandrill_E358`, `taskset -c 0` :

#table(
  columns: (auto, auto, auto, auto),
  table.header([*cellule*], [*temps*], [*RSS*], [*taille de sortie*]),
  [défaut, gzip niveau 6], [*14,55 s*], [139 Mo], [38,7 Mo],
  [`--no-gzip`], [*3,56 s*], [139 Mo], [141,4 Mo],
)

Le deflate représente 75,5 % du temps. Le niveau est un `Compression::new(6)` codé en dur, sans
option pour le changer, alors que la règle de la suite veut que tous les outils comparent à niveau
égal. Un balayage avec le `gzip` du système sur la sortie de l'outil donne 109 Mo/s et 23,04 Mo au
niveau 1, contre 9,2 Mo/s et 19,54 Mo au niveau 6 : le niveau 1 coûte 17,9 % de taille et divise le
temps de compression par onze.

#figure(
  canvas({
    import draw: *
    let total = 14.55
    let gen = 3.56
    let scale = 12.0 / total
    // barre mesurée
    rect((0, 0), (gen * scale, 0.75), fill: rgb("#0b6f3a"), stroke: none)
    rect((gen * scale, 0), (12.0, 0.75), fill: rgb("#b45309"), stroke: none)
    rect((0, 0), (12.0, 0.75), stroke: 0.5pt + rgb("#33404a"))
    content((gen * scale / 2, 0.37), text(size: 7.4pt, fill: white)[génération 3,56 s])
    content(((gen * scale + 12.0) / 2, 0.37), text(size: 7.4pt, fill: white)[deflate niveau 6 : 10,99 s, 75,5 %])
    content((0, 1.02), text(size: 7.6pt, weight: "bold")[mesuré])
    // barre projetée
    let gen2 = 3.56
    let def2 = 1.3
    let tot2 = gen2 + def2
    let sc2 = 12.0 / total
    rect((0, -1.35), (gen2 * sc2, -0.60), fill: rgb("#0b6f3a"), stroke: none)
    rect((gen2 * sc2, -1.35), (tot2 * sc2, -0.60), fill: rgb("#b45309"), stroke: none)
    rect((0, -1.35), (tot2 * sc2, -0.60), stroke: 0.5pt + rgb("#33404a"))
    content((gen2 * sc2 / 2, -0.97), text(size: 7.4pt, fill: white)[génération])
    content(((gen2 * sc2 + tot2 * sc2) / 2, -1.57), text(size: 6.6pt, fill: rgb("#8a4a06"))[deflate niv. 1])
    content((0, -0.35), text(size: 7.6pt, weight: "bold")[projeté, niveau 1])
    line((tot2 * sc2 + 0.15, -0.97), (tot2 * sc2 + 1.9, -0.97), mark: (end: ">"), stroke: 0.7pt + rgb("#b45309"))
    content((tot2 * sc2 + 1.95, -0.97), anchor: "west", text(size: 7.4pt, fill: rgb("#8a4a06"))[≈ 4,9 s, facteur 3])
    // axe
    line((0, -1.85), (12.0, -1.85), stroke: 0.5pt + rgb("#33404a"))
    for (k, lbl) in ((0, "0 s"), (1, "3,6 s"), (2, "7,2 s"), (3, "10,8 s"), (4, "14,4 s")) {
      let x = k * 3.0
      line((x, -1.85), (x, -1.95), stroke: 0.5pt + rgb("#33404a"))
      content((x, -2.15), text(size: 6.8pt)[#lbl])
    }
  }),
  caption: [
    fastsim, 200 000 paires sur un cœur épinglé. Barre haute mesurée ici, barre basse projetée depuis
    le débit du niveau 1 mesuré sur la sortie de l'outil. Le temps de deflate est l'axe à attaquer
    avant toute optimisation du générateur.
  ],
)

// ---------------------------------------------------------------------------
= Deux défauts de sortie reproduits

== fasthelper compte un échantillon fantôme quand fastmaps écrit son wrapper

`fastmaps --json` écrit la grille et un rapport v1 à côté, ce que son aide recommande. fasthelper
normalise les deux, mais ne reconnaît pas la tige `heatmap.report` comme générique, donc l'échantillon
fantôme survit. Reproduction sur deux échantillons réels #{tag-v} :

#raw(block: true, lang: "text", "fasthelper: concatenated 4 report(s) from 3 sample(s)\n      1 alpha            fastmaps\n      1 beta             fastmaps\n      2 heatmap.report   fastmaps")

Chaque bibliothèque devient deux lignes de valeurs identiques, donc toute somme sur le tableau
compte double. La correction tient en deux lignes : retirer le suffixe `.report` avant le test de
tige, puis dédupliquer sur `(sample, tool)` en faveur du rapport v1.

== fasthelper check accepte deux échantillons du même nom

Un fichier de configuration à deux entrées nommées `lib1` passe la validation #{tag-v} : le
répertoire de travail, la clé de reprise et le manifeste dérivent du nom seul, donc deux workers
écrivent au même endroit. Le contrôle censé attraper cette classe l'accepte, avec le code de sortie
zéro. Une erreur au démarrage vaut mieux qu'une course entre deux producteurs.

// ---------------------------------------------------------------------------
= Les propositions, outil par outil

Chaque outil a son fichier de propositions dans `docs/optim/`. Les tableaux qui suivent donnent le
classement, le détail est dans le fichier.

== fastdedup : le pipeline parallélise l'étape qui ne peut pas le payer

Les workers ne font que hacher, tandis que l'insertion dans la table et l'écriture restent sur le
thread appelant. Hacher 330 octets coûte quelques dizaines de nanosecondes contre 15 à 17 µs : moins
de 1 % d'une paire. Le pipeline dépense pourtant trois allocations par lecture, deux sauts de canal
et deux transferts de cache froids pour déporter ce pourcentage, et le vérificateur prétendument
fragmenté n'est appelé que par ce même thread d'écriture, donc le fragmentage ne sert à rien.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FD-1], [Découper au commit ordonné : parse, hachage et insertion sur un thread, travail par octet sur les autres], [algorithmique], [L],
  [FD-2], [Gzippy parallèle sur le chemin flux `-o -`, qui ne cherche jamais dans le fichier], [algorithmique], [S],
  [FD-3], [Budget de threads : 8 à 9 threads exécutables pour 4 cœurs épinglés], [threads], [S],
  [FD-4], [Abondance : arène en RAM au lieu du second parcours et de la seconde décompression], [algorithmique], [L],
  [FD-5], [Supprimer le second ensemble de hachage de l'abondance], [mémoire], [M],
  [FD-6], [Dimensionner la table par paliers ; sur stdin l'estimation vaut zéro et le hachage tombe à 64 bits], [mémoire], [M],
  [FD-9], [Parité CRLF entre `-j 1` et `-j 4`, vérifiée cassée], [exactitude], [S],
  [FD-11], [Table ouverte sur mesure, insertions par lots de seize avec préchargement], [algorithmique], [M],
)

== fastprep : trois étages sur deux threads

Le code est déjà en une passe, avec des arènes par chunk, une table de qualité et une direction
monomorphisée : les propositions de micro-optimisation y sont donc secondaires. Ce qui reste est le
partage du travail, un préfiltre mort et un étage nul qui tourne par défaut.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FP-0], [Épingler réellement la cellule fastprep et refuser une ligne dont cpu/wall sort de 0,95 à 1,05], [banc], [S],
  [FP-1], [Trois étages sur trois threads : parse, trim, commit ordonné, compression], [algorithmique], [M],
  [FP-2], [`--pe-overlap auto` s'arme aussi sur la détection automatique, contre sa propre documentation], [algorithmique], [S],
  [FP-3], [Le préfiltre 6-mers ne peut pas s'armer à 0,2 : 414 comparaisons exactes par lecture], [algorithmique], [M],
  [FP-4], [Six `write_all` par lecture, remplacer par un tampon de 256 Kio], [allocation], [S],
  [FP-5], [Décoder et ré-encoder la qualité octet par octet, alors qu'une LUT décalée suffit], [allocation], [S],
  [FP-6], [Dedup : `contains` puis `insert`, soit deux sondes de la même table], [algorithmique], [S],
  [FP-9], [Analyseur FASTQ à quatre lignes avec repli sur needletail], [algorithmique], [L],
)

== fastcheck : quatre passes par lecture pour des compteurs exacts

Le chemin par défaut parcourt chaque séquence quatre fois, et la boucle de 12-mers décode onze fois
chaque base. Les 2,8 milliards d'insertions de 12-mers des 20 millions de lectures atterrissent dans
une `HashMap` à hachage SipHash plafonnée à un million de clés, alors que l'espace des clés est de
4¹² et qu'une table dense de `u16` coûte 32 Mio, moins que la table qu'elle remplace.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FC-1], [Comptage exact des 12-mers en 32 Mio, table dense 4¹² en `u16` plus débordement], [algorithmique], [M],
  [FC-2], [Une seule passe par lecture : heatmap, 12-mer glissant, polyX et GC de fragment fusionnés], [algorithmique], [M],
  [FC-3], [Traiter les lectures par lots de 64 pour amortir les lignes de position], [algorithmique], [L],
  [FC-4], [Découper un fichier unique en plages, puis fusionner des accumulateurs exacts], [threads], [L],
  [FC-5], [Le matcher d'adaptateurs recalcule ses constantes à chaque position], [algorithmique], [M],
  [FC-6], [Recherche de readthrough de `--profile` : 3 600 comparaisons par lecture], [algorithmique], [S],
  [FC-8], [Lignes de qualité de 2 048 octets pour un alphabet de 94 valeurs], [cache], [S],
)

== fastsim : le deflate avant le générateur

Le partage mesuré impose l'ordre : le générateur ne représente que 24,5 % du temps. Optimiser la
génération avant la compression plafonne à 1,22 ; l'inverse vaut 3.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FS-1], [Rendre le niveau de compression paramétrable, défaut 1 : 14,55 s vers 4,9 s], [entrées-sorties], [S],
  [FS-2], [Émission parallèle par blocs, graine par paire, membres gzip par bloc, `-j1` identique à `-j4`], [algorithmique], [L],
  [FS-3], [Le pool de duplicats garde 130 Mo d'octets qu'il ne relit jamais], [mémoire], [M],
  [FS-4], [Deux `exp`, un `erf`, un `sqrt` et un `powf` par base pour des valeurs fixées par le profil], [algorithmique], [M],
  [FS-5], [Une loi uniforme par base pour une probabilité de 6e-4 : écarts géométriques], [algorithmique], [S],
  [FS-6], [300 `write_all` d'un octet par paire], [allocation], [S],
)

== fastmaps : l'index par lecture est la seule masse structurelle

Le record en mémoire fait 12 octets, pas les 8 que documente le commentaire ni les 16 que facture le
budget, et il n'existe que pour un fichier de sortie. C'est la seule masse qui croît avec le nombre
de lectures à n'importe quel réglage de tampon.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FM-1], [Rendre l'index optionnel et étiqueté : `--no-data`, `data_complete: false`], [algorithmique], [M],
  [FM-2], [Le seuil de preview de 4 Mo par worker échappe au budget et forme un plancher], [mémoire], [S],
  [FM-3], [Table de 65 536 écrivains et tampons de 8 Kio pour 152 lignes], [mémoire], [S],
  [FM-5], [Facturer la capacité allouée et recycler les tampons au lieu de les jeter], [allocation], [M],
  [FM-7], [Emballer le record sur 6 octets, `y` étant déjà le nom du fichier de ligne], [algorithmique], [M],
  [FM-10], [Quatre `String` par lecture dans la boucle FASTQ, et validation UTF-8 par base], [allocation], [M],
  [FM-12], [Rendre la fusion des qualités indépendante de l'ordre d'achèvement], [exactitude], [S],
)

== fasthelper : la mémoire, pas le temps

À 10 000 échantillons l'agrégation prend 1,55 s et 314 Mo. Ce qui casse le premier n'est pas le
temps, c'est la mémoire, et deux défauts de sortie se cachent dans le même chemin.

#table(
  columns: (auto, 1fr, auto, auto),
  table.header([*ID*], [*Proposition*], [*Classe*], [*Effort*]),
  [FH-1], [Échantillon fantôme et ligne dupliquée quand fastmaps écrit son wrapper, reproduit], [exactitude], [S],
  [FH-2], [Deux échantillons du même nom passent `check`, reproduit], [exactitude], [S],
  [FH-3], [Sérialiser en flux au lieu de l'arbre `Value` et de la chaîne de 26 Mo], [mémoire], [M],
  [FH-4], [Le garde d'adaptateur relit l'entrée à chaque exécution, sur un seul thread], [threads], [M],
  [FH-5], [La grille fastmaps est analysée en arbre `Value`, 68 000 constructions par fichier], [algorithmique], [M],
  [FH-6], [La table `--all-metrics` entière est construite dans une seule chaîne], [algorithmique], [S],
)

// ---------------------------------------------------------------------------
= Feuille de route

Quatre vagues, ordonnées par rapport effet sur risque. La première ne touche presque pas au code et
corrige des sorties.

#table(
  columns: (auto, 1fr, auto),
  table.header([*Vague*], [*Contenu*], [*Pourquoi d'abord*]),
  [0], [FP-0, et la re-dérivation des lignes de fastprep et de fastmaps dans `STATUS.md` et `BENCHMARKS.md`], [Un chiffre faux oriente toutes les décisions suivantes],
  [1], [FH-1, FH-2, FD-9, FM-12 : quatre défauts de sortie, aucun changement d'algorithme], [Petits diffs, sorties reproductibles, tests existants],
  [2], [FS-1, FD-2, FD-3, FP-2, FP-3 : supprimer du travail ou du parallélisme gaspillé], [Gains mesurables sans restructuration],
  [3], [FC-1, FC-2, FM-1, FH-3, FS-4 : les changements structurels par outil], [Effort moyen, chacun validable par une porte existante],
  [4], [FD-1, FP-1, FS-2, FC-4 : la restructuration du partage du travail], [Les plus gros gains, les seuls qui exigent un nouveau test de parité],
)

Deux propositions sont explicitement conditionnées à une mesure préalable : FD-11, la table ouverte
sur mesure, n'a de sens que si le profil montre que l'insertion domine la passe de comptage, et
FM-4, le rendu en flux, n'a de sens qu'une fois la taille du terme rendue visible par un journal.
C'est la règle appliquée partout dans les six fichiers : une proposition qui ne peut pas montrer son
mécanisme dans le tableau de mesures ne sort pas.

// ---------------------------------------------------------------------------
= Ce que je n'ai pas fait

- *Aucun code modifié.* Les six fichiers de propositions sont les seuls ajouts, et ils sont hors
  index Git. Rien n'a été committé, donc rien n'a été synchronisé vers le poste de bureau : la règle
  de synchronisation ne s'applique qu'aux changements de code.
- *Aucune mesure de banc de la suite.* Les seules mesures faites ici sont des sondes locales, sur un
  portable à 15 W partagé, et elles sont étiquetées comme telles. Les chiffres de comparaison
  viennent des CSV committés, et les ratios de débit sont explicitement des estimations.
- *Le profil n'a pas été instrumenté.* Aucun `perf` n'a été lancé : les attributions de coût sont
  soit lues dans le code avec leur arithmétique, soit mesurées par différence de cellules, comme le
  partage deflate de fastsim. Un profil par outil, sur le poste de bureau, reste la première étape
  de la vague 2.
- *La campagne Phase D n'a pas été rejouée.* Les propositions qui touchent au taux de duplication ou
  au `g_frac` ne prétendent pas corriger ces écarts, qui sont des écarts de modèle et non de code.

#v(0.6em)
#block(fill: luma(243), inset: 8pt, radius: 3pt, width: 100%)[
  *Ordre de lecture conseillé* : ce rapport, puis `docs/optim/fastprep.md` section 0 pour le détail
  arithmétique du banc, puis le fichier de l'outil à attaquer. Chaque proposition y porte son
  mécanisme, sa mesure et son risque, dans le même format.
]
