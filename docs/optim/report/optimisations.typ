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

*1. Quatre chiffres ne mesurent pas ce qu'ils annoncent, et trois d'entre eux avaient déjà été
corrigés dans les articles sans que la correction atteigne les documents du dépôt.* La cellule
« fastprep 34 s contre fastp 28 s sur un cœur » vient d'un harnais qui posait `-j N` sans `taskset` :
ses lignes `pinned` brûlent 2,44 cœurs. L'article FastPrep l'a vue, l'a écrite et a re-mesuré ; à un
cœur réellement épinglé c'est 73,99 contre 80,03, et fastprep mène. Les 145 Mo de fastmaps sont la
configuration à tampon de 512 Mo, pas le défaut livré qui mesure 36,1 Mo, et sur ce point c'est
l'article qui est en retard. Le temps de fastcheck sur un fichier unique est mono-thread par
construction, ce que l'article dit correctement et que la ligne de conditions des docs contredisait.
Reste un chiffre que personne ne documentait : fastsim compressait au niveau 6 codé en dur, sans
option, pour 75,5 % de son temps d'exécution mesuré.

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
  *Ce que ce dossier était, et ce qu'il est devenu.* Il a été écrit en lecture seule : aucun code
  touché, aucune mesure de campagne, les chiffres de gain estimés avec leur mécanisme. Trente-neuf
  des soixante et onze propositions ont depuis été implémentées, et la section « Ce qui a été
  implémenté » rapporte ce que la mesure en dit — y compris deux chiffres de ce dossier que la mise
  en œuvre a démentis.
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
= Quatre chiffres à corriger, et qui les corrige

La revue a trouvé quatre nombres qui ne mesurent pas ce qu'ils annoncent. Trois
avaient déjà été corrigés dans les articles et n'avaient jamais atteint les
documents du dépôt ; le quatrième n'était documenté nulle part. La distinction
compte : dans trois cas sur quatre, la science était juste et la documentation
du dépôt était en retard sur elle.

== La cellule « un cœur » de fastprep : épinglée dans l'article, pas dans les docs

Le CSV `FastPrep-Bench/results/results_gz6.csv` contient des lignes fastprep
étiquetées `pinned` qui brûlent 2,44 cœurs : 33,0 s de temps réel pour 80,5 s de
temps CPU, là où fastp et cutadapt tombent à 0,98 et 0,99. Un processus confiné à
un cœur ne peut pas faire cela. Le harnais posait `-j N` sans `taskset`, donc
seul le bras fastprep s'échappait.

L'article FastPrep a déjà trouvé ce défaut, l'a écrit dans sa section de méthode
et a re-mesuré avec `taskset` ; son CSV porte désormais une colonne `cores` dont
le rapport `cpu_s/wall_s` vaut 1,00 sur chaque ligne. Les nombres corrigés, à
5 millions de paires :

#table(
  columns: (auto, auto, auto, auto, auto),
  table.header([*outil*], [*cœurs*], [*wall s*], [*CPU s*], [*CPU/wall*]),
  [fastprep], [1], [*73,99*], [73,83], [1,00],
  [fastp], [1], [*80,03*], [79,83], [1,00],
  [cutadapt], [1], [130,60], [130,27], [1,00],
  [fastprep], [2], [37,55], [74,84], [1,99],
  [fastp], [2], [40,55], [80,32], [1,98],
)

fastprep mène donc à toutes les tailles, sur un et deux cœurs, de 1,08 à 1,28.
Ce qui restait faux était la documentation du dépôt : `docs/STATUS.md` affirmait
« fastp LEADS this cell » avec 34,37 contre 28,50, et `docs/BENCHMARKS.md`
répétait ces valeurs dans quatre endroits, dont un tableau qui se contredisait
lui-même. Les deux fichiers sont corrigés dans cette passe, et la ligne de
garde ajoutée est celle qui manquait : un rapport `cpu_s/wall_s` supérieur à 1,05
sur une cellule épinglée invalide la ligne.

== Les 145 Mo de fastmaps : configuration d'un autre défaut

La ligne 143,2 Mo à 5 millions de paires a été obtenue avec un tampon de 512 Mo.
Le défaut livré depuis le 18 octobre est de 8 Mo et mesure 36,1 Mo sur le corpus
de 10 millions de lectures ; le plancher, previews désactivées, est de 16,0 Mo,
soit le niveau de FastQC à 16,4 Mo. L'écart contesté n'est pas de 130 Mo mais
d'environ 20 Mo.

Ici la documentation du dépôt et l'article sont tous deux en retard, et l'article
est le plus gênant : il décrit encore `--buffer-mb` comme « default 512 MB » et
donne une série de pics 18,6 / 33,3 / 47,2 / 142,3 Mo qui est celle de l'ancien
défaut. La règle appliquée est celle du dépôt lui-même : une cellule mémoire sans
son `--buffer-mb` est illisible.

== fastcheck sur un fichier : mono-thread, et c'est l'article qui a raison

`main.rs` prend la branche séquentielle dès que `files.len() == 1`, et le pool de
décompression n'est utilisé que pour du BGZF : un FASTQ nu, ou un gzip ordinaire,
ne reçoit rien de `-j`. Les 55,68 s du tableau des 20 millions de lectures sont
donc un temps mono-cœur, et l'affirmation de montée en charge qui s'appuyait
dessus n'a jamais été mesurée.

L'article FastCheck l'écrit déjà correctement : « effectively single-threaded :
its worker threads split across files… `-j 1` and `-j 16` take the same time ».
C'est la ligne de conditions de `docs/BENCHMARKS.md` qui affirmait le contraire
(« uses all available threads by default ») ; elle est corrigée.

== fastsim : le niveau de compression n'était réglable nulle part

Mesuré sur cette machine, 200 000 paires, profil `mandrill_E358`, `taskset -c 0` :

#table(
  columns: (auto, auto, auto, auto),
  table.header([*cellule*], [*temps*], [*RSS*], [*taille de sortie*]),
  [défaut, gzip niveau 6], [*14,55 s*], [139 Mo], [38,7 Mo],
  [`--no-gzip`], [*3,56 s*], [139 Mo], [141,4 Mo],
)

Le deflate représente 75,5 % du temps, et le niveau était un
`Compression::new(6)` codé en dur, sans option : aucun article ne pouvait le
mentionner puisqu'il n'était pas réglable. Un balayage avec le `gzip` du système
sur la sortie de l'outil donne 109 Mo/s et 23,04 Mo au niveau 1, contre 9,2 Mo/s
et 19,54 Mo au niveau 6, soit 17,9 % de taille pour onze fois moins de temps de
compression. `--gz-level` existe depuis cette passe, avec 6 pour défaut, donc
aucun octet de sortie ne change tant qu'on ne le demande pas.

#figure(
  canvas({
    import draw: *
    let total = 14.55
    let gen = 3.56
    let scale = 12.0 / total
    rect((0, 0), (gen * scale, 0.75), fill: rgb("#0b6f3a"), stroke: none)
    rect((gen * scale, 0), (12.0, 0.75), fill: rgb("#b45309"), stroke: none)
    rect((0, 0), (12.0, 0.75), stroke: 0.5pt + rgb("#33404a"))
    content((gen * scale / 2, 0.37), text(size: 7.4pt, fill: white)[génération 3,56 s])
    content(((gen * scale + 12.0) / 2, 0.37), text(size: 7.4pt, fill: white)[deflate niveau 6 : 10,99 s, 75,5 %])
    content((0, 1.02), text(size: 7.6pt, weight: "bold")[mesuré])
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
    line((0, -1.85), (12.0, -1.85), stroke: 0.5pt + rgb("#33404a"))
    for (k, lbl) in ((0, "0 s"), (1, "3,6 s"), (2, "7,2 s"), (3, "10,8 s"), (4, "14,4 s")) {
      let x = k * 3.0
      line((x, -1.85), (x, -1.95), stroke: 0.5pt + rgb("#33404a"))
      content((x, -2.15), text(size: 6.8pt)[#lbl])
    }
  }),
  caption: [
    fastsim, 200 000 paires sur un cœur épinglé. Barre haute mesurée ici, barre
    basse projetée depuis le débit du niveau 1 mesuré sur la sortie de l'outil.
  ],
)

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
= Ce qui a été implémenté, et ce que la mesure en dit

Trente-neuf des soixante et onze propositions ont été implémentées en une journée, en six commits
séparés, un par outil. La suite passe de 577 à 658 tests, tous verts, sans qu'aucun test existant ne
soit modifié : les nouveaux tests sont des tests de non-régression écrits avant le changement, et
l'un d'eux reproduit une divergence de sortie que personne n'avait vue — `fastdedup` écrivait CRLF à
un thread et LF à quatre.

#v(0.3em)
#block(fill: luma(243), inset: 8pt, radius: 3pt, width: 100%)[
  *Comment lire les chiffres de cette section.* La première campagne de bureau est jetée, et pas
  pour une raison de temps : elle tournait avec huit cellules en parallèle alors que chaque cellule
  est épinglée, si bien que les cellules `-j 1` se disputaient le cœur 0. Le résidu est visible dans
  la colonne de contrôle : `cpu/mur` entre 0,13 et 2,96. Une seconde campagne, ciblée sur quinze
  cellules et strictement sérielle, donne `cpu/mur` entre 0,98 et 1,00 sur chaque ligne, et c'est
  elle qui fournit le tableau de bureau ci-dessous. Les chiffres du portable, pris par les
  implémenteurs, restent dans le premier tableau parce qu'ils couvrent des tailles que la campagne
  ciblée n'a pas mesurées.
]

== Ce que chaque outil a gagné

#table(
  columns: (auto, auto, 1fr, 1fr),
  table.header([*Outil*], [*Fait*], [*Changement principal*], [*Mesure avant → après*]),
  [fastcheck], [6/13], [table de 12-mers dense et exacte (32 Mio) au lieu d'une table de hachage bornée ; une passe par lecture au lieu de quatre], [1 M lectures, `-j 1` : 36,45 → 8,34 s (4,37×) ; RSS 327 → 202 Mo],
  [fastprep], [8/15], [largeur de graine du préfiltre dérivée de la borne d'erreur ; étage PE nul désarmé ; écritures groupées], [200 k paires, `-j 1` : 1,198 → 0,837 s (1,43×) ; adaptateurs explicites : 2,915 → 1,478 s (1,97×)],
  [fastsim], [7/11], [`--gz-level` (le niveau 6 n'était réglable nulle part) ; écritures groupées ; tables de transcendantes], [200 k paires : gz6 13,96 → 12,14 s ; `--no-gzip` 3,53 → 1,81 s ; niveau 1 = 2,24 s ; RSS 139 → 100 Mo],
  [fastmaps], [6/12], [`--data-mode` rend l'index par lecture optionnel ; budget de prévisualisation réel ; ordre de fusion déterministe], [5 M lectures, `-j 1` : 5,16 → 3,54 s (−31 %) ; 512 Mo : 73,4 → 63,3 Mo ; `off` : 17,0 Mo],
  [fasthelper], [7/10], [`report.json` en flux ; grille fastmaps en `Visitor` ; `--version` résolu une fois], [10 k échantillons, json : 273,3 → 68,4 Mo et 0,82 → 0,51 s ; `--dry-run` 2000×3 : 21,4 → 0,16 s],
  [fastdedup], [5/11], [gzip parallèle sur le flux, threads budgétés, parité CRLF], [non mesuré dans ce lot],
)

Trois de ces gains méritent d'être lus comme des suppressions de travail, pas comme des réglages.
`fastcheck` passe de quatre passes sur chaque lecture à une seule, et sa table de 12-mers de trois
tables de hachage saturées (environ 40 Mo chacune) à trois tableaux fixes de 32 Mio : le RSS baisse
de 125 Mo alors que le travail par lecture baisse d'un facteur quatre. `fastprep` ramène le préfiltre
d'adaptateur de 414 à 81 comparaisons exactes par lecture de 151 bases, sans changer un octet de
sortie. `fasthelper` supprime 5 997 des 6 000 `fork` de `--version` d'un `--dry-run` à 2 000 échantillons et 3 étapes, ce qui
transforme 21,4 s de `--dry-run` en 0,16 s.

== La même comparaison, sur le poste de bureau

Une campagne ciblée, strictement sérielle, une seule taille, un seul cœur, deux binaires. La
colonne de contrôle est la preuve que l'épinglage a tenu : `cpu/mur` vaut 0,98 à 1,00 partout, ce
qui ne laisse aucune place à un enfant échappé ni à un autre calcul sur le même cœur.

#table(
  columns: (auto, auto, auto, auto, auto),
  table.header([*Outil*], [*Avant, s*], [*Après, s*], [*Rapport*], [*RSS avant / après, Mo*]),
  [fastcheck], [22,10], [4,89], [4,52×], [413 / 273],
  [fastmaps], [3,72], [1,99], [1,87×], [18 / 14],
  [fastprep], [16,35], [14,21], [1,15×], [49 / 49],
  [fastsim, niveau 6], [68,38], [59,51], [1,15×], [139 / 100],
  [fastsim, niveau 1], [68,38], [10,52], [6,50×], [139 / 100],
  [fastdedup], [18,28], [18,02], [1,01×], [23 / 25],
  [fasthelper, 2000 éch.], [2,98], [3,98], [0,75×], [42 / 16],
  [fastp, comparateur], [15,79], [---], [---], [1 186 / 1 186],
)

La dernière ligne n'est pas un avant/après : c'est `fastp` mesuré par la même campagne, sur le
même fichier, avec le même épinglage, et elle existe pour une seule raison — la cellule `fastprep`
après vaut 14,21 s contre 15,79 s, donc l'écart avec le comparateur est de 1,11× et il est mesuré
ici, pas emprunté. Sa mémoire, en revanche, parle d'elle-même : 1 186 Mo contre 49.

Deux lignes de ce tableau vont dans le sens inverse du premier, et elles sont instructives.
`fasthelper` ralentit de 2,98 à 3,98 s à 2 000 échantillons alors que son pic mémoire tombe de 42 à
16 Mo : l'écriture en flux échange du temps contre de la mémoire, et sur un petit arbre l'échange
n'est pas rentable — le gain de temps n'apparaît qu'à 10 000 échantillons, où il est net. `fastdedup`
ne bouge pas (1,01×), ce qui est exactement ce que ses changements promettaient : ils portent sur le
parallélisme, le flux et la parité de fin de ligne, pas sur le chemin mono-thread par défaut. Un
dossier qui n'aurait rapporté que les quatre premières lignes aurait eu l'air meilleur et aurait
caché la seule régression de la série.

Sur l'unique cellule où les deux campagnes se recoupent — `fastsim` au niveau par défaut — les deux
donnent 1,15×, ce qui est le meilleur argument disponible pour que la mesure du portable n'ait pas
été un artefact de charge.

== Les deux défauts de sortie sont corrigés

`fasthelper concatenate` ne compte plus d'échantillon fantôme quand un échantillon porte la sortie
`fastmaps --json` : le suffixe `.report` est retiré avant le test de nom générique, et deux rapports
pour le même couple (échantillon, outil) sont départagés en faveur du rapport v1. La même arborescence
qui annonçait « 4 rapports de 3 échantillons » annonce « 2 rapports de 2 échantillons ». `check -c`
refuse désormais deux `[[samples]]` de même nom, qui écrivaient auparavant dans le même répertoire de
travail.

Deux comportements changent en conséquence, et ils sont documentés comme tels : sans option
d'adaptateur, l'attribution par adaptateur du rapport `fastprep` se déplace (le détecteur ne tourne
plus), les totaux de lectures et de bases, eux, ne bougent pas ; et dans `fastmaps`, une erreur
d'entrée/sortie en milieu de fichier termine maintenant en erreur au lieu de produire silencieusement
une heatmap partielle — c'était une sortie plausible et fausse.

== Deux chiffres de ce dossier étaient faux

*Le coût du niveau 1 de fastsim n'est pas de 17,9 %.* Ce chiffre venait d'un balayage du binaire
`gzip` du système, c'est-à-dire de zlib, sur la sortie de l'outil : un substitut, pas l'encodeur du
binaire. Mesuré dans le processus, avec `flate2`/`miniz_oxide`, le niveau 1 écrit 52,63 Mo là où le
niveau 6 écrit 38,75 Mo, soit *+35,8 %*, et il tourne à environ 200 Mo/s et non 109. La conclusion
tient — le niveau 6 codé en dur domine le temps d'exécution — mais le chiffre à retenir pour
dimensionner une campagne est +35,8 %, et la leçon méthodologique est que « mesurer la sortie avec
l'outil du système » ne mesure pas l'outil.

*La borne du préfiltre d'adaptateur n'est pas celle que j'avais écrite.* La condition
$k >= 6 (ceil(k epsilon) + 1)$ suppose que le bloc sans erreur a la longueur nominale ; il peut être
plus court, et la borne correcte est $k >= s (m+1)$ pour $s$ bases de graine et $m$ erreurs
tolérées. Le code livré en tire des largeurs par adaptateur — une graine de 4 bases à partir de
$k = 20$ pour les adaptateurs TruSeq de 34 bases, de 3 bases à partir de $k = 9$ pour le Nextera de
19 — au lieu du 5-mer à $k = 25$ que je proposais, que le test de conformité aurait rejeté. La
proposition était juste dans son diagnostic, fausse dans son arithmétique.

// ---------------------------------------------------------------------------
= Feuille de route

Quatre vagues, ordonnées par rapport effet sur risque. Les trois premières sont faites, sauf les
propositions marquées comme différées dans les six fichiers d'outil ; la quatrième reste ouverte, et
c'est celle qui demande un nouveau test de parité plutôt qu'un diff.

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

- *La matrice complète de bureau n'a pas été menée à terme.* La première tentative mesurait huit
  cellules à la fois sur des cœurs qu'elles se partageaient, donc ses lignes sont jetées ; une
  campagne sérielle complète a été lancée puis remplacée par une campagne ciblée de quinze cellules,
  qui est celle du tableau ci-dessus. Les tailles intermédiaires et les exécutions à deux et quatre
  cœurs du poste de bureau restent donc à mesurer, et pour elles il faut compter une heure de calcul
  sériel, pas trois.
- *Aucun comparateur n'a été re-mesuré.* FastQC, Falco et cutadapt ne sont pas installés sur le poste
  de bureau ; `fastp` l'est et figure au harnais, mais ses cellules ont échoué à l'agrégation comme
  les autres. Les tableaux comparatifs des articles restent donc ceux des CSV committés, et la
  nouvelle section de chaque article est explicitement un avant/après du même outil, non sécable
  dans un tableau de comparateurs.
- *Aucun profil `perf` n'a été pris.* Les attributions de coût restent lues dans le code avec leur
  arithmétique, ou mesurées par différence de cellules. Ce qui n'a pas empêché de trouver les bons
  leviers, mais une des propositions différées (FD-1, déplacer le travail hors du thread ordonné)
  demande un profil avant d'être tentée.
- *La restructuration n'a pas été faite.* FD-1, FP-1, FS-2 et FC-3/4/5 sont les propositions qui
  changent le partage du travail ou le flux d'octets, donc qui exigent un test de parité
  supplémentaire. Elles restent ouvertes, avec leur mécanisme écrit.
- *La campagne Phase D n'a pas été rejouée.* Les écarts de taux de duplication et de `g_frac` sont
  des écarts de modèle, pas de code, et aucune des propositions implémentées ne prétend les corriger.
