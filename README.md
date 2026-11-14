# 3D Pore Boundary Benchmark — MATLAB

**Code source de l'étude comparative des sept définitions de frontière d'un espace poral 3D.**

Ce dépôt réunit **34 fichiers MATLAB `.m` d'origine**, retrouvés dans les fichiers du projet (travaux d'octobre–novembre 2026). Les programmes MATLAB sont conservés sans réécriture algorithmique. Un fichier archivé sous `analyser_tore_128_3_points(1).m` a simplement été renommé `analyser_tore_128_3_points.m` pour respecter le nom de sa fonction principale.

> **Statut :** archive de recherche / code expérimental. L'intégration a fait l'objet d'un contrôle statique des fichiers et des noms de fonctions, **pas d'une exécution MATLAB complète**. Les données d'étude ne sont pas fournies. Les différents `main_...` sont des variantes historiques de protocole et ne donnent pas nécessairement les mêmes résultats.

## Objet scientifique

Une image 3D binaire représente le pore (`1`, `true` ou `255`) et le solide (`0`). Sept ensembles de points frontières sont comparés :

| Méthode | Définition |
| --- | --- |
| `I6`, `I18`, `I26` | Centres des voxels **poreux** adjacents au solide sous le voisinage 6, 18 ou 26 |
| `E6`, `E18`, `E26` | Centres des voxels **solides** adjacents au pore sous le voisinage 6, 18 ou 26 |
| `F` | Centres des faces entre pore et solide |

**Convention de bord :** l'extérieur de l'image est assimilé au solide. Une couche virtuelle est ajoutée lors de l'extraction ; des points externes peuvent donc avoir une coordonnée `0` ou `dim+1`, et des centres de faces `0.5` ou `dim+0.5`.

Chaîne principale : extraction des frontières → tétraédrisation de Delaunay → exclusion des tétraèdres inadmissibles → calcul des sphères circonscrites → filtrage TD (débordement) → sélection TC (chevauchement) → reconstruction et évaluation géométrique.

## Organisation des fichiers

Tous les `.m` restent à la racine : plusieurs fonctions appellent leurs dépendances par leur nom MATLAB, sans chemin de package.

### Étapes fondamentales

- `extraire_7_frontieres.m` — extraction des sept ensembles de points.
- `analyser_delaunay_7_frontieres.m` — Delaunay global, statistiques et qualité des tétraèdres.
- `nettoyer_tetra_7_frontieres.m` — nettoyage des tétraèdres (qualité, géométrie, interface).
- `calculer_boules_7_frontieres.m` — sphères circonscrites, dédoublonnage, classification.
- `analyser_boules_7_frontieres_etape3.m` — analyse des banques de sphères.
- `filtrer_boules_TD_7_frontieres.m` — filtrage des sphères sur le débordement TD.
- `selection_TC_volumique_7_frontieres.m` — sélection par chevauchement TC.
- `couverture_union_avant_TC_7_frontieres.m` — unions et couverture avant la sélection TC.

### Benchmarks et variantes conservées

- `main_benchmark_7_frontieres_64.m` — pipeline 64³, grille sous-voxelique `k=4`.
- `main_benchmark_7_frontieres_64_geoTD.m` — variante 64³ avec traitement géométrique de TD.
- `main_benchmark_7_frontieres_512.m` — variante 512³ historique (très coûteuse en mémoire).
- `main_benchmark_7_frontieres_512_reprise.m` — reprise des calculs 512³.
- `main_benchmark_7_frontieres_step1_rapide.m` — pipeline sans subdivision `k=1`, TD par centres de voxels, variantes TD0/5/10 et TC0/5/10 ; possède une option de contrôle géométrique TD0.
- `main_benchmark_7_frontieres_step1_sequentiel_audit.m` — pipeline séquentiel avec audits TD.
- `main_benchmark_7_frontieres_step1_cubes_intersectes.m` — variante d'évaluation volumique.
- `main_benchmark_7_frontieres_k1_zero.m`, `main_benchmark_7_frontieres_k1_zero_TD0_5_10.m`, `main_benchmark_7_frontieres_k2_zero_TD0_5_10_TC0_5_10.m` — autres protocoles expérimentaux.

### Audits, géométrie et évaluation

- `verifier_contenance_geometrique_boules_retenues.m`
- `auditer_TD_volume_boules_retenues.m`
- `auditer_TD_continu_7_frontieres.m`
- `auditer_TD_continu_7_frontieres_v2.m`
- `budget_fixe_7_frontieres.m`
- `lancer_budget_fixe_BD12_DRP372.m`
- `calculer_MSD_RSAE_7_frontieres.m`
- `01_surface_BD12.m`, `02_surface_DRP372.m`, `03_surface_Tore.m`
- `validation_physique_boules_DRP372.m`

### Tore synthétique et topologie

- `generer_tore_128_benchmark.m`
- `analyser_tore_128_3_points.m`
- `analyser_topologie_tore_3_configurations.m`
- `analyser_surface_tore_TD0_TC0.m`
- `analyser_surface_tore_TD5_TC5.m`

## Prérequis

- MATLAB récent, avec **Image Processing Toolbox**.
- **Statistics and Machine Learning Toolbox** pour les scripts utilisant `knnsearch`.
- **Parallel Computing Toolbox** pour certaines variantes utilisant `parfor`. Consulter les options des scripts ; toutes les variantes ne sont pas séquentielles.
- Mémoire vive et espace disque importants pour 256³ et surtout 512³ : les triangulations globales et certaines grilles sous-voxelisées peuvent être très coûteuses.

Aucune version MATLAB précise n'a été validée lors de la préparation de cette archive.

## Démarrage rapide

Depuis MATLAB, ouvrir la racine du dépôt et ajouter ce dossier au path :

```matlab
addpath(pwd);

% Petit volume synthétique ; 1 = pore, 0 = solide.
img = false(16,16,16);
img(4:13,4:13,4:13) = true;

% Vérifier l'extraction, sans lancer le pipeline lourd :
[P, stats] = extraire_7_frontieres(img);
disp(stats);
```

Pour lancer l'ensemble de la chaîne sur un **véritable volume 64³** déjà chargé comme `img`, avec la convention `1/255 = pore` :

```matlab
assert(isequal(size(img), [64 64 64]));
resultats = main_benchmark_7_frontieres_64(img, 'resultats_delaunay_64');
```

Pour un volume 3D d'une taille quelconque avec le protocole `k=1` :

```matlab
% Le dossier de sortie doit être inexistant ou vide pour cet entrypoint.
resultats = main_benchmark_7_frontieres_step1_rapide(img, 'resultats_essai_neuf');
```

Pour générer le tore analytique 128³ sans lancer le benchmark :

```matlab
[imgTore, ~] = generer_tore_128_benchmark(false);
```

Pour le lancer avec la variante utilisée dans le générateur de tore :

```matlab
[imgTore, resultatsTore] = generer_tore_128_benchmark(true);
```

## Attention à la convention des images BD12 et DRP-372

Les scripts du benchmark attendent `0 = solide`, `1/255 = pore`, **mais certains fichiers RAW/MAT d'origine codent `0 = pore`**. Il faut inverser le masque au chargement :

```matlab
img = uint8(imageOriginale == 0) * 255;
```

Pour DRP-372, le script `02_surface_DRP372.m` lit la variable `bin` du fichier `374_02_024_256.mat` et utilise `bin == 0` pour le pore. Pour BD12, `01_surface_BD12.m` attend par défaut `bd12.mat` et une variable `img` ; les noms de fichiers et de variables doivent être adaptés à l'environnement local. Les scripts de surface supposent les dossiers de benchmark déjà calculés.

**Jeux de données attendus par les scripts historiques :** BD12 512³ (23 056 832 voxels poreux), DRP-372 256³ (688 309) et tore 128³ (182 384). Les images d'origine et sorties volumineuses sont exclues du dépôt.

## Reproductibilité, limitations et vigilance

- Les paramètres et les détails d'algorithme sont dans les en-têtes des fichiers et dans leurs options (notamment seuil de qualité `q=10^-3` et tolérance de dédoublonnage des sphères `10^-7` voxel pour les programmes concernés).
- **TD0 voxelique n'est pas TD0 géométrique.** La variante `step1_rapide` l'indique explicitement : aucun centre de voxel solide dans la sphère n'implique pas une inclusion géométrique exacte.
- Les résultats des méthodes changent selon la variante (`k=1`, `k=2`, `k=4`, audit continu, etc.). Préciser dans les publications ou issues quel `main_...` et quelles options ont produit les résultats.
- Les scripts `01_surface_...` et `lancer_budget_fixe_...` contiennent des chemins relatifs et noms de jeux de données à adapter. **Attention :** le chargeur BD12 du lanceur à budget fixe cherche une variable avec `nnz(x) == 23056832` ; si votre fichier original code `0 = pore`, convertissez ce fichier ou adaptez le chargeur, faute de quoi il ne reconnaîtra pas la bonne variable.
- Aucune exécution MATLAB ni reproduction des tableaux de publication n'a pu être effectuée dans l'environnement de préparation de l'archive ; valider les résultats avant citation ou publication définitive.
- Ne pas ajouter au dépôt les grands fichiers `.mat`, `.raw`, images et résultats générés sans vérifier leurs licences et leur provenance.

## Mise en ligne sur GitHub

Créer un dépôt GitHub (par exemple `3d-pore-boundary-benchmark`), décompresser cette archive localement, ouvrir un terminal **dans le dossier `Frontieres3D_GitHub`**, puis :

```bash
git init
git add .
git commit -m "Initial MATLAB boundary benchmark source"
git branch -M main
git remote add origin https://github.com/<utilisateur>/<depot>.git
git push -u origin main
```

Ne publier une `LICENSE` que lorsque les droits de distribution des contributions et des données ont été clarifiés. Ajouter la référence bibliographique définitive et les crédits des contributeurs avant la mise en ligne publique.

---

**English summary:** Original MATLAB research scripts for comparing seven 3D pore/solid boundary point sets (I6/I18/I26, E6/E18/E26, F), Delaunay circumsphere reconstruction, deviation/overlap filtering, geometry and topology analyses, and BD12/DRP-372/torus benchmarks. Research archive; no datasets or execution verification included.
