use assert_cmd::Command;
use assert_cmd::cargo;
use assert_fs::prelude::*;
use flate2::Compression;
use flate2::write::GzEncoder;
use predicates::prelude::*;
use std::fs;
use std::fs::File;
use std::io::Write;

// ================================================================
// UTILITAIRE : BUILDER DE COMMANDE
// ================================================================

fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fastdedup")
}

// ================================================================
// UTILITAIRE : GÉNÉRATION DE SÉQUENCES
// ================================================================

/// Génère une séquence ACGT déterministe de 150 bases via LCG.
///
/// Remplace l'encodage en base 5 précédent :
///   - N'utilise que ACGT (pas de N, conforme NovaSeq haute qualité)
///   - Distribution équilibrée (~25% par base) pour toutes les plages d'indices
///   - Garantit l'unicité en pratique pour tous les indices utilisés dans les tests
///     (max ~200 ; probabilité de collision dans 4^150 séquences ≈ 0)
///
/// Convention d'offset :
///   offset = 0    → séquences R1 / Single-End
///   offset = 9999 → séquences R2 (distinctes pour le même index)
fn generer_sequence(index: usize, offset: usize) -> String {
    const LONGUEUR: usize = 150;
    const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    // Graine dérivée de (index, offset) ; +1 pour éviter l'état nul
    let mut state = index
        .wrapping_mul(2_654_435_761)
        .wrapping_add(offset.wrapping_mul(40_503))
        .wrapping_add(1) as u64;
    let mut seq = Vec::with_capacity(LONGUEUR);
    for _ in 0..LONGUEUR {
        // LCG de Knuth (période maximale sur u64)
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        seq.push(NUCLEOTIDES[((state >> 33) & 3) as usize]);
    }
    String::from_utf8(seq).expect("Séquence UTF-8 invalide")
}

/// Génère des scores de qualité NovaSeq réalistes sur 150 cycles.
///
/// NovaSeq utilise exclusivement 4 bins de qualité :
///   Q2  → '#' (ASCII 35) — bases très basses qualité
///   Q12 → '-' (ASCII 45) — bases dégradées
///   Q23 → '8' (ASCII 56) — qualité intermédiaire
///   Q37 → 'F' (ASCII 70) — haute qualité (majorité)
///
/// Distribution approximative : 5% Q2, 5% Q12, 15% Q23, 75% Q37
/// (profil typique de lecture NovaSeq 2×150 sur génome entier).
fn generer_qualite(index: usize) -> String {
    const BINS: [u8; 4] = [b'#', b'-', b'8', b'F'];
    let mut state = index.wrapping_mul(2_654_435_761).wrapping_add(99_991) as u64;
    let mut qual = Vec::with_capacity(150);
    for _ in 0..150 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let r = (state >> 33) % 100;
        let bin = if r < 5 {
            0
        } else if r < 10 {
            1
        } else if r < 25 {
            2
        } else {
            3
        };
        qual.push(BINS[bin]);
    }
    String::from_utf8(qual).unwrap()
}

// ================================================================
// UTILITAIRE : GÉNÉRATEURS DE FICHIERS
// ================================================================

fn generer_fastq(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    for i in 0..nb_uniques {
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
    }
    for i in 0..nb_duplications {
        // Même séquence qu'à l'index i, qualité différente (duplication PCR réaliste)
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:5000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i,
            generer_sequence(i, 0),
            generer_qualite(i + 1000)
        ));
    }

    fichier.write_str(&contenu).unwrap();
    fichier
}

/// Crée un fichier FASTQ compressé (.fastq.gz) utilisable comme entrée.
fn generer_fastq_gz(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    for i in 0..nb_uniques {
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
    }
    for i in 0..nb_duplications {
        contenu.push_str(&format!(
            "@A00123:456:HFWV2DSXX:1:1101:5000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
            1000 + i,
            generer_sequence(i, 0),
            generer_qualite(i + 1000)
        ));
    }

    let file = File::create(fichier.path()).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(contenu.as_bytes()).unwrap();
    encoder.finish().unwrap();
    fichier
}

fn generer_fasta(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();

    for i in 0..nb_uniques {
        contenu.push_str(&format!(
            ">A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n",
            1000 + i,
            generer_sequence(i, 0)
        ));
    }
    for i in 0..nb_duplications {
        contenu.push_str(&format!(
            ">A00123:456:HFWV2DSXX:1:1101:5000:{} 1:N:0:ATGC\n{}\n",
            1000 + i,
            generer_sequence(i, 0)
        ));
    }

    fichier.write_str(&contenu).unwrap();
    fichier
}

fn generer_fastq_paire(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();

    for i in 0..nb_uniques {
        // L'ID de base est identique pour R1 et R2 ; seul le champ direction (1: vs 2:) diffère.
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:1000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            "@{} 1:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
        contenu_r2.push_str(&format!(
            "@{} 2:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 9999),
            generer_qualite(i + 500)
        ));
    }
    for i in 0..nb_duplications {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:5000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            "@{} 1:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 0),
            generer_qualite(i + 1000)
        ));
        contenu_r2.push_str(&format!(
            "@{} 2:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 9999),
            generer_qualite(i + 1500)
        ));
    }

    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

fn generer_fasta_paire(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_uniques: usize,
    nb_duplications: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();

    for i in 0..nb_uniques {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:1000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            ">{} 1:N:0:ATGC\n{}\n",
            base_cluster,
            generer_sequence(i, 0)
        ));
        contenu_r2.push_str(&format!(
            ">{} 2:N:0:ATGC\n{}\n",
            base_cluster,
            generer_sequence(i, 9999)
        ));
    }
    for i in 0..nb_duplications {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:5000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            ">{} 1:N:0:ATGC\n{}\n",
            base_cluster,
            generer_sequence(i, 0)
        ));
        contenu_r2.push_str(&format!(
            ">{} 2:N:0:ATGC\n{}\n",
            base_cluster,
            generer_sequence(i, 9999)
        ));
    }

    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

/// Génère une paire de fichiers FASTQ avec un nombre de reads délibérément
/// différent entre R1 et R2, pour tester la détection de désynchronisation
/// en fin de fichier.
fn generer_fastq_paire_desynchronisee(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_r1: usize,
    nb_r2: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();

    // Portion commune : IDs identiques, boucle synchronisée
    let nb_common = nb_r1.min(nb_r2);
    for i in 0..nb_common {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:1000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            "@{} 1:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
        contenu_r2.push_str(&format!(
            "@{} 2:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 9999),
            generer_qualite(i + 500)
        ));
    }
    // Reads supplémentaires côté R1 (si R1 > R2)
    for i in nb_common..nb_r1 {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:1000:{}", 1000 + i);
        contenu_r1.push_str(&format!(
            "@{} 1:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
    }
    // Reads supplémentaires côté R2 (si R2 > R1)
    for i in nb_common..nb_r2 {
        let base_cluster = format!("A00123:456:HFWV2DSXX:1:1101:1000:{}", 1000 + i);
        contenu_r2.push_str(&format!(
            "@{} 2:N:0:ATGC\n{}\n+\n{}\n",
            base_cluster,
            generer_sequence(i, 9999),
            generer_qualite(i + 500)
        ));
    }

    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

// ================================================================
// UTILITAIRE : HELPERS D'ANALYSE DE SORTIE
// ================================================================

fn decompresser_gz(chemin: &std::path::Path) -> String {
    let octets = fs::read(chemin).unwrap();
    let mut decoder = flate2::read::GzDecoder::new(&octets[..]);
    let mut contenu = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut contenu).unwrap();
    contenu
}

fn compter_sequences_fastq(contenu: &str) -> usize {
    contenu.lines().filter(|l| l.starts_with('@')).count()
}

fn compter_sequences_fasta(contenu: &str) -> usize {
    contenu.lines().filter(|l| l.starts_with('>')).count()
}

/// Extrait la liste des identifiants de base (avant l'espace) depuis un contenu FASTQ.
/// Permet une vérification de contenu, pas seulement de compte.
fn extraire_ids_fastq(contenu: &str) -> Vec<String> {
    contenu
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| {
            l.trim_start_matches('@')
                .split_whitespace()
                .next()
                .unwrap_or(l)
                .to_string()
        })
        .collect()
}

// ================================================================
// TESTS SINGLE-END — ERREURS D'ENTRÉE
// ================================================================

#[test]
fn test_fichier_entree_inexistant() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    cmd()
        .current_dir(temp_dir.path())
        .arg("-1")
        .arg("fichier_fantome.fastq")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Input file not found"));
}

/// Ne teste pas de chaîne exacte issue de needletail (fragile aux mises à jour
/// de la bibliothèque) — seul l'échec est vérifié.
#[test]
fn test_fichier_entree_vide() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_vide = temp_dir.child("vide.fastq");
    fichier_vide.touch().unwrap();
    cmd()
        .arg("-1")
        .arg(fichier_vide.path())
        .arg("-s")
        .assert()
        .failure();
}

#[test]
fn test_fichier_entree_mauvais_format() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_texte = temp_dir.child("mauvais_format.txt");
    fichier_texte
        .write_str("Ceci n'est pas un fichier FASTQ\nAvec plusieurs lignes.")
        .unwrap();
    cmd()
        .arg("-1")
        .arg(fichier_texte.path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("Invalid sequence data")
                .or(predicate::str::contains("Failed to read input file")),
        );
}

#[test]
fn test_fichier_sortie_chemin_inexistant() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    cmd()
        .arg("-1")
        .arg(fichier_entree.path())
        .arg("-o")
        .arg("/chemin/vers/un/dossier/inexistant/sortie.fastq")
        .assert()
        .failure();
}

#[test]
fn test_fichier_sortie_gz_corrompu() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier_entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let fichier_sortie_corrompu = temp_dir.child("sortie_corrompue.fastq.gz");
    fichier_sortie_corrompu
        .write_str("Ceci n'est pas un gzip valide.")
        .unwrap();
    cmd()
        .arg("-1")
        .arg(fichier_entree.path())
        .arg("-o")
        .arg(fichier_sortie_corrompu.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("Error opening preload file"));
}

// ================================================================
// TESTS SINGLE-END — TAUX DE DUPLICATION
// ================================================================

#[test]
fn test_taux_duplication_zero_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "100_uniques.fastq", 100, 0);
    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 100"))
        .stdout(predicate::str::contains("Duplicates removed: 0.00%"));
}

#[test]
fn test_taux_duplication_9_09_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 100 uniques + 10 duplicats = 110 lectures, 10/110 ≈ 9.09%
    let entree = generer_fastq(&temp_dir, "110_seq.fastq", 100, 10);
    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 110"))
        .stdout(predicate::str::contains("Duplicates removed: 9.09%"));
}

#[test]
fn test_taux_duplication_50_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 10 uniques + 10 duplicats = 20 lectures, 10/20 = 50.00%
    let entree = generer_fastq(&temp_dir, "50_pourcent.fastq", 10, 10);
    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 20"))
        .stdout(predicate::str::contains("Duplicates removed: 50.00%"));
}

/// Taux de duplication typique NovaSeq sur librairie amplicon/faible-input (~30%).
#[test]
fn test_taux_duplication_30_pourcent_novaseq() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 70 uniques + 30 duplicats = 100 lectures, 30/100 = 30.00%
    let entree = generer_fastq(&temp_dir, "30_pourcent.fastq", 70, 30);
    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 100"))
        .stdout(predicate::str::contains("Duplicates removed: 30.00%"));
}

#[test]
fn test_taux_duplication_entree_toutes_identiques() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 10 lectures identiques : 1 unique + 9 doublons = 90.00%
    let fichier = temp_dir.child("tous_identiques.fastq");
    let seq = generer_sequence(42, 0);
    let mut contenu = String::new();
    for i in 0..10 {
        contenu.push_str(&format!(
            "@read_{} 1:N:0:ATGC\n{}\n+\n{}\n",
            i,
            seq,
            generer_qualite(i)
        ));
    }
    fichier.write_str(&contenu).unwrap();

    cmd()
        .arg("-1")
        .arg(fichier.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 10"))
        .stdout(predicate::str::contains("Duplicates removed: 90.00%"));
}

// ================================================================
// TESTS SINGLE-END — ENTRÉE COMPRESSÉE (GZ)
// ================================================================

/// Vérifie que needletail lit correctement un fichier .fastq.gz en entrée.
#[test]
fn test_entree_fastq_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq_gz(&temp_dir, "entree.fastq.gz", 20, 5);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(
        compter_sequences_fastq(&contenu),
        20,
        "Doit contenir 20 séquences uniques depuis une entrée .fastq.gz"
    );
}

/// Vérifie le taux de duplication en dry-run sur une entrée .fastq.gz.
/// 50 uniques + 10 duplicats = 60 lectures, 10/60 ≈ 16.67%
#[test]
fn test_entree_fastq_gz_dry_run() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq_gz(&temp_dir, "entree.fastq.gz", 50, 10);

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Processed fragments: 60"))
        .stdout(predicate::str::contains("Duplicates removed: 16.67%"));
}

// ================================================================
// TESTS SINGLE-END — FORMAT ET COMPRESSION DE SORTIE
// ================================================================

#[test]
fn test_sortie_fastq_non_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(
        contenu.starts_with('@'),
        "La sortie FASTQ doit commencer par '@'"
    );
    assert!(
        contenu.contains('+'),
        "La sortie FASTQ doit contenir la ligne '+'"
    );
    assert_eq!(compter_sequences_fastq(&contenu), 10);
}

#[test]
fn test_sortie_fastq_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 3);
    let sortie = temp_dir.child("sortie.fastq.gz");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f, "Magic byte gzip attendu");
    assert_eq!(octets[1], 0x8b, "Magic byte gzip attendu");

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu), 10);
}

#[test]
fn test_sortie_fq_extension() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 5, 0);
    let sortie = temp_dir.child("sortie.fq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(
        contenu.starts_with('@'),
        "L'extension .fq doit produire du FASTQ"
    );
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_sortie_fq_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 5, 0);
    let sortie = temp_dir.child("sortie.fq.gz");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f);
    assert_eq!(octets[1], 0x8b);

    let contenu = decompresser_gz(sortie.path());
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_sortie_fasta_extension_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 2);
    let sortie = temp_dir.child("sortie.fasta");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(
        contenu.starts_with('>'),
        "La sortie .fasta doit commencer par '>'"
    );
    assert!(
        !contenu.contains('+'),
        "La sortie FASTA ne doit pas contenir de ligne '+'"
    );
    assert_eq!(compter_sequences_fasta(&contenu), 10);
}

#[test]
fn test_sortie_fasta_extension_fa() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 7, 0);
    let sortie = temp_dir.child("sortie.fa");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(
        contenu.starts_with('>'),
        "L'extension .fa doit produire du FASTA"
    );
    assert_eq!(compter_sequences_fasta(&contenu), 7);
}

#[test]
fn test_sortie_fasta_extension_fna() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 7, 0);
    let sortie = temp_dir.child("sortie.fna");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(
        contenu.starts_with('>'),
        "L'extension .fna doit produire du FASTA"
    );
    assert_eq!(compter_sequences_fasta(&contenu), 7);
}

#[test]
fn test_sortie_fasta_gz_compresse() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 8, 2);
    let sortie = temp_dir.child("sortie.fasta.gz");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f);
    assert_eq!(octets[1], 0x8b);

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 8);
}

#[test]
fn test_entree_fasta_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fasta(&temp_dir, "entree.fasta", 10, 3);
    let sortie = temp_dir.child("sortie.fasta");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = fs::read_to_string(sortie.path()).unwrap();
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 10);
}

/// La conversion FASTA → FASTQ est impossible (pas de scores de qualité en FASTA).
/// Le binaire doit refuser explicitement cette opération.
#[test]
fn test_entree_fasta_sortie_fastq_non_supporte() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fasta(&temp_dir, "entree.fasta", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "FASTA → FASTQ conversion not supported",
        ));
}

// ================================================================
// TESTS SINGLE-END — COMPORTEMENT --force ET REPRISE
// ================================================================

#[test]
fn test_force_overwrites_existing_file() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree_1 = generer_fastq(&temp_dir, "entree_1.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree_1.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie.path()).unwrap()),
        10
    );

    // Second passage avec --force et un fichier plus petit : doit écraser, pas ajouter
    let entree_2 = generer_fastq(&temp_dir, "entree_2.fastq", 4, 0);
    cmd()
        .arg("-1")
        .arg(entree_2.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("--force")
        .assert()
        .success();

    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie.path()).unwrap()),
        4,
        "--force doit écraser le fichier existant, pas y ajouter"
    );
}

#[test]
fn test_without_force_resumes_from_existing_output() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree_initiale = generer_fastq(&temp_dir, "entree_initiale.fastq", 5, 0);
    let sortie = temp_dir.child("sortie_reprise.fastq");

    // Premier passage : 5 séquences uniques
    cmd()
        .arg("-1")
        .arg(entree_initiale.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    // Deuxième passage : les 5 mêmes + 5 nouvelles
    let entree_complete = generer_fastq(&temp_dir, "entree_complete.fastq", 10, 0);
    cmd()
        .arg("-1")
        .arg(entree_complete.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu_final = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(
        compter_sequences_fastq(&contenu_final),
        10,
        "Après reprise, la sortie doit contenir 10 séquences uniques au total"
    );

    // Vérification de contenu : aucun ID ne doit apparaître deux fois
    let ids = extraire_ids_fastq(&contenu_final);
    let mut ids_tries = ids.clone();
    ids_tries.sort();
    ids_tries.dedup();
    assert_eq!(
        ids.len(),
        ids_tries.len(),
        "Aucun ID de lecture ne doit être dupliqué dans la sortie de reprise"
    );
}

// ================================================================
// TESTS SINGLE-END — SÉLECTION AUTOMATIQUE DU MODE DE HACHAGE
// ================================================================

/// Un threshold très petit (→ 0) force la sélection du 128-bit pour toute taille de fichier.
/// Vérifie le chemin de sélection automatique sans nécessiter un fichier de 67 Go.
#[test]
fn test_selection_automatique_hash_128_bits_via_threshold() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("-t")
        .arg("0.000000000000000000001") // threshold ≈ 0 → seuil = ~0 → 128-bit
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("128-bit Hash"));
}

/// Un threshold de 1.0 (100% de collisions acceptées) force le 64-bit pour toute taille pratique.
#[test]
fn test_selection_automatique_hash_64_bits_via_threshold() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("-t")
        .arg("1.0") // seuil ≈ 6 milliards → fichiers de test bien en dessous
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("64-bit Hash"));
}

/// Vérifie que l'avertissement est émis quand --hash et --threshold sont fournis ensemble.
#[test]
fn test_avertissement_hash_et_threshold_ignores() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 0);
    let sortie = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("-H")
        .arg("64")
        .arg("-t")
        .arg("0.5")
        .assert()
        .success()
        .stderr(predicate::str::contains("threshold").and(predicate::str::contains("ignored")));
}

// ================================================================
// TESTS SINGLE-END — MODE DE HACHAGE MANUEL
// ================================================================

#[test]
fn test_hachage_64_bits_explicite() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 20, 5);
    let sortie = temp_dir.child("sortie_64.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("-H")
        .arg("64")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("64-bit Hash"));

    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie.path()).unwrap()),
        20
    );
}

#[test]
fn test_hachage_128_bits_explicite() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 20, 5);
    let sortie = temp_dir.child("sortie_128.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .arg("-H")
        .arg("128")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("128-bit Hash"));

    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie.path()).unwrap()),
        20
    );
}

#[test]
fn test_hachage_64_et_128_bits_produisent_le_meme_resultat() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 15, 5);
    let sortie_64 = temp_dir.child("sortie_64.fastq");
    let sortie_128 = temp_dir.child("sortie_128.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie_64.path())
        .arg("-H")
        .arg("64")
        .assert()
        .success();

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie_128.path())
        .arg("-H")
        .arg("128")
        .assert()
        .success();

    let count_64 = compter_sequences_fastq(&fs::read_to_string(sortie_64.path()).unwrap());
    let count_128 = compter_sequences_fastq(&fs::read_to_string(sortie_128.path()).unwrap());
    assert_eq!(
        count_64, count_128,
        "Les deux modes de hachage doivent produire le même nombre de séquences uniques"
    );
}

// ================================================================
// TESTS SINGLE-END — INTÉGRITÉ APRÈS COMPRESSION
// ================================================================

#[test]
fn test_compression_gz_valide_fastq() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 50, 10);
    let sortie = temp_dir.child("sortie.fastq.gz");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu), 50);
}

#[test]
fn test_compression_gz_valide_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 50, 10);
    let sortie = temp_dir.child("sortie.fasta.gz");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie.path())
        .assert()
        .success();

    let contenu = decompresser_gz(sortie.path());
    assert!(contenu.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu), 50);
}

#[test]
fn test_deduplication_preservee_apres_compression() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 10, 10);
    let sortie_gz = temp_dir.child("sortie.fastq.gz");
    let sortie_non_gz = temp_dir.child("sortie.fastq");

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie_gz.path())
        .assert()
        .success();

    cmd()
        .arg("-1")
        .arg(entree.path())
        .arg("-o")
        .arg(sortie_non_gz.path())
        .arg("--force")
        .assert()
        .success();

    let count_gz = compter_sequences_fastq(&decompresser_gz(sortie_gz.path()));
    let count_non_gz = compter_sequences_fastq(&fs::read_to_string(sortie_non_gz.path()).unwrap());

    assert_eq!(
        count_gz, count_non_gz,
        "Le nombre de séquences doit être identique entre .gz et non-.gz"
    );
    assert_eq!(
        count_gz, 10,
        "Doit contenir exactement 10 séquences uniques"
    );
}

// ================================================================
// TESTS PAIRED-END — ERREURS ET VALIDATION
// ================================================================

#[test]
fn test_paire_erreur_sortie_r2_manquante() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 0);
    let sortie_r1 = temp_dir.child("sortie_r1.fastq");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        // -p intentionnellement absent
        .assert()
        .failure()
        .stderr(predicate::str::contains("--output-r2").and(predicate::str::contains("required")));
}

#[test]
fn test_paire_desynchronisee_en_tete_bloquee() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1_desync.fastq", "R2_desync.fastq", 10, 0);

    // Insérer une lecture en tête de R2 : l'ID du premier record R2 ne correspond plus à R1
    let mut contenu_r2_corrompu =
        String::from("@A00123:456:HFWV2DSXX:1:1101:9999:9999 2:N:0:ATGC\nATGC\n+\nFFFF\n");
    contenu_r2_corrompu.push_str(&fs::read_to_string(r2.path()).unwrap());
    r2.write_str(&contenu_r2_corrompu).unwrap();

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-p")
        .arg(temp_dir.child("dummy_r2.fastq").path())
        .arg("-s")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Critical desynchronization detected",
        ));
}

/// Régression : les en-têtes pré-CASAVA (`/1` `/2`) doivent être acceptés par le
/// chemin d'écriture sérial. Le chemin parallèle les acceptait déjà ; les deux
/// doivent produire exactement la même sortie.
#[test]
fn test_paire_entetes_pre_casava_serial_et_parallele() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    // 2 paires uniques + 1 paire dupliquée = 3 paires, 1 duplicat PCR.
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();
    for (i, dup) in [(0usize, false), (1, false), (0, true)] {
        contenu_r1.push_str(&format!(
            "@read{}/1\n{}\n+\n{}\n",
            i + 1,
            generer_sequence(i, 0),
            generer_qualite(i + if dup { 1000 } else { 0 })
        ));
        contenu_r2.push_str(&format!(
            "@read{}/2\n{}\n+\n{}\n",
            i + 1,
            generer_sequence(i, 9999),
            generer_qualite(i + 500)
        ));
    }
    let r1_in = temp_dir.child("precasava_R1.fastq");
    let r2_in = temp_dir.child("precasava_R2.fastq");
    r1_in.write_str(&contenu_r1).unwrap();
    r2_in.write_str(&contenu_r2).unwrap();

    let serial_r1_path = temp_dir.child("serial_R1.fastq");
    let serial_r2_path = temp_dir.child("serial_R2.fastq");
    let par_r1_path = temp_dir.child("par_R1.fastq");
    let par_r2_path = temp_dir.child("par_R2.fastq");

    // Chemin sérial (-j 1) : c'est lui qui rejetait ces en-têtes comme désync.
    cmd()
        .arg("-1")
        .arg(r1_in.path())
        .arg("-2")
        .arg(r2_in.path())
        .arg("-o")
        .arg(serial_r1_path.path())
        .arg("-p")
        .arg(serial_r2_path.path())
        .arg("-j")
        .arg("1")
        .assert()
        .success();

    // Chemin parallèle (-j 4) : référence (il acceptait déjà ces en-têtes).
    cmd()
        .arg("-1")
        .arg(r1_in.path())
        .arg("-2")
        .arg(r2_in.path())
        .arg("-o")
        .arg(par_r1_path.path())
        .arg("-p")
        .arg(par_r2_path.path())
        .arg("-j")
        .arg("4")
        .assert()
        .success();

    let serial_r1 = fs::read_to_string(serial_r1_path.path()).unwrap();
    let serial_r2 = fs::read_to_string(serial_r2_path.path()).unwrap();
    let par_r1 = fs::read_to_string(par_r1_path.path()).unwrap();
    let par_r2 = fs::read_to_string(par_r2_path.path()).unwrap();

    // Le duplicat PCR est éliminé : 2 paires uniques écrites de chaque côté.
    assert_eq!(compter_sequences_fastq(&serial_r1), 2);
    assert_eq!(compter_sequences_fastq(&serial_r2), 2);

    // Le suffixe /1 /2 de l'entrée est préservé tel quel en sortie.
    assert!(serial_r1.contains("@read1/1"));
    assert!(serial_r2.contains("@read1/2"));

    // Sérial et parallèle doivent produire exactement le même résultat.
    assert_eq!(serial_r1, par_r1, "R1 sérial != R1 parallèle");
    assert_eq!(serial_r2, par_r2, "R2 sérial != R2 parallèle");
}

/// Le chemin d'écriture sérial doit toujours refuser une VRAIE désynchronisation
/// (IDs de base réellement différents) avec des en-têtes pré-CASAVA `/1` `/2` :
/// les deux premiers appariements passent, le troisième diverge.
#[test]
fn test_paire_entetes_pre_casava_desync_reelle_bloquee() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();
    for i in 0..5 {
        let n = i + 1;
        contenu_r1.push_str(&format!(
            "@read{}/1\n{}\n+\n{}\n",
            n,
            generer_sequence(i, 0),
            generer_qualite(i)
        ));
        let id_r2 = if n == 3 {
            format!("other{n}")
        } else {
            format!("read{n}")
        };
        contenu_r2.push_str(&format!(
            "@{}/2\n{}\n+\n{}\n",
            id_r2,
            generer_sequence(i, 9999),
            generer_qualite(i + 500)
        ));
    }
    let r1_in = temp_dir.child("precasava_desync_R1.fastq");
    let r2_in = temp_dir.child("precasava_desync_R2.fastq");
    r1_in.write_str(&contenu_r1).unwrap();
    r2_in.write_str(&contenu_r2).unwrap();

    cmd()
        .arg("-1")
        .arg(r1_in.path())
        .arg("-2")
        .arg(r2_in.path())
        .arg("-o")
        .arg(temp_dir.child("desync_out_R1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("desync_out_R2.fastq").path())
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("Critical desynchronization detected at pair #3")
                .and(predicate::str::contains("read3"))
                .and(predicate::str::contains("other3")),
        );
}

/// Vérifie la détection de désynchronisation de fin de fichier quand R1 > R2.
#[test]
fn test_paire_desync_fin_fichier_r1_plus_long() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // R1 : 12 reads, R2 : 10 reads — les 10 premiers s'apparient correctement
    let (r1, r2) = generer_fastq_paire_desynchronisee(&temp_dir, "R1.fastq", "R2.fastq", 12, 10);

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(temp_dir.child("out_r1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("out_r2.fastq").path())
        .arg("--force")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Desynchronization detected at the end",
        ));
}

/// Vérifie la détection de désynchronisation de fin de fichier quand R2 > R1.
#[test]
fn test_paire_desync_fin_fichier_r2_plus_long() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // R1 : 10 reads, R2 : 12 reads
    let (r1, r2) = generer_fastq_paire_desynchronisee(&temp_dir, "R1.fastq", "R2.fastq", 10, 12);

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(temp_dir.child("out_r1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("out_r2.fastq").path())
        .arg("--force")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Desynchronization detected at the end",
        ));
}

#[test]
fn test_paire_formats_incoherents() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 0);

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(temp_dir.child("sortie_R1.fasta").path())
        .arg("-p")
        .arg(temp_dir.child("sortie_R2.fastq").path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("must have the same format"));
}

/// La conversion FASTA → FASTQ doit être bloquée en mode Paired-End.
#[test]
fn test_paire_entree_fasta_sortie_fastq_non_supporte() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fasta_paire(&temp_dir, "R1.fasta", "R2.fasta", 10, 0);

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(temp_dir.child("sortie_R1.fastq").path())
        .arg("-p")
        .arg(temp_dir.child("sortie_R2.fastq").path())
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "FASTA → FASTQ conversion not supported",
        ));
}

// ================================================================
// TESTS PAIRED-END — TAUX DE DUPLICATION
// ================================================================

#[test]
fn test_paire_taux_duplication_33_pourcent() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // 100 paires uniques + 50 dupliquées = 150 paires, 50/150 ≈ 33.33%
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 100, 50);

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-p")
        .arg(temp_dir.child("dummy_r2.fastq").path())
        .arg("-s")
        .arg("-v")
        .assert()
        .success()
        .stdout(predicate::str::contains("Paired-End"))
        .stdout(predicate::str::contains("Processed fragments: 150"))
        .stdout(predicate::str::contains("Duplicates removed: 33.33%"));
}

// ================================================================
// TESTS PAIRED-END — FORMAT ET COMPRESSION DE SORTIE
// ================================================================

#[test]
fn test_paire_sortie_fastq() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('@'));
    assert!(contenu_r2.starts_with('@'));
    assert_eq!(
        compter_sequences_fastq(&contenu_r1),
        10,
        "R1 doit contenir 10 paires uniques"
    );
    assert_eq!(
        compter_sequences_fastq(&contenu_r2),
        10,
        "R2 doit contenir 10 paires uniques"
    );
    assert_eq!(
        compter_sequences_fastq(&contenu_r1),
        compter_sequences_fastq(&contenu_r2),
        "R1 et R2 doivent être synchronisés"
    );
}

#[test]
fn test_paire_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('>'));
    assert!(contenu_r2.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(compter_sequences_fasta(&contenu_r2), 10);
}

#[test]
fn test_paire_sortie_fastq_gz() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq.gz");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq.gz");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    for sortie in [&sortie_r1, &sortie_r2] {
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f);
        assert_eq!(octets[1], 0x8b);
    }

    let contenu_r1 = decompresser_gz(sortie_r1.path());
    let contenu_r2 = decompresser_gz(sortie_r2.path());

    assert!(contenu_r1.starts_with('@'));
    assert_eq!(compter_sequences_fastq(&contenu_r1), 10);
    assert_eq!(
        compter_sequences_fastq(&contenu_r1),
        compter_sequences_fastq(&contenu_r2)
    );
}

#[test]
fn test_paire_sortie_fasta_gz() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta.gz");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta.gz");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    for sortie in [&sortie_r1, &sortie_r2] {
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f);
        assert_eq!(octets[1], 0x8b);
    }

    let contenu_r1 = decompresser_gz(sortie_r1.path());
    let contenu_r2 = decompresser_gz(sortie_r2.path());

    assert!(contenu_r1.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(
        compter_sequences_fasta(&contenu_r1),
        compter_sequences_fasta(&contenu_r2)
    );
}

#[test]
fn test_paire_entree_fasta_sortie_fasta() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fasta_paire(&temp_dir, "R1.fasta", "R2.fasta", 10, 3);
    let sortie_r1 = temp_dir.child("sortie_R1.fasta");
    let sortie_r2 = temp_dir.child("sortie_R2.fasta");

    cmd()
        .arg("-1")
        .arg(r1.path())
        .arg("-2")
        .arg(r2.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert!(contenu_r1.starts_with('>'));
    assert!(contenu_r2.starts_with('>'));
    assert_eq!(compter_sequences_fasta(&contenu_r1), 10);
    assert_eq!(
        compter_sequences_fasta(&contenu_r1),
        compter_sequences_fasta(&contenu_r2)
    );
}

// ================================================================
// TESTS PAIRED-END — REPRISE (RESUME)
// ================================================================

/// Vérifie que la reprise sans --force fonctionne en mode paired-end :
/// après le second passage, R1 et R2 doivent contenir exactement le total
/// de paires uniques sans duplication de contenu.
#[test]
fn test_paire_reprise_depuis_sortie_existante() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    // Premier passage : 5 paires uniques
    let (r1_init, r2_init) = generer_fastq_paire(&temp_dir, "R1_init.fastq", "R2_init.fastq", 5, 0);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq");

    cmd()
        .arg("-1")
        .arg(r1_init.path())
        .arg("-2")
        .arg(r2_init.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie_r1.path()).unwrap()),
        5
    );
    assert_eq!(
        compter_sequences_fastq(&fs::read_to_string(sortie_r2.path()).unwrap()),
        5
    );

    // Deuxième passage : 5 mêmes paires + 5 nouvelles, sans --force
    let (r1_full, r2_full) =
        generer_fastq_paire(&temp_dir, "R1_full.fastq", "R2_full.fastq", 10, 0);

    cmd()
        .arg("-1")
        .arg(r1_full.path())
        .arg("-2")
        .arg(r2_full.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let final_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let final_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    assert_eq!(
        compter_sequences_fastq(&final_r1),
        10,
        "R1 doit contenir 10 paires uniques après la reprise"
    );
    assert_eq!(
        compter_sequences_fastq(&final_r2),
        10,
        "R2 doit contenir 10 paires uniques après la reprise"
    );
    assert_eq!(
        compter_sequences_fastq(&final_r1),
        compter_sequences_fastq(&final_r2),
        "R1 et R2 doivent rester synchronisés après la reprise"
    );
}

/// Vérifie que la synchronisation R1/R2 est préservée après une reprise :
/// le read à la position N dans R1 et le read à la position N dans R2
/// doivent avoir le même identifiant de cluster.
#[test]
fn test_paire_reprise_synchronisation_r1_r2() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    let (r1_init, r2_init) = generer_fastq_paire(&temp_dir, "R1_init.fastq", "R2_init.fastq", 8, 0);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq");

    cmd()
        .arg("-1")
        .arg(r1_init.path())
        .arg("-2")
        .arg(r2_init.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let (r1_full, r2_full) =
        generer_fastq_paire(&temp_dir, "R1_full.fastq", "R2_full.fastq", 15, 0);

    cmd()
        .arg("-1")
        .arg(r1_full.path())
        .arg("-2")
        .arg(r2_full.path())
        .arg("-o")
        .arg(sortie_r1.path())
        .arg("-p")
        .arg(sortie_r2.path())
        .assert()
        .success();

    let contenu_r1 = fs::read_to_string(sortie_r1.path()).unwrap();
    let contenu_r2 = fs::read_to_string(sortie_r2.path()).unwrap();

    // L'ID de base est la partie avant l'espace (commun à R1 et R2 pour chaque paire)
    let base_ids_r1: Vec<&str> = contenu_r1
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| {
            l.trim_start_matches('@')
                .split_whitespace()
                .next()
                .unwrap_or(l)
        })
        .collect();

    let base_ids_r2: Vec<&str> = contenu_r2
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| {
            l.trim_start_matches('@')
                .split_whitespace()
                .next()
                .unwrap_or(l)
        })
        .collect();

    assert_eq!(
        base_ids_r1.len(),
        base_ids_r2.len(),
        "R1 et R2 doivent avoir le même nombre de reads"
    );
    for (pos, (id_r1, id_r2)) in base_ids_r1.iter().zip(base_ids_r2.iter()).enumerate() {
        assert_eq!(
            id_r1, id_r2,
            "Désynchronisation à la position {} : R1={}, R2={}",
            pos, id_r1, id_r2
        );
    }
}

// ================================================================
// TESTS — VALIDATION DES SORTIES .GZ
// ================================================================

/// Vérifie qu'une sortie .gz en mode SE a les bons magic bytes et contient
/// le bon nombre de séquences après décompression.
#[test]
fn test_validation_sortie_gz_single_end() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 20, 5);
    let sortie = temp_dir.child("sortie.fastq.gz");

    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie.path())
        .arg("--threads").arg("4")
        .assert().success();

    // Vérifier les magic bytes gzip (0x1f 0x8b)
    let octets = fs::read(sortie.path()).unwrap();
    assert_eq!(octets[0], 0x1f, "Premier octet doit être le magic byte gzip 0x1f");
    assert_eq!(octets[1], 0x8b, "Second octet doit être le magic byte gzip 0x8b");

    // Vérifier que flate2 peut décompresser
    let contenu = decompresser_gz(sortie.path());
    let count = compter_sequences_fastq(&contenu);
    assert_eq!(count, 20, "Sortie .gz doit contenir 20 séquences uniques");
}

/// Vérifie qu'une sortie .gz en mode multi-thread est valide et
/// que le compte de séquences est correct.
#[test]
fn test_validation_sortie_gz_threads_multiple() {
    let temp_dir = assert_fs::TempDir::new().unwrap();

    for threads in &[2usize, 4, 8] {
        let entree = generer_fastq(&temp_dir, &format!("entree_t{threads}.fastq"), 50, 10);
        let sortie = temp_dir.child(&format!("sortie_t{threads}.fastq.gz"));

        cmd()
            .arg("-1").arg(entree.path())
            .arg("-o").arg(sortie.path())
            .arg("--threads").arg(threads.to_string())
            .assert().success();

        // Magic bytes
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f, "threads={threads}: magic byte gzip");
        assert_eq!(octets[1], 0x8b, "threads={threads}: magic byte gzip");

        // Décompresser avec flate2 et vérifier le contenu
        let contenu = decompresser_gz(sortie.path());
        let count = compter_sequences_fastq(&contenu);
        assert_eq!(
            count, 50,
            "threads={threads}: sortie .gz doit contenir 50 séquences uniques"
        );
    }
}

/// Vérifie qu'une sortie .gz en mode PE est valide, que R1 et R2 sont
/// synchronisés et que le contenu est correct.
#[test]
fn test_validation_sortie_gz_paired_end() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1, r2) = generer_fastq_paire(&temp_dir, "R1.fastq", "R2.fastq", 30, 8);
    let sortie_r1 = temp_dir.child("sortie_R1.fastq.gz");
    let sortie_r2 = temp_dir.child("sortie_R2.fastq.gz");

    cmd()
        .arg("-1").arg(r1.path())
        .arg("-2").arg(r2.path())
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .arg("--threads").arg("4")
        .assert().success();

    // Vérifier les magic bytes
    for (nom, sortie) in [("R1", &sortie_r1), ("R2", &sortie_r2)] {
        let octets = fs::read(sortie.path()).unwrap();
        assert_eq!(octets[0], 0x1f, "{nom}: premier magic byte gzip");
        assert_eq!(octets[1], 0x8b, "{nom}: second magic byte gzip");
    }

    // Décompresser et vérifier la synchronisation
    let contenu_r1 = decompresser_gz(sortie_r1.path());
    let contenu_r2 = decompresser_gz(sortie_r2.path());

    let count_r1 = compter_sequences_fastq(&contenu_r1);
    let count_r2 = compter_sequences_fastq(&contenu_r2);

    assert_eq!(count_r1, 30, "R1 doit contenir 30 séquences uniques");
    assert_eq!(count_r2, 30, "R2 doit contenir 30 séquences uniques");

    // Vérifier la synchronisation R1/R2
    let ids_r1: Vec<&str> = contenu_r1
        .lines().filter(|l| l.starts_with('@'))
        .map(|l| l.trim_start_matches('@').split_whitespace().next().unwrap_or(l))
        .collect();
    let ids_r2: Vec<&str> = contenu_r2
        .lines().filter(|l| l.starts_with('@'))
        .map(|l| l.trim_start_matches('@').split_whitespace().next().unwrap_or(l))
        .collect();

    for (pos, (id1, id2)) in ids_r1.iter().zip(ids_r2.iter()).enumerate() {
        assert_eq!(id1, id2, "Désynchronisation .gz à la position {pos}: R1={id1}, R2={id2}");
    }
}

/// Vérifie que le décompte de séquences est identique entre une sortie .gz
/// et une sortie non compressée, en mode multi-thread.
#[test]
fn test_validation_sortie_gz_vs_non_gz_multi_thread() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq(&temp_dir, "entree.fastq", 100, 25);
    let sortie_gz = temp_dir.child("sortie.fastq.gz");
    let sortie_plain = temp_dir.child("sortie.fastq");

    // Sortie compressée (multi-thread)
    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_gz.path())
        .arg("--threads").arg("4")
        .assert().success();

    // Sortie non compressée
    cmd()
        .arg("-1").arg(entree.path())
        .arg("-o").arg(sortie_plain.path())
        .arg("--force")
        .assert().success();

    let count_gz = compter_sequences_fastq(&decompresser_gz(sortie_gz.path()));
    let count_plain = compter_sequences_fastq(&fs::read_to_string(sortie_plain.path()).unwrap());

    assert_eq!(
        count_gz, count_plain,
        "Sortie .gz et non compressée doivent avoir le même nombre de séquences"
    );
    assert_eq!(count_gz, 100, "Doit contenir exactement 100 séquences uniques");
}

// ================================================================
// TESTS PAIRED-END — REPRISE / RÉSILIANCE
// ================================================================
