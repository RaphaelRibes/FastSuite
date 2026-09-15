use assert_cmd::Command;
use assert_cmd::cargo;
use assert_fs::prelude::*;
use predicates::prelude::*;
use std::fs;
use std::io::Write;

// ================================================================
// HELPERS
// ================================================================

fn cmd() -> Command {
    cargo::cargo_bin_cmd!("fastdedup")
}

fn generer_sequence(index: usize, offset: usize) -> String {
    const LONGUEUR: usize = 150;
    const NUCLEOTIDES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut state = index
        .wrapping_mul(2_654_435_761)
        .wrapping_add(offset.wrapping_mul(40_503))
        .wrapping_add(1) as u64;
    let mut seq = Vec::with_capacity(LONGUEUR);
    for _ in 0..LONGUEUR {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        seq.push(NUCLEOTIDES[((state >> 33) & 3) as usize]);
    }
    String::from_utf8(seq).expect("Séquence UTF-8 invalide")
}

fn generer_qualite(index: usize) -> String {
    const BINS: [u8; 4] = [b'#', b'-', b'8', b'F'];
    let mut state = index.wrapping_mul(2_654_435_761).wrapping_add(99_991) as u64;
    let mut qual = Vec::with_capacity(150);
    for _ in 0..150 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let r = (state >> 33) % 100;
        let bin = if r < 5 { 0 } else if r < 10 { 1 } else if r < 25 { 2 } else { 3 };
        qual.push(BINS[bin]);
    }
    String::from_utf8(qual).unwrap()
}

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

/// Generates a FASTQ where N unique sequences each appear `repetitions` times.
fn generer_fastq_repetitions(
    dir: &assert_fs::TempDir,
    nom_fichier: &str,
    nb_uniques: usize,
    repetitions: usize,
) -> assert_fs::fixture::ChildPath {
    let fichier = dir.child(nom_fichier);
    let mut contenu = String::new();
    for _rep in 0..repetitions {
        for i in 0..nb_uniques {
            contenu.push_str(&format!(
                "@A00123:456:HFWV2DSXX:1:1101:1000:{} 1:N:0:ATGC\n{}\n+\n{}\n",
                1000 + i,
                generer_sequence(i, 0),
                generer_qualite(i)
            ));
        }
    }
    fichier.write_str(&contenu).unwrap();
    fichier
}

/// Generates PE FASTQ where N unique pairs each appear `repetitions` times.
fn generer_fastq_paire_repetitions(
    dir: &assert_fs::TempDir,
    nom_r1: &str,
    nom_r2: &str,
    nb_uniques: usize,
    repetitions: usize,
) -> (assert_fs::fixture::ChildPath, assert_fs::fixture::ChildPath) {
    let fichier_r1 = dir.child(nom_r1);
    let fichier_r2 = dir.child(nom_r2);
    let mut contenu_r1 = String::new();
    let mut contenu_r2 = String::new();
    for _rep in 0..repetitions {
        for i in 0..nb_uniques {
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
    }
    fichier_r1.write_str(&contenu_r1).unwrap();
    fichier_r2.write_str(&contenu_r2).unwrap();
    (fichier_r1, fichier_r2)
}

// ================================================================
// TESTS
// ================================================================

#[test]
fn test_abundance_single_end_annotates_size() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq_repetitions(&temp_dir, "entree.fastq", 5, 3);
    let sortie = temp_dir.child("out.fq");
    cmd()
        .arg("-1").arg(entree.path())
        .arg("--abundance")
        .arg("-o").arg(sortie.path())
        .assert().success();
    let contenu = fs::read_to_string(sortie.path()).unwrap();
    for ligne in contenu.lines() {
        if ligne.starts_with('@') {
            assert!(ligne.ends_with(";size=3"),
                "L'en-tête doit se terminer par ;size=3, obtenu : {}", ligne);
        }
    }
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_abundance_paired_end_annotates_size() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let (r1_in, r2_in) = generer_fastq_paire_repetitions(&temp_dir, "R1.fastq", "R2.fastq", 5, 3);
    let sortie_r1 = temp_dir.child("out_R1.fq");
    let sortie_r2 = temp_dir.child("out_R2.fq");
    cmd()
        .arg("-1").arg(r1_in.path())
        .arg("-2").arg(r2_in.path())
        .arg("--abundance")
        .arg("-o").arg(sortie_r1.path())
        .arg("-p").arg(sortie_r2.path())
        .assert().success();
    for (nom, sortie) in [("R1", &sortie_r1), ("R2", &sortie_r2)] {
        let contenu = fs::read_to_string(sortie.path()).unwrap();
        for ligne in contenu.lines() {
            if ligne.starts_with('@') {
                assert!(ligne.ends_with(";size=3"),
                    "{nom} : en-tête doit finir par ;size=3, obtenu : {}", ligne);
            }
        }
        assert_eq!(compter_sequences_fastq(&contenu), 5, "{nom} : doit contenir 5 séquences uniques");
    }
}

#[test]
fn test_abundance_custom_separator() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq_repetitions(&temp_dir, "entree.fastq", 5, 3);
    let sortie = temp_dir.child("out.fq");
    cmd()
        .arg("-1").arg(entree.path())
        .arg("--abundance")
        .arg("--abundance-sep").arg(";count=")
        .arg("-o").arg(sortie.path())
        .assert().success();
    let contenu = fs::read_to_string(sortie.path()).unwrap();
    for ligne in contenu.lines() {
        if ligne.starts_with('@') {
            assert!(ligne.contains(";count="), "Doit contenir ;count=, obtenu : {}", ligne);
            assert!(!ligne.contains(";size="), "Ne doit PAS contenir ;size=, obtenu : {}", ligne);
        }
    }
    assert_eq!(compter_sequences_fastq(&contenu), 5);
}

#[test]
fn test_abundance_temp_dir_matches_classic() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let entree = generer_fastq_repetitions(&temp_dir, "entree.fastq", 10, 3);
    let sortie_classic = temp_dir.child("classic.fq");
    let sortie_temp = temp_dir.child("tempdir.fq");
    let temp_subdir = temp_dir.child("mon_tempdir");
    fs::create_dir_all(temp_subdir.path()).unwrap();

    cmd().arg("-1").arg(entree.path()).arg("--abundance")
        .arg("-o").arg(sortie_classic.path()).assert().success();
    cmd().arg("-1").arg(entree.path()).arg("--abundance")
        .arg("--temp-dir").arg(temp_subdir.path())
        .arg("-o").arg(sortie_temp.path()).assert().success();

    let c_classic = fs::read_to_string(sortie_classic.path()).unwrap();
    let c_temp = fs::read_to_string(sortie_temp.path()).unwrap();
    let mut h_classic: Vec<&str> = c_classic.lines().filter(|l| l.starts_with('@')).collect();
    let mut h_temp: Vec<&str> = c_temp.lines().filter(|l| l.starts_with('@')).collect();
    h_classic.sort_unstable();
    h_temp.sort_unstable();
    assert_eq!(h_classic, h_temp, "Les en-têtes triés doivent être identiques");
    assert_eq!(compter_sequences_fastq(&c_classic), compter_sequences_fastq(&c_temp));
}

#[test]
fn test_abundance_temp_dir_reports_duplicates() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    // Use the built-in generer_fastq: 70 uniques + 30 dups of first 30 indices = 30/100 = 30%
    let fichier = temp_dir.child("entree.fastq");
    let mut contenu = String::new();
    for i in 0..70 {
        contenu.push_str(&format!("@read_{} 1:N:0:ATGC\n{}\n+\n{}\n",
            i, generer_sequence(i, 0), generer_qualite(i)));
    }
    for i in 0..30 {
        contenu.push_str(&format!("@read_{}_dup 1:N:0:ATGC\n{}\n+\n{}\n",
            i, generer_sequence(i, 0), generer_qualite(i + 1000)));
    }
    fichier.write_str(&contenu).unwrap();
    let sortie = temp_dir.child("out.fq");
    cmd()
        .arg("-1").arg(fichier.path())
        .arg("--abundance").arg("--temp-dir").arg("/tmp").arg("-v")
        .arg("-o").arg(sortie.path())
        .assert().success()
        .stdout(predicate::str::contains("Duplicates removed: 30.00%"));
}

#[test]
fn test_pre_count_matches_estimate() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let fichier = temp_dir.child("entree.fastq");
    let mut contenu = String::new();
    for i in 0..10 {
        contenu.push_str(&format!("@read_{} 1:N:0:ATGC\n{}\n+\n{}\n",
            i, generer_sequence(i, 0), generer_qualite(i)));
    }
    // 2 dups
    for i in 0..2 {
        contenu.push_str(&format!("@read_{}_dup 1:N:0:ATGC\n{}\n+\n{}\n",
            i, generer_sequence(i, 0), generer_qualite(i + 1000)));
    }
    fichier.write_str(&contenu).unwrap();
    let sortie = temp_dir.child("out.fq");
    cmd()
        .arg("-1").arg(fichier.path())
        .arg("--pre-count").arg("-v")
        .arg("-o").arg(sortie.path())
        .assert().success()
        .stdout(predicate::str::contains("Found"));
    let contenu_sortie = fs::read_to_string(sortie.path()).unwrap();
    assert_eq!(compter_sequences_fastq(&contenu_sortie), 10);
}

#[test]
fn test_pe_abundance_detects_desync() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let r1 = temp_dir.child("R1.fastq");
    let r2 = temp_dir.child("R2.fastq");
    let mut cr1 = String::new();
    let mut cr2 = String::new();
    // Pair 1: synchronized
    cr1.push_str(&format!("@SEQ_001 1:N:0:ATGC\n{}\n+\n{}\n",
        generer_sequence(0, 0), generer_qualite(0)));
    cr2.push_str(&format!("@SEQ_001 2:N:0:ATGC\n{}\n+\n{}\n",
        generer_sequence(0, 9999), generer_qualite(500)));
    // Pair 2: DESYNCHRONIZED (SEQ_002 vs SEQ_999)
    cr1.push_str(&format!("@SEQ_002 1:N:0:ATGC\n{}\n+\n{}\n",
        generer_sequence(1, 0), generer_qualite(1)));
    cr2.push_str(&format!("@SEQ_999 2:N:0:ATGC\n{}\n+\n{}\n",
        generer_sequence(1, 9999), generer_qualite(1500)));
    r1.write_str(&cr1).unwrap();
    r2.write_str(&cr2).unwrap();

    let o1 = temp_dir.child("out_R1.fq");
    let o2 = temp_dir.child("out_R2.fq");
    cmd().arg("-1").arg(r1.path()).arg("-2").arg(r2.path())
        .arg("--abundance").arg("-o").arg(o1.path()).arg("-p").arg(o2.path())
        .assert().failure();

    let o3 = temp_dir.child("out2_R1.fq");
    let o4 = temp_dir.child("out2_R2.fq");
    cmd().arg("-1").arg(r1.path()).arg("-2").arg(r2.path())
        .arg("--abundance").arg("--temp-dir").arg("/tmp")
        .arg("-o").arg(o3.path()).arg("-p").arg(o4.path())
        .assert().failure();
}
