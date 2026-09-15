#!/bin/bash

# Couleurs pour un affichage lisible dans le terminal
VERT='\033[0;32m'
BLEU='\033[0;34m'
JAUNE='\033[1;33m'
NC='\033[0m' # Pas de couleur

echo -e "${BLEU}======================================================${NC}"
echo -e "${BLEU}           DÉMONSTRATION DE L'OUTIL FDEDUP            ${NC}"
echo -e "${BLEU}======================================================${NC}\n"

# 1. Compilation de l'outil
echo -e "${JAUNE}[1/4] Compilation du projet Rust en mode release...${NC}"
cargo build --release
echo -e "${VERT}Compilation terminée !${NC}\n"

# 2. Préparation des données de test
echo -e "${JAUNE}[2/4] Génération des données FASTQ de démonstration...${NC}"
mkdir -p demo_data

# ---> Fichier Single-End (3 reads: 2 identiques, 1 unique)
cat <<EOF > demo_data/single.fastq
@SEQ_1_ORIGINAL
ATGCATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_2_DUPLICAT
ATGCATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_3_UNIQUE
CGTACGTACGTACGTACGTA
+
FFFFFFFFFFFFFFFFFFFF
EOF

# ---> Fichiers Paired-End (3 paires: 2 identiques, 1 unique)
cat <<EOF > demo_data/R1.fastq
@SEQ_1 1:N:0:ATGC
ATGCATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_2 1:N:0:ATGC
ATGCATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_3 1:N:0:ATGC
CGTACGTACGTACGTACGTA
+
FFFFFFFFFFFFFFFFFFFF
EOF

cat <<EOF > demo_data/R2.fastq
@SEQ_1 2:N:0:ATGC
TTAATTAATTAATTAATTAA
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_2 2:N:0:ATGC
TTAATTAATTAATTAATTAA
+
FFFFFFFFFFFFFFFFFFFF
@SEQ_3 2:N:0:ATGC
GCCCGCCCGCCCGCCCGCCC
+
FFFFFFFFFFFFFFFFFFFF
EOF

echo -e "${VERT}Fichiers générés dans le dossier ./demo_data/${NC}\n"

# 3. Lancement du mode Single-End
echo -e "${JAUNE}[3/5] Test du Mode Single-End...${NC}"
echo "Attendu : 3 fragments traités, 1 doublon supprimé (33.33%)"
./target/release/fdedup \
    -1 demo_data/single.fastq \
    -o demo_data/single_dedup.fastq \
    -v
echo -e "\n"

# 4. Lancement du mode Paired-End
echo -e "${JAUNE}[4/5] Test du Mode Paired-End...${NC}"
echo "Attendu : 3 fragments traités, 1 doublon supprimé (33.33%)"
./target/release/fdedup \
    -1 demo_data/R1.fastq \
    -2 demo_data/R2.fastq \
    -o demo_data/R1_dedup.fastq \
    -p demo_data/R2_dedup.fastq \
    -v
echo -e "\n"

# 5. Démonstration du dry-run sans fichier de sortie
echo -e "${JAUNE}[5/5] Test du dry-run sans fichier de sortie (Single-End)...${NC}"
echo "Attendu : taux de duplication affiché, aucun fichier écrit"
./target/release/fdedup \
    -1 demo_data/single.fastq \
    -s -v
echo -e "\n"

echo -e "${VERT}======================================================${NC}"
echo -e "${VERT}                DÉMONSTRATION RÉUSSIE !               ${NC}"
echo -e "${VERT}======================================================${NC}"

rm -rf demo_data