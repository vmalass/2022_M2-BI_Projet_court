# Projet court n°2 : Calcul de la surface accessible au solvant d’une protéine
## Description : 
L'objectif de ce programme python est de calculer la surface accessible d'une protéine au slovant (modélisé par une molécule d'eau)

## Prerequis :
### Systeme :
Noyaux Unix 

### Environnement :
Le fichier .yml contient l'environnement conda, voici les lignes de code à faire un terminal pour recréer l'environnement :
```bash
$ conda env create -f env_conda.yml
$ conda activate env_conda
```

** N'hésitez pas à utiliser mamba pour plus de rapidité :
```bash
$ conda install mamba -n base -c conda-forge
$ mamba env create -f env_conda.yml
$ conda activate env_conda
```

###  Utilisation :
1. Cloner le dépôt :
```bash
$ git clone https://github.com/vmalass/2022_M2-BI_Projet_court.git
```
vous pouvez également télécharger le fichier zip que le décompresser

2. Exécution du code :
```bash
$ python main.py file.pdb
```

### Procedure et exmple :
En résumé voici l'ordre d'exécution des commandes dans le terminal avec la protéine 3i40:
```bash
$ git clone https://github.com/vmalass/2022_M2-BI_Projet_court.git
$ conda env create -f env_conda.yml
$ conda activate env_conda
$ python main.py 3i40.pdb
```
Temps d'execution ~5s
Résultat obtenu :
```bash
Solvent surface protein accessible per atom : 3297.35 Å
Exposed surface per residue : 51.00 Å
Percentage of accessible surface : 5.94 %
```
Avec une autre protéine 6a5j:
```bash
$ python main.py 6a5j.pdb
```
Temps d'execution ~5s
Résultat obtenu :
```bash
Solvent surface protein accessible per atom : 1710.57 Å
Exposed surface per residue : 13.00 Å
Percentage of accessible surface : 5.46 %
```
Avec une autre protéine 7u52 de 135 résidues:
```bash
$ python main.py 7u52.pdb
```
Temps d'execution ~1min30s
Résultat obtenu :
```bash
Solvent surface protein accessible per atom : 33308.80 Å
Exposed surface per residue : 747.00 Å
Percentage of accessible surface : 4.01 %
```

## Auteurs :
Malassigné Victor