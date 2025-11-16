# 💀 lateNDthal_demography

Data and Code for: “Archaeogenetic insights into the demographic history of Late Neanderthals”  

### Authors: 
Charoula M. Fotiadou¹,²,  
Jesper Borre Pedersen³,  
Hélène Rougier⁴,  
Mirjana Roksandic⁵,  
Maria A. Spyrou¹,²,  
Kathrin Nägele⁶,  
Ella Reiter¹,  
Hervé Bocherens²,⁷,  
Andrew W. Kandel³,  
Miriam N. Haidle³,⁸,  
Timo P. Streicher³,  
Nicholas J. Conard²,⁸,  
Flora Schilt⁹,¹⁰,  
Ricardo Miguel Godinho¹⁰,  
Thorsten Uthmeier¹¹,  
Luc Doyon¹²,  
Patrick Semal¹³,  
Johannes Krause¹,⁶,  
Alvise Barbieri¹⁰,  
Dušan Mihailović¹⁴,  
Isabelle Crevecoeur¹²,  
Cosimo Posth¹,² 

### Affiliation: 
1. Archaeo- and Paleogenetics, Institute for Archaeological Sciences, Department of Geosciences, University of Tübingen, Tübingen 72074, Germany;
2. Senckenberg Centre for Human Evolution and Palaeoenvironment at the University of Tübingen, Tübingen 72074, Germany;
3. The Role of Culture in Early Expansions of Humans (ROCEEH), Heidelberg Academy of Sciences and Humanities, University of Tübingen, Hölderlinstrasse 12, Tübingen 72074, Germany;
4. Department of Anthropology, California State University Northridge, Northridge, CA 91330, USA;
5. Department of Anthropology, University of Winnipeg, Winnipeg, MB R3T 3C7, Canada;
6. Department of Archaeogenetics, Max Planck Institute for Evolutionary Anthropology, 04103 Leipzig, Germany;
7. Biogeology, Department of Geosciences, University of Tübingen, 72074 Tübingen, Germany;
8. Department of Early Prehistory and Quaternary Ecology, University of Tübingen, 72070 Tübingen, Germany;
9. Department of Art and Culture, History and Antiquity, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands;
10. Interdisciplinary Center for Archaeology and the Evolution of Human Behavior, University of Algarve, Faro, Portugal;
11. Department of Classical World and Asian Cultures, Institute of Prehistory and Protohistory, Friedrich-Alexander Universität Erlangen–Nürnberg, 91054 Erlangen, Germany;
12. PACEA UMR 5199, CNRS, Université de Bordeaux, Ministère de la Culture, Pessac, France;
13. Service of Scientific Heritage, Royal Belgian Institute of Natural Sciences, 1000 Brussels, Belgium;
14. Department of Archaeology, Faculty of Philosophy, University of Belgrade, 11000 Belgrade, Serbia;

### Date: YYYY-MM-DD

### Publication DOI: *[Add DOI or URL here]*

### Compendium DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16737433.svg)](https://doi.org/10.5281/zenodo.16737433)

---

## 📖 Content

- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Data Description](#-data-description)
- [Getting started with the code](#-getting-started-with-the-code)
- [Software Requirements](#-software-requirements)
- [License](#-license)
- [Acknowledgements (from paper)](#-acknowledgements-(from-paper))
- [Contact](#-contact)


---

## 🧭 Overview

This repository contains the datasets, R scripts, and documentation for the analyses presented in the academic paper titled:

**“Archaeogenetic insights into the demographic history of Late Neanderthals”**

---

## 📁 Repository Structure

```
├── 1_data
│   ├── dating
│   ├── map
│   ├── pairwise_distance
│   ├── road_analysis
│   ├── yaworsky_et_al_2024
│   └── yaworsky_extended
├── 2_scripts
│   ├── dating
│   ├── map
│   ├── pairwise_distance
│   ├── road_analysis
│   └── yaworsky_extended
├── 3_output
│   ├── dating
│   ├── map
│   ├── pairwise_distance
│   ├── road_analysis
│   └── yaworsky_extended
├── lateNDthal_demography.Rproj
└── README.md
```


> ⚠️ The `1_data/yaworsky_extended/raw_data/Climate/` folder is not tracked by Git due to its large size. This folder is created when running the code.

---

## 📊 Data Description
- `1_data/dating/` – All data used for creating Fig. 3B and Fig. S18
- `1_data/map/` – All data used for creating Fig. 1
- `1_data/pairwise_distance/` – All data used for checking the pairwise distance (Fig. 2B)
- `1_data/road_analysis/` – All data used creating distribution maps.
- `1_data/yaworsky_extended/` – All data used for running the same analysis as Yaworsky et al. 2024 with an extended dataset. These data were all downloaded using the PHP scripts provided in the markdown document accompanying the original paper by Yaworsky et al. (2024). They were retrieved using the same copy-paste method described in the paper, except for the ROCEEH Neanderthal Data, which were too large for this approach and were therefore downloaded using the *0.1.load_prep_save_neanderthal_data.R* script.
- `1_data/yaworsky_et_al_2024/` – The orignial data from Yaworsky et al. 2024 for comparison.

---
## ⚙️ Getting started with the code

---

## 💾 Software Requirements

---

## 📜 License

This project is licensed under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to:
- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- **Attribution** — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

For full legal terms, see the [CC-BY 4.0 License](https://creativecommons.org/licenses/by/4.0/legalcode).

---

## 🤝 Acknowledgements (from paper)

We thank the Archeo- and Paleogenetics team at the University of Tübingen for comments, P. Yaworsky for help in reproducing a published analysis with an expanded dataset, S. Anastasios for work on an earlier version of the Pešturina 3 mtDNA, and C. Schwab (Musée d’Archéologie nationale) for access to the Saint-Césaire collection in her care. MR and DM are funded by NSERC grant (no RGPAS-2019-00039), by SSHRC Partnership grant (no 895-2024-1005), and NEEMO grant (no 7746827) funded by the Science Fund of the Republic of Serbia. AB is funded by the Portuguese Ministry of Science (2002.08622.CEECIND) and has received funding for the analysis of the Sesselfelsgrotte individual by the National Geographic Society (NGS-96087R-22). JBP, AWK and MNH are funded through the research center ROCEEH of the Heidelberg Academy of Sciences and Humanities (https://www.hadw-bw.de/) which is promoted by the Joint Science Conference of the Federal Government and the state governments of the Federal Republic of Germany in the Academies‘ Programme of the Union of the German Academies (https://www.akademienunion.de/forschung/akademienprogramm/). HR received support from the CSUN Competition for RSCA Awards. The Collective Research Project at La Roche-à-Pierrot (IC, dir.) is funded by the Direction Régionale des Affaires Culturelles (DRAC) of the Région Nouvelle-Aquitaine, by the Département de Charente-Maritime (CG 17, France) and by the University of Bordeaux's IdEx “Investments for the Future” program / GPR “Human Past. Research at the Tourtoirac rock shelter (LD, dir.) benefited from the financial support from the Direction des Affaires Culturelles - Nouvelle-Aquitaine, the Shandong University 111 Project (no 111-2-20), the University of Bordeaux via its IdEx “Talent” (191022-001), “Bordeaux International Support” (no 191203-003) and “Investments for the Future / GPR Human Past” programs as well as the European Research Council Starting Grant for the ExOsTech project (no 101161065). RMG is funded by Fundação para a Ciência e a Tecnologia (FCT; contract reference 2020.00499.CEECIND; https://doi.org/10.54499/2020.00499.CEECIND/CP1613/CT0002 and by the FCT R&D research project “ParaFunction” (project reference 2022.07737.PTDC; https://doi.org/10.54499/2022.07737.PTDC).

---

## 📬 Contact

For questions or data requests, please contact:

**Jesper Borre Pedersen**  
Email: <jesper-borre.pedersen@ifu.uni-tuebingen.de>  
GitHub: https://github.com/JesperBorrePedersen  
ORCID: [![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--3468--0986-green.svg)](https://orcid.org/0000-0002-3468-0986)

**Charoula M. Fotiadou**

Email: <charoula.fotiadou@zv.uni-tuebingen.de>

GitHub: https://github.com/charoulafotiadou

