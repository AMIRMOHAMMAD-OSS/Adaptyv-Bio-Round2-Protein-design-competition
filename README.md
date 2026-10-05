# EGF Binders: Round 2 of Adaptyv Bio's Protein Design Competition

**AmirMohammad (@a.m.m.h2003)**

<p align="center">
  <a href="https://www.linkedin.com/in/amir-mohammad-mohammad-hosseini-9027b7229" target="_blank">
    <img alt="LinkedIn" src="https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white">
  </a>
  <a href="https://github.com/AMIRMOHAMMAD-OSS" target="_blank">
    <img alt="GitHub" src="https://img.shields.io/badge/GitHub-AMIRMOHAMMAD--OSS-181717?style=for-the-badge&logo=github&logoColor=white">
  </a>
  <a href="mailto:a.m.m.h2003@gmail.com">
    <img alt="Email" src="https://img.shields.io/badge/Email-a.m.m.h2003%40gmail.com-EA4335?style=for-the-badge&logo=gmail&logoColor=white">
  </a>
</p>

Adaptyv Bio accepted submissions for Round 2 of their Protein Design Competition.  
I designed *20* binders targeting **Epidermal Growth Factor (EGF)**; **6** were selected after in‑silico filters and proceeded to **wet‑lab evaluation**.

<p align="center">
  <img src="Picture1.svg" alt="Kinetic curves (BLI) for representative designs" width="720px"/>
</p>

<p align="center">
  <img src="Screenshot 2025-10-10 215832.png" alt="Kinetic curves (BLI) for representative designs" width="720px"/>
</p>

## Methods (brief)
1. **Target prep** – EGF structure and interface context with EGFR DIII for hotspot guidance.
2. **Sequence generation** – multiple independent seeds; length and secondary‑structure constraints varied.
3. **In‑silico filters** – structural confidence and interface quality; basic developability screens (e.g., charge, hydrophobics, motifs).
4. **Complex modeling** – binder:EGF models; interface sanity checks and contact maps.
5. **Down‑selection** – top 6 forwarded for wet‑lab.
6. **Wet‑lab** – BLI runs at several concentrations; global fitting of association/dissociation to estimate k_on, k_off, and K_D.


## Methods: template‑guided design from an EGFR:TGF‑α complex

I started from the **EGFR–TGF‑α** co‑crystal structure (TGF‑α is a natural ligand of EGFR). I cropped the complex to the **binding region** (positions **30–100** in my design numbering) and identified **anchor residues** at the interface using a distance cutoff and hydrogen‑bond criteria. Those anchor residues were **fixed**, the remaining positions were **masked**, and designs were generated in two stages:

1. **Backbone/shape generation** with **RFdiffusion** using the fixed‑residue constraints to preserve the native contact geometry at the interface.
2. **Sequence design** with **ProteinMPNN**, seeded on the RFdiffusion backbones while keeping the fixed residues immutable.

The final models were re‑scored with standard structure/interaction filters and down‑selected for wet‑lab testing (BLI).

## Reference to the competition analysis

> Cotet, T.-S.; Krawczuk, I.; Stocco, F.; Ferruz, N.; Gitter, A.; Kurumida, Y.; de Almeida Machado, L.; Paesani, F.; Calia, C. N.; Challacombe, C. A.; Haas, N.; Qamar, A.; Correia, B. E.; Pacesa, M.; Nickel, L.; Subr, K.; Castorina, L. V.; Campbell, M. J.; Ferragu, C.; Kidger, P.; Hallee, L.; Wood, C. W.; Stam, M. J.; Kluonis, T.; Ünal, S. M.; Belot, E.; Naka, A.; Adaptyv Competition Organizers. **Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition.** *bioRxiv* (2025). https://doi.org/10.1101/2025.04.17.648362

## Runnable Ubuntu pipeline

The [pipeline](pipeline/README.md) contains the Ubuntu setup, resumable model downloads, RFdiffusion backbone generation, ProteinMPNN sequence design, ColabFold complex prediction, scoring, selection, and optional BLI fitting.

The supplied example targets **EGFR** using a TGF-alpha motif from PDB 1MOX. It is a new implementation of the method description; the original competition configuration and exact submitted designs cannot be recovered from the report alone. CPU workflow tests pass, while full GPU model inference remains to be verified.
