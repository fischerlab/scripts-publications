This repository holds scripts used for modeling and analysis in Slabicki M*, Yoon H*, Koeppel J*, Nitsch L, Roy Burman SS, Di Genua CA, Donovan KA, Sperling AS, Hunkeler M, Tsai JM, Sharma R, Guirguis A, Zou C, Chudasama P, Gasser JA, Miller PG, Scholl C, Fröhling S, Nowak RP, Fischer ES & Ebert BL (2020). “Small molecule-induced polymerization triggers degradation of BCL6.”

<ol>
<li>`tmt4degrader_basic.R` was used for analysis of TMT-based whole proteome mass spectrometry experiments using limma moderated t-statistics. The inputs required are infile, sc, pre, w, contrast, and tmt.</li>
<li>`docking_script.sh` shows how to use the Rosetta local docking protocol to dock the BCL6 dimers with BI-3802 at the interface. The options used are provided in `docking_protocol_flags`. Prerequisites are the installation of Rosetta and its dependencies. The inputs required are a BCL6 tetramer starting structure and a params file describing BI-3802.</li>
<li>`oligomerize_from_local_dock.py` was used to create create structures of BCL6 polymers with docking modes copied from the docked structures. Prerequisites are Python 3+ and PyMOL python library. Inputs required are the docked structures.
</ol>

If using any part of the scripts, please cite the paper above.