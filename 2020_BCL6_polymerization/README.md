This repository holds scripts used for modeling and analysis in <b>Slabicki M<sup>\*</sup>, Yoon H<sup>\*</sup>, Koeppel J<sup>\*</sup>, Nitsch L, Roy Burman SS, Di Genua CA, Donovan KA, Sperling AS, Hunkeler M, Tsai JM, Sharma R, Guirguis A, Zou C, Chudasama P, Gasser JA, Miller PG, Scholl C, Fröhling S, Nowak RP, Fischer ES<sup>+</sup> & Ebert BL<sup>+</sup> (2020). “Small molecule-induced polymerization triggers degradation of BCL6.”</b>

<ol>
<li><code>tmt4degrader_basic.R</code> was used for analysis of TMT-based whole proteome mass spectrometry experiments using limma moderated t-statistics. The inputs required are infile, sc, pre, w, contrast, and tmt.</li>
<li><code>docking_script.sh</code> shows how to use the Rosetta local docking protocol to dock the BCL6 dimers with BI-3802 at the interface. The options used are provided in <code>docking_protocol_flags</code>. Prerequisites are the installation of Rosetta and its dependencies. The inputs required are a BCL6 tetramer starting structure and a params file describing BI-3802.</li>
<li><code>oligomerize_from_local_dock.py</code> was used to create create structures of BCL6 polymers with docking modes copied from the docked structures. Prerequisites are Python 3+ and PyMOL python library. Inputs required are the docked structures.
</ol>

If using any part of the scripts, please cite the paper above.