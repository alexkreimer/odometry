(TeX-add-style-hook "text"
 (lambda ()
    (LaTeX-add-environments
     "exmp")
    (LaTeX-add-labels
     "eq:ml"
     "fig:poses"
     "fig:proj"
     "fig:graphs")
    (TeX-run-style-hooks
     "tikz"
     "subcaption"
     "caption"
     "graphicx"
     "mathtools"
     "amssymb"
     "amsfonts"
     "amsmath"
     ""
     "latex2e"
     "art10"
     "article"
     "10pt")))

