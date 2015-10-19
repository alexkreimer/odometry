(TeX-add-style-hook "survey"
 (lambda ()
    (LaTeX-add-bibliographies)
    (TeX-run-style-hooks
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

