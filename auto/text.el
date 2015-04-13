(TeX-add-style-hook
 "text"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "10pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "amsfonts"
    "amssymb"
    "mathtools"
    "graphicx"
    "caption"
    "subcaption"
    "tikz")
   (LaTeX-add-labels
    "eq:prob"
    "eq:sparsej"
    "sec: rot_mat")
   (LaTeX-add-environments
    "exmp")))

