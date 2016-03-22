(TeX-add-style-hook "bmvc_review"
 (lambda ()
    (LaTeX-add-bibliographies
     "egbib")
    (LaTeX-add-labels
     "fig:pipeline"
     "fig:5ptm"
     "eq:central_projection"
     "eq:central_projection1"
     "eq:intrinsic"
     "eq:point_motion"
     "eq:h_inf")
    (TeX-add-symbols
     "eg"
     "Eg"
     "etal")
    (TeX-run-style-hooks
     "wrapfig"
     "smartdiagram"
     "etoolbox"
     "color"
     "amssymb"
     "mathtools"
     "graphicx"
     "latex2e"
     "bmvc2k10"
     "bmvc2k"
     ""
     "r_err"
     "t_err")))

