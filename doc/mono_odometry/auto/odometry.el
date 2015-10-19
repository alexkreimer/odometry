(TeX-add-style-hook "odometry"
 (lambda ()
    (LaTeX-add-bibliographies
     "sample")
    (LaTeX-add-environments
     "exmp")
    (LaTeX-add-labels
     "sec:mono_odo"
     "fig:two_views"
     "sec:cross_ratio"
     "eq:puret"
     "eq:cr"
     "fig:cross_ratio"
     "algo"
     "fig:stereo_rig")
    (TeX-run-style-hooks
     "cite"
     "tikz"
     "subcaption"
     "caption"
     "margin=1cm"
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

