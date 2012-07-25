(TeX-add-style-hook "main"
 (lambda ()
    (TeX-add-symbols
     '("footlineextra" 1)
     '("software" 1)
     '("Rcode" 1)
     '("Rclass" 1)
     '("Rfunarg" 1)
     '("Rmethod" 1)
     '("Rpackage" 1)
     '("Robject" 1)
     '("Rfunction" 1)
     "R"
     "Bioc"
     "IRanges"
     "biovizBase"
     "ggbio"
     "visnab"
     "ggplot"
     "gridExtra"
     "qplot"
     "autoplot"
     "gr"
     "newblock"
     "frame")
    (TeX-run-style-hooks
     "verbatim"
     "natbib"
     "longtable"
     "overpic"
     "Sweave"
     "hyperref"
     "graphicx"
     "latex2e"
     "beamer10"
     "beamer")))

