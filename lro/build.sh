#!/bin/bash

set -eux
pdflatex --shell-escape test
bibtex test
pdflatex --shell-escape test
#dvipdf test.dvi test.pdf

pandoc -f latex -s -o test.pre.html test.tex --mathjax=https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS-MML_HTMLorMML --bibliography=test.bib

date_iso=$(date --iso-8601)

#cat > test.html <<EOF
#---
#title:  "Counting integer points in a convex rational polytope"
#categories: math algorithm
#excerpt: "...a rather informal presentation of the algorithm skimming over many proofs, focusing instead on the intuition."
#date: $date_iso
#---
#EOF
#
#cat test.pre.html >> test.html
#
#filename=barvinok-1
#PDFPATH=assets/files/$filename.pdf
#set -eu
#perl -0777 -i -pe 's#<h1 class="title">Counting integer points in a convex rational polytope</h1>\s+</div>#$&'"The latex version can be found <a href='/$PDFPATH'>here</a> (though Pandoc seems to have done its job rather well)#gs" test.html
#
#
#cp test.pdf /home/benoit/programmation/repo/jambonkapa.github.io/$PDFPATH
#cp test.html /home/benoit/programmation/repo/jambonkapa.github.io/_posts/${date_iso}-$filename.html


