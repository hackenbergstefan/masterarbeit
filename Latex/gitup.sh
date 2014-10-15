#!/bin/sh
#cp ./latexmk_files/document.pdf Masterarbeit.pdf
if [ -z "$1" ]; then
  git commit -a -m "tex weiter"
else
  git commit -a -m "$1"
fi
git push
