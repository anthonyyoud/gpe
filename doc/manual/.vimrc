au BufWritePre *.tex :silent !sed -i "$ s/\(\%\%\) .* \(\%\%\)/\1 `date` \2/" manual.tex
