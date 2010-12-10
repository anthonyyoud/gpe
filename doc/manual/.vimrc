au BufWritePre *.tex :silent !sed -i -r "s/[A-Za-z]{3} [0-9]{2} [A-Za-z]{3} [0-9]{2}:[0-9]{2}:[0-9]{2} 20[0-9]{2}/`date '+\%a \%d \%b \%T \%Y'`/g" manual.tex
