unset border
set xtics nomirror
set ytics nomirror
unset xtics
unset ytics

set param
set samples 1000
RR=50

### Parametros ajustaveis para visualizar (range de x e y)
YRANGE=130
XX=120
XMIN = XX-YRANGE/2.0
XMAX = XX+YRANGE/2.0
set xrange [XMIN:XMAX]
set yrange [0:YRANGE]
set size square

FIRST=int((XX-1.3*RR)/(w+b))
LAST=int((XX+1.3*RR)/(w+b))
TOP=BASE+h


X0=120.004  ## usado apenas para visualizar as SIMULACOES. 


## isto nao precisa mudar. é a parte que vai desenhar os pilares
###SUBSTRATE###
i=1
do for [n=FIRST:LAST] { 

x1=n*(w+b)
x2=n*(w+b)+w-1
x3=(n+1)*(w+b)

if (n != FIRST){
set arrow i nohead lc rgb '#808080' lw 4 from x1,TOP  to x2,TOP;  i=i+1 
set arrow i nohead lc rgb '#808080' lw 4 from x2,BASE to x2,TOP;  i=i+1 }
set arrow i nohead lc rgb '#808080' lw 4 from x2,BASE to x3,BASE; i=i+1
if (n != LAST){
set arrow i nohead lc rgb '#808080' lw 4 from x3,BASE to x3,TOP;    i=i+1 }

}
###SUBSTRATE###

## apenas para simu

file="file.dsf_0"
plot file u 2:($5==1? (int($2) % int(w+b) < w ? ($4<=TOP? 1/0:$4) : $4):1/0) ev 50 w p pt 7 ps 0.5 lc rgb "#1D84FF" t '',\
file u 2:($5==2? (int($2) % int(w+b) < w ? ($4<=TOP? 1/0:$4) : $4):1/0) ev 50 w p pt 7 ps 0.5 lc rgb "orange" t ''

## fim

do for [n=1:i-1] {
unset arrow n
}



