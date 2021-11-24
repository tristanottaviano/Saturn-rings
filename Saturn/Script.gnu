# Petit script a lancer avec gnuplot: >load 'script.gnu'
# pour tracer sur une meme figure
# - la carte de couleurs de la fonction f(x,y) dont les données sont enregistrées dans le fichiers 'vagues.txt'
# - les contours associes
# - la trajectoire de maximisation enregistree dans le fichiers 'grad.txt'

# Mode 3D, carte de couleurs
set pm3d map

# Pour avoir la carte et les contours
#set contour surface
# Les niveaux de contours
#set cntrparam levels auto 50

# La palette de couleurs
set palette rgbformulae 33,13,10

# Les axes
set size ratio 1
set xlabel 'Time(s)'
set ylabel 'Radius(m)'
set title 'Asteroid density as a function of radius and time'
set cbrange[0:70]

# Pas de legende
unset key

# Tracer la carte uniquement
splot 'Hist_Formation.txt' using 1:2:3
