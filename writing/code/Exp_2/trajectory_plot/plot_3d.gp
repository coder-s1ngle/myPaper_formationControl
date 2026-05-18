# ============================================================
#  3D Trajectory Plot — gnuplot 6.0
#  Strictly matches matplotlib plot_3d_trajectory parameters
# ============================================================

data_dir   = '../output/exp_2_target_capture/gnuplot_data'
out_png    = '../output/exp_2_target_capture/trajectory_3d_gnuplot.png'

load data_dir . '/ranges.cfg'

# ---- terminal (wider to accommodate Y-dominant box_aspect) ---
set terminal pngcairo size 2800,2000 enhanced font 'Arial,10'
set output out_png

# ---- view: match matplotlib elev=18 azim=28 -----------------
set view 72, 28
set view equal xy

# ---- margins ------------------------------------------------
set lmargin 5
set rmargin 4
set tmargin 3
set bmargin 8

# ---- axis ranges --------------------------------------------
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set zrange [zmin:zmax]

# ---- ticks: exact positions from matplotlib MaxNLocator ------
set xtics @xtics_vals font ',8'
set ytics @ytics_vals font ',8'
set ztics @ztics_vals font ',8'
set tics in scale 1.5

# ---- axis labels --------------------------------------------
set xlabel 'X / m' font ',10'
set ylabel 'Y / m' font ',10'
set zlabel 'Z / m' font ',10' offset 1,-1

# ---- clean style --------------------------------------------
unset grid
unset colorbox

# ---- border: thin light-gray ---------------------------------
set border lc '#cccccc' lw 0.8

# ---- key (legend): horizontal, compact, below ----------------
set key horizontal font ',7' spacing 0.8 samplen 2.0 width -4
set key at screen 0.5, screen 0.05 center

# ---- xy plane at data bottom --------------------------------
set xyplane at zmin

# ==============================================================
#   Plot
# ==============================================================
splot \
    data_dir . '/traj_1.dat' u 2:3:4 w l lc '#2166AC' lw 0.8 dt 2 t 'UAV 1', \
    data_dir . '/traj_2.dat' u 2:3:4 w l lc '#D6604D' lw 0.8 dt 2 t 'UAV 2', \
    data_dir . '/traj_3.dat' u 2:3:4 w l lc '#4DAF4A' lw 0.8 dt 2 t 'UAV 3', \
    data_dir . '/traj_4.dat' u 2:3:4 w l lc '#9467BD' lw 0.8 dt 2 t 'UAV 4', \
    data_dir . '/snap_edges.dat'  u 1:2:3 w l lc '#F60606' lw 1.5 notitle, \
    data_dir . '/snap_labels.dat' u 1:2:3:4 w labels font ',8' tc '#333333' notitle

set output
