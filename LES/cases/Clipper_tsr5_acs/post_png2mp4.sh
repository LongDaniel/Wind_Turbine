# see http://trac.ffmpeg.org/wiki/Slideshow
module load ffmpeg
source ~/env_utilities.source

rm -r ffmpeg_export
mkdir -p ffmpeg_export
for i in $(seq -f "%05g" 1 1 50)
do
  n1=$(echo $i/1 | bc)
  echo $n1
  #mv export/export_slice_y_0.5_iso_q_2.5_ti_$n1.png ffmpeg_export/export_slice_y_0.5_iso_q_2.5_ti_$i.png
  cp export/export_slice_y_0.5_iso_q_2.5_ti_$i.png ffmpeg_export/temp-$n1.png
  #cp export_q_iso_30_ti_$i.png ffmpeg_export/export_q_iso_30_$n1.png
  #cp export_u_q_iso_30_ti_$i.png ffmpeg_export/export_u_q_iso_30_$n1.png
  #cp export_yslice_u_q_iso_30_ti_$i.png ffmpeg_export/export_yslice_u_q_iso_30_$n1.png
done

cd ffmpeg_export

ffmpeg -framerate 5.0 -qscale:v 0 -i temp-%d.png output_1.avi
#ffmpeg -framerate 2.69 -i ./ffmpeg_export/export_q_iso_30_%d.png ./ffmpeg_export/output_1.mp4
#ffmpeg -framerate 2.69 -i ./ffmpeg_export/export_u_q_iso_30_%d.png ./ffmpeg_export/output_2.mp4
#ffmpeg -framerate 2.69 -i ./ffmpeg_export/export_yslice_u_q_iso_30_%d.png ./ffmpeg_export/output_3.mp4

cd ..

#!/bin/bash
DELAY=0.372
LOOP=0
DISPOSE=Background

INPUT=export_yslice_u_q_iso_30_*.png
OUTPUT=output_3.gif

mogrify -format gif ffmpeg_export/*.png
gifsicle --delay=20  ffmpeg_export/*.gif > anim.gif 

#cd $(dirname $0)

#convert -delay $DELAY -loop $LOOP -alpha set -dispose $DISPOSE $INPUT $OUTPUT

# If you prefer a transparent background,

# use `-transparent` option.
