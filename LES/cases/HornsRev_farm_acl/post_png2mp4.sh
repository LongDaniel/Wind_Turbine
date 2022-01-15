# see http://trac.ffmpeg.org/wiki/Slideshow
#module load ffmpeg
#source ~/env_ffmpeg.source
#source ~/env_utilities.source
rm -r ffmpeg_export
mkdir -p ffmpeg_export
for i in $(seq -f "%05g" 1 1 100)
do
  n1=$(echo $i/1 | bc)
  echo $n1
  #cp export/hr_style1_xy_u_$i.png ffmpeg_export/temp-$n1.png
  cp export/hr_style3_turbine_$i.png ffmpeg_export/temp-$n1.png
  #cp export/hr_style2_bottom_$i.png ffmpeg_export/temp-$n1.png
  #cp export/hr_style3_q_$i.png ffmpeg_export/temp-$n1.png
  #cp export/hr_style4_q_$i.png ffmpeg_export/temp-$n1.png
  #cp export/hr_style5_xz_u_$i.png ffmpeg_export/temp-$n1.png
done

cd ffmpeg_export

#ffmpeg -framerate 5.0 -qscale:v 0 -i temp-%d.png output_1.avi
ffmpeg -start_number 1 -i temp-%d.png -c:v libx264 -strict -2 -preset slow -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -f mp4 output.mp4
#ffmpeg  -start_number 1 -i temp-%d.png -c:v libx264 -crf 22 -f mp4 output.mp4
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

#mogrify -format gif ./ffmpeg_export/*.png
#gifsicle --delay=20  ./ffmpeg_export/*.gif > anim.gif 

#cd $(dirname $0)

#convert -delay $DELAY -loop $LOOP -alpha set -dispose $DISPOSE $INPUT $OUTPUT

# If you prefer a transparent background,

# use `-transparent` option.
