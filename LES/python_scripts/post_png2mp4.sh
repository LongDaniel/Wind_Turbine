# see http://trac.ffmpeg.org/wiki/Slideshow

pre=""

rm -r ffmpeg_export
mkdir -p ffmpeg_export
for i in $(seq -f "%05g" 1 1 200)
do
  n1=$(echo $i/1 | bc)
  echo $n1
  cp export/export_slice_y_0.5_iso_q_2.5_ti_$i.png ffmpeg_export/temp-$n1.png
done

cd ffmpeg_export

ffmpeg -start_number 1 -i temp-%d.png -c:v libx264 -strict -2 -preset slow -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -r 10 -f mp4 output.mp4

cd ..

##!/bin/bash
## Following lines are used for GIF generation.

#DELAY=0.372
#LOOP=0
#DISPOSE=Background

#INPUT=export_yslice_u_q_iso_30_*.png
#OUTPUT=output_3.gif

#mogrify -format gif ./ffmpeg_export/*.png
#gifsicle --delay=20  ./ffmpeg_export/*.gif > anim.gif 

#cd $(dirname $0)

#convert -delay $DELAY -loop $LOOP -alpha set -dispose $DISPOSE $INPUT $OUTPUT

# If you prefer a transparent background,

# use `-transparent` option.
