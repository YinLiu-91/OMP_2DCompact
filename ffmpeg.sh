ffmpeg -framerate 16 -pattern_type glob -i 'P.*.png' -r 30 -b:v 60000k -vcodec mpeg4 ffmpeg_P.mp4
