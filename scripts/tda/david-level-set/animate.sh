ffmpeg -framerate 2 -i frames/%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		level-set-animation.mp4
