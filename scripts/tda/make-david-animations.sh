ffmpeg -framerate 2 -i david-subject-trajectories/frames/B%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		david-subject-trajectories/B.mp4
