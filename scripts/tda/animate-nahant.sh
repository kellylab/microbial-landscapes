#!/bin/bash
ffmpeg -framerate 2 -i nahant-trajectories/frames/Bacteria%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		nahant-trajectories/Bacteria.mp4
ffmpeg -framerate 2 -i nahant-trajectories/frames/Eukaryota%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		nahant-trajectories/Eukaryota.mp4
