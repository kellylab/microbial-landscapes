#!/bin/bash
ffmpeg -framerate 2 -i frames/B%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		B.mp4
ffmpeg -framerate 2 -i frames/A%d.png \
		-start_number 1 \
		-pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
		A.mp4
